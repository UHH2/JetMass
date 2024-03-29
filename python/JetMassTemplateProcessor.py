#!/usr/bin/env pythonJMS.sh
import awkward as ak
import numpy as np
from coffea import processor
from coffea.nanoevents import BaseSchema
from coffea.analysis_tools import PackedSelection
import coffea.lookup_tools
import correctionlib
import hist
from coffea.util import save
import os
import glob
from coffea_util import CoffeaWorkflow
from utils import jms_correction_files
from copy import deepcopy

jetmass_path = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass"
ddtmaps_n2_path = f"{jetmass_path}/Histograms/ddtmaps_n2.npy"
ddtmaps_particlenet_path = f"{jetmass_path}/Histograms/ddtmaps_particlenet.npy"
kfactor_path = f"{jetmass_path}/NLOweights/"


class JMSTemplates(processor.ProcessorABC):
    def __init__(
            self,
            year: str = "2017",
            jec: str = "nominal",
            variation_weight: str = "nominal",
            trigger_sf_var: str = "nominal",
            tagger: str = "substructure",
            gen_sort: str = "pt",
    ):
        self._year = year
        self._jec = jec
        self._variation_weight = variation_weight
        self._trigger_sf_variation = trigger_sf_var
        self._tagger_approach = tagger
        self._gen_sort = gen_sort
        tagger_approaches = ["substructure", "particlenet", "particlenetDDT"]
        if self._tagger_approach not in tagger_approaches:
            raise NotImplementedError(
                "You chose a tagger ({}) approach that is not implemented. Choose among:".format(self._tagger_approach),
                tagger_approaches
            )

        dataset_ax = hist.axis.StrCategory([], name="dataset", growth=True)
        # fakes_ax = hist.axis.StrCategory([], name="fakes", growth=True)
        fakes_ax = hist.axis.Boolean(name="fakes")
        # shift_ax = hist.axis.StrCategory([], name="shift", growth=True)
        jec_applied_ax = hist.axis.StrCategory(
            [], name="jecAppliedOn", label="JEC applied on", growth=True
        )
        # jec_ax = hist.axis.StrCategory(
        #     ["raw", "pt", "pt_up", "pt_down","pt_","pt_mJ_up", "down"], name="JEC", label="JEC"
        # )
        mJ_ax = hist.axis.Regular(50, 0.0, 500.0, name="mJ", label=r"$m_{SD}$ [GeV]")
        pT_ax = hist.axis.Regular(300, 0.0, 3000.0, name="pt", label=r"$p_{T}$ [GeV]")
        eta_ax = hist.axis.Regular(100, -6.5, 6.5, name="eta", label=r"$\eta$")
        eta_regions_ax = hist.axis.Variable([0, 1.3, 2.5], name="abs_eta_regions", label=r"$|\eta|$")
        phi_ax = hist.axis.Regular(100, -4, 4, name="phi", label=r"$\Phi$")

        chf_ax = hist.axis.Regular(51, 0, 1.02, name="chf", label="CHF")
        nhf_ax = hist.axis.Regular(51, 0, 1.02, name="nhf", label="NHF")

        mJ_fit_ax = hist.axis.Regular(
            500, 0.0, 500.0, name="mJ", label=r"$m_{SD}$ [GeV]"
        )

        mPnet_fit_ax = hist.axis.Regular(
            500, 0.0, 500.0, name="mPnet", label=r"$m_{\mathrm{ParticleNet}}$ [GeV]"
        )

        rho_ax = hist.axis.Regular(100, -10.0, 0, name="rho", label=r"$\rho$")

        self._unfolding_ax = {
            "vjets": {
                "ptgen": hist.axis.Variable(
                    np.array([0, 650, 800, 1200, np.inf]),
                    name="ptgen",
                    label=r"$p_{T,\mathrm{gen}}$ [GeV]",
                ),
                "mJgen": hist.axis.Variable(
                    np.array([0.0, 70., 80, 90, np.inf]),
                    name="mJgen",
                    label=r"$m_{SD,\mathrm{gen}}$ [GeV]",
                ),
                "ptreco": hist.axis.Variable(
                    np.array([500, 575, 650, 725, 800, 1000, 1200, np.inf]),
                    name="ptreco",
                    label=r"$p_{T,\mathrm{reco}}$ [GeV]",
                ),
                "mJreco": hist.axis.Regular(
                    500,
                    0.0,
                    500.0,
                    name="mJreco",
                    label=r"$m_{SD,\mathrm{reco}}$ [GeV]",
                ),
            },
            "ttbar": {
                "ptgen": hist.axis.Variable(
                    np.array([0, 250, 400, 650, np.inf]),
                    name="ptgen",
                    label=r"$p_{T,\mathrm{gen}}$ [GeV]",
                ),
                "mJgen": hist.axis.Variable(
                    np.array(
                        [0, 55, 80, 87.5,  np.inf]
                    ),
                    name="mJgen",
                    label=r"$m_{SD,\mathrm{gen}}$ [GeV]",
                ),
                "ptreco": hist.axis.Variable(
                    np.array([200, 300, 400, 500, 650, np.inf]),
                    name="ptreco",
                    label=r"$p_{T,\mathrm{reco}}$ [GeV]",
                ),
                "mJreco": hist.axis.Regular(
                    500, 0, 500, name="mJreco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"
                ),
            },
        }

        self._pT_fit_ax = {
            "vjets": hist.axis.Variable(
                np.array([500, 650, 800, 1200, np.inf]),
                name="pt",
                label=r"$p_{T}$ [GeV]",
            ),
            "ttbar": hist.axis.Variable(
                np.array([200, 300, 400, 500, 650, np.inf]),
                name="pt",
                label=r"$p_{T}$ [GeV]",
            ),
        }

        hists = {}

        # create dense_lookup from custom n2ddt map
        corrected_str = {
            "none": "",
            "pt": "_corrected_pt",
            "pt&mJ": "_corrected_pt_mass",
        }

        n2ddtmap = np.load(ddtmaps_n2_path, allow_pickle=True).item()
        self._n2ddtmaps = {
            jec_applied_on: coffea.lookup_tools.dense_lookup.dense_lookup(
                n2ddtmap[
                    f"discddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                ][0],
                dims=(
                    n2ddtmap[
                        f"discddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                    ][1],
                    n2ddtmap[
                        f"discddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                    ][2],
                ),
            )
            for jec_applied_on in ["none", "pt", "pt&mJ"]
        }
        particlenetddtmap = np.load(ddtmaps_particlenet_path, allow_pickle=True).item()
        self._pNetMDWvsQCDddtmaps = {
            jec_applied_on: coffea.lookup_tools.dense_lookup.dense_lookup(
                particlenetddtmap[
                    f"discddt_map_{year}_smooth_4_0p975{corrected_str[jec_applied_on]}"
                ][0],
                dims=(
                    particlenetddtmap[
                        f"discddt_map_{year}_smooth_4_0p975{corrected_str[jec_applied_on]}"
                    ][1],
                    particlenetddtmap[
                        f"discddt_map_{year}_smooth_4_0p975{corrected_str[jec_applied_on]}"
                    ][2],
                ),
            )
            for jec_applied_on in ["none", "pt", "pt&mJ"]
        }

        self.trigger_scalefactors = correctionlib.CorrectionSet.from_file(
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/notebooks/data/"
            + "HLT_AK8PFJet_MC_trigger_sf_c2e731345f.json"
        )

        self.mjet_reco_correction = correctionlib.CorrectionSet.from_file(
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/"
            + jms_correction_files["notagger"]
        )

        # get some corrections and pack them into dense_lookups
        corrections_extractor = coffea.lookup_tools.extractor()
        for boson in ["W", "Z"]:
            fname = f"{kfactor_path}/{boson}JetsCorr.root"
            corrections_extractor.import_file(fname)
            corrections_extractor.add_weight_sets(
                [
                    f"{boson}_kfactor kfactor {fname}",
                    f"{boson}_ewcorr ewcorr {fname}",
                ]
            )

        corrections_extractor.finalize()

        self.corrections = corrections_extractor.make_evaluator()

        self._vjets_corrections = correctionlib.CorrectionSet.from_file(
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/"
            "ULvjets_corrections.json"
        )

        self._selections = ["vjets", "ttbar"]

        # these are selections that are linked as one specifies later in self._regions!
        self._cuts = {"vjets": ["n2ddt", "rhocut"], "ttbar": ["tau32", "tau21"]}

        # these are sample specific matching criteria
        # these are handled via 'any' of PackedSelection, so in case of multiple requirements,
        # be aware, that they are linked with OR!!
        self._matching_mappings = {
            "vjets_WJetsMatched": {"IsMergedWZ": 1},
            "vjets_WJetsUnmatched": {"IsMergedWZ": 0},
            "vjets_ZJetsMatched": {"IsMergedWZ": 1},
            "vjets_ZJetsUnmatched": {"IsMergedWZ": 0},
            "ttbar_TTToSemiLeptonic_mergedTop": {"IsMergedTop": 1},
            "ttbar_TTToSemiLeptonic_mergedW": {"IsMergedWZ": 1},
            "ttbar_TTToSemiLeptonic_mergedQB": {"IsMergedQB": 1},
            "ttbar_TTToSemiLeptonic_semiMergedTop": {"IsMergedWZ": 1, "IsMergedQB": 1},
            "ttbar_TTToSemiLeptonic_notMerged": {"IsNotMerged": 1},
        }

        # common control-plots
        hists.update(
            {
                "pt": hist.Hist(pT_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()),
                "eta": hist.Hist(eta_ax, dataset_ax, storage=hist.storage.Weight()),
                "phi": hist.Hist(phi_ax, dataset_ax, storage=hist.storage.Weight()),
                "mjet": hist.Hist(mJ_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()),
                "rho": hist.Hist(rho_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()),
                "npv": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(80, 0, 80, name="npv", label=r"$N_{PV}$"),
                    storage=hist.storage.Weight(),
                ),
                "ntrueint": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(80, 0, 80, name="ntrueint", label=r"$N_{TrueInt}$"),
                    storage=hist.storage.Weight(),
                ),
                "chf": hist.Hist(
                    dataset_ax,
                    chf_ax,
                    storage=hist.storage.Weight(),
                ),
                "nhf": hist.Hist(
                    dataset_ax,
                    nhf_ax,
                    storage=hist.storage.Weight(),
                ),
                "n2": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, -2, 2, name="n2", label="$N_{2}$"),
                    storage=hist.storage.Weight(),
                ),
                "n2ddt": hist.Hist(
                    dataset_ax,
                    jec_applied_ax,
                    hist.axis.Regular(51, -2, 2, name="n2ddt", label=r"$N_{2}^{\mathrm{DDT}}$"),
                    storage=hist.storage.Weight(),
                ),
                "tau21": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, 0, 1, name="tau21", label=r"$\frac{\tau_{2}}{\tau_{1}}$"),
                    storage=hist.storage.Weight(),
                ),
                "tau32": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, 0, 1, name="tau32", label=r"$\frac{\tau_{3}}{\tau_{2}}$"),
                    storage=hist.storage.Weight(),
                ),
                "pNet_TvsQCD": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, 0, 1, name="pNet_TvsQCD", label=r"pNet TvsQCD"),
                    storage=hist.storage.Weight(),
                ),
                "pNet_WvsQCD": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, 0, 1, name="pNet_WvsQCD", label=r"pNet WvsQCD"),
                    storage=hist.storage.Weight(),
                ),
                "pNet_MD_WvsQCD": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(51, 0, 1, name="pNet_MD_WvsQCD", label=r"pNet MD WvsQCD"),
                    storage=hist.storage.Weight(),
                ),
            }
        )

        # define regions in terms of selection-bits
        # for "old" approach using energy correlation and n-subjettiness substructure variables
        # new approach using particlenet (working points from https://indico.cern.ch/event/1152827/contributions/
        # 4840404/attachments/2428856/4162159/ParticleNet_SFs_ULNanoV9_JMAR_25April2022_PK.pdf
        tagger = {
            "vjets": {
                # "substructure": {"pass": {"n2": True}, "fail": {"n2": False}},
                "substructure": {"pass": {"n2ddt": True}, "fail": {"n2ddt": False}},
                "particlenet": {"pass": {"particlenetMDWvsQCD": True}, "fail": {"particlenetMDWvsQCD": False}},
                "particlenetDDT": {
                    "pass": {"particlenetMDWvsQCD_DDT": True},
                    "fail": {"particlenetMDWvsQCD_DDT": False}
                },
            },
            "ttbar": {
                "substructure": {
                    "pass": {"tau32": True},
                    "passW": {"tau32": False, "tau21": True},
                    "fail": {"tau32": False, "tau21": False},
                },
                "particlenet": {
                    "pass": {"particlenetTvsQCD": True},
                    "passW": {"particlenetTvsQCD": False, "particlenetWvsQCD": True},
                    "fail": {"particlenetTvsQCD": False, "particlenetWvsQCD": False},
                },
                "particlenetDDT": {
                    "pass": {"particlenetTvsQCD": True},
                    "passW": {"particlenetTvsQCD": False, "particlenetWvsQCD": True},
                    "fail": {"particlenetTvsQCD": False, "particlenetWvsQCD": False},
                },
            },
        }

        self._regions = {
            "vjets": {
                "inclusive": {"rhocut": True, "pt500cut": True, "trigger": True},
                "pass": {
                    **tagger["vjets"][self._tagger_approach]["pass"],
                    "rhocut": True, "pt500cut": True, "trigger": True
                },
                "fail": {
                    **tagger["vjets"][self._tagger_approach]["fail"],
                    "rhocut": True, "pt500cut": True, "trigger": True
                },
                # "inclusive_trigger": {
                #     "rhocut": True,
                #     "pt500cut": True,
                #     "trigger": True,
                # },
                # "pass_trigger": {
                #     **tagger["vjets"][self._tagger_approach]["pass"],
                #     "rhocut": True,
                #     "pt500cut": True,
                #     "trigger": True,
                # },
                # "fail_trigger": {
                #     **tagger["vjets"][self._tagger_approach]["fail"],
                #     "rhocut": True,
                #     "pt500cut": True,
                #     "trigger": True,
                # },
            },
            "ttbar": {
                "inclusive": {"pt200cut": True},
                "pass": {
                    **tagger["ttbar"][self._tagger_approach]["pass"],
                    "pt200cut": True
                },
                "passW": {
                    **tagger["ttbar"][self._tagger_approach]["passW"],
                    "pt200cut": True
                },
                "fail": {
                    **tagger["ttbar"][self._tagger_approach]["fail"],
                    "pt200cut": True
                },
            },
        }

        self._variations = [
            "0_0_all",
            "0_0_chargedH",
            "0_0_gamma",
            "0_0_neutralH",
            "0_0_other",
        ]
        for selection in self._selections:
            for region in self._regions[selection].keys():
                hists.update(
                    {
                        f"{selection}_mjet_{region}": hist.Hist(
                            mJ_fit_ax,
                            self._pT_fit_ax[selection],
                            dataset_ax,
                            jec_applied_ax,
                            eta_regions_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_mPnet_{region}": hist.Hist(
                            mPnet_fit_ax,
                            self._pT_fit_ax[selection],
                            dataset_ax,
                            jec_applied_ax,
                            eta_regions_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_pt_{region}": hist.Hist(
                            pT_ax,
                            dataset_ax,
                            jec_applied_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_eta_{region}": hist.Hist(
                            eta_ax, dataset_ax, storage=hist.storage.Weight()
                        ),
                        f"{selection}_chf_{region}": hist.Hist(
                            dataset_ax,
                            chf_ax,
                            self._pT_fit_ax[selection],
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_nhf_{region}": hist.Hist(
                            dataset_ax,
                            nhf_ax,
                            self._pT_fit_ax[selection],
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_rho_{region}": hist.Hist(
                            rho_ax,
                            dataset_ax,
                            jec_applied_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_mjet_v_jecfactor_{region}": hist.Hist(
                            mJ_ax,
                            eta_regions_ax,
                            hist.axis.Regular(
                                300, 0, 3, name="jecfactor", label="jecfactor",
                            ),
                            dataset_ax,
                            storage=hist.storage.Weight(),
                        ),
                    }
                )

                hists.update(
                    {
                        f"{selection}_mjet_unfolding_{region}": hist.Hist(
                            self._unfolding_ax[selection]["mJgen"],
                            self._unfolding_ax[selection]["ptgen"],
                            self._unfolding_ax[selection]["mJreco"],
                            self._unfolding_ax[selection]["ptreco"],
                            dataset_ax,
                            jec_applied_ax,
                            fakes_ax,
                            storage=hist.storage.Weight(),
                        ),
                    }
                )

                hists.update(
                    {
                        f"{selection}_mjetgen_unfolding_{region}": hist.Hist(
                            self._unfolding_ax[selection]["mJgen"],
                            self._unfolding_ax[selection]["ptgen"],
                            dataset_ax,
                            storage=hist.storage.Weight(),
                        ),
                    }
                )

                for variation in self._variations:
                    hists.update(
                        {
                            f"{selection}_mjet_{variation}_variation_{region}__up": hist.Hist(
                                mJ_fit_ax,
                                self._pT_fit_ax[selection],
                                eta_regions_ax,
                                dataset_ax,
                                jec_applied_ax,
                                storage=hist.storage.Weight(),
                            ),
                        }
                    )

                    hists.update(
                        {
                            f"{selection}_mjet_{variation}_variation_{region}__down": hist.Hist(
                                mJ_fit_ax,
                                self._pT_fit_ax[selection],
                                eta_regions_ax,
                                dataset_ax,
                                jec_applied_ax,
                                storage=hist.storage.Weight(),
                            ),
                        }
                    )
                    hists.update(
                        {
                            f"{selection}_mjet_unfolding_{variation}_variation_{region}__up": hist.Hist(
                                self._unfolding_ax[selection]["mJgen"],
                                self._unfolding_ax[selection]["ptgen"],
                                self._unfolding_ax[selection]["mJreco"],
                                self._unfolding_ax[selection]["ptreco"],
                                dataset_ax,
                                jec_applied_ax,
                                storage=hist.storage.Weight(),
                            ),
                        }
                    )
                    hists.update(
                        {
                            f"{selection}_mjet_unfolding_{variation}_variation_{region}__down": hist.Hist(
                                self._unfolding_ax[selection]["mJgen"],
                                self._unfolding_ax[selection]["ptgen"],
                                self._unfolding_ax[selection]["mJreco"],
                                self._unfolding_ax[selection]["ptreco"],
                                dataset_ax,
                                jec_applied_ax,
                                storage=hist.storage.Weight(),
                            ),
                        }
                    )

                hists.update(
                    {
                        f"{selection}_mPnet_0_0_all_variation_{region}__up": hist.Hist(
                            mPnet_fit_ax,
                            self._pT_fit_ax[selection],
                            eta_regions_ax,
                            dataset_ax,
                            jec_applied_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_mPnet_0_0_all_variation_{region}__down": hist.Hist(
                            mPnet_fit_ax,
                            self._pT_fit_ax[selection],
                            eta_regions_ax,
                            dataset_ax,
                            jec_applied_ax,
                            storage=hist.storage.Weight(),
                        )
                    }
                )

        self._hists = lambda: {
            **hists,
            "nevents": processor.defaultdict_accumulator(float),
            "sumw": processor.defaultdict_accumulator(float),
            "sumw2": processor.defaultdict_accumulator(float),
        }

        self._triggerbits = [
            "HLT_PFJet320_v*",
            "HLT_PFJet400_v*",
            "HLT_PFJet450_v*",
            "HLT_PFJet500_v*",
            "HLT_PFJet550_v*",
            "HLT_AK8PFJet320_v*",
            "HLT_AK8PFJet400_v*",
            "HLT_AK8PFJet450_v*",
            "HLT_AK8PFJet500_v*",
            "HLT_AK8PFJet550_v*",
        ]

    def passes_trigger(self, events, trigger_name):
        return events["trigger_bits"][:, self._triggerbits.index(trigger_name)]

    @property
    def accumulator(self):
        return self._hists

    def n2ddt(self, pt, rho, n2, corrected="none"):
        quantile = self._n2ddtmaps[corrected](rho, pt)
        return n2 - quantile

    def pNetMDWvsQCDddt(self, pt, rho, val, corrected="none"):
        quantile = self._pNetMDWvsQCDddtmaps[corrected](rho, pt)
        return quantile - val

    def vjets_syst_weights(self, syst, vpt, boson="W"):
        evaluator = self._vjets_corrections[f"{boson}_FixedOrderComponent"]
        nominal = evaluator.evaluate("nominal", vpt)
        return evaluator.evaluate(syst, vpt)/nominal

    def vjets_syst_weights_envelope(self, systs, vpt, boson="W", edge="upper"):
        vals = np.array([self.vjets_syst_weights(f"{syst}", vpt, boson) for syst in systs])
        envelope_edge = np.max(vals, axis=0) if edge == "upper" else np.min(vals, axis=0)
        return envelope_edge

    def postprocess(self, accumulator):
        return accumulator

    def process(self, events):
        out = self.accumulator()

        dataset = events.metadata["dataset"]
        selection = events.metadata["selection"]

        isMC = "data" not in dataset.lower()

        # HEM15/16 Treatment
        # https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
        # courtesy of C.Matthies
        # https://github.com/MatthiesC/LegacyTopTagging/blob/master/include/Utils.h#L325-L328
        HEM_affected_lumi_fraction = 0.64844705699  # (Run 319077 (17.370008/pb) + Run C + Run D) / all 2018
        HEM_eta_min, HEM_eta_max = -3.2, -1.3
        HEM_phi_min, HEM_phi_max = -1.57, -0.87
        HEM_affected_event_jets = (
            (events.eta < HEM_eta_max)
            & (events.eta > HEM_eta_min)
            & (events.phi < HEM_phi_max)
            & (events.phi > HEM_phi_min)
        )
        treat_HEM = False
        event_weight_for_gen = deepcopy(events.weight)
        if isMC:
            if self._year == "UL18":
                events["weight"] = ak.where(
                    HEM_affected_event_jets, events.weight * (1 - HEM_affected_lumi_fraction), events.weight
                )
        else:
            if self._year == "UL18" and any(run in events.metadata["filename"] for run in ["RunC", "RunD"]):
                treat_HEM = True

        # evaluate matching criteria and created new masked events dataframe
        matching_mask = np.ones(len(events), dtype="bool")

        if dataset in self._matching_mappings.keys():

            matching_selection = PackedSelection()
            for branch_name, branch_value in self._matching_mappings[dataset].items():

                matching_selection.add(branch_name, events[branch_name] == branch_value)
            # BE AWARE OF THE FOLLOWING OR!!!!
            matching_mask = matching_selection.any(
                *self._matching_mappings[dataset].keys()
            )

        events = events[matching_mask]
        event_weight_for_gen = event_weight_for_gen[matching_mask]

        out["nevents"][dataset] = len(events)
        out["sumw"][dataset] = ak.sum(events.weight)
        out["sumw2"][dataset] = ak.sum(events.weight * events.weight)

        # right now its commented since i apply the kfactor when making the flat tree in UHH2
        # genpt = events.genjetpt
        # vpt = events.V_pt
        # if('WJets' in dataset):
        #     events.weight = events.weight/self.corrections['W_kfactor'](genpt)*self.corrections[f'W_ewcorr'](vpt)
        #     events.weight = events.weight*self.corrections['W_kfactor'](genpt)*self.corrections[f'W_ewcorr'](genpt)
        # if('ZJets' in dataset):
        #     events.weight = events.weight/self.corrections['Z_kfactor'](genpt)*self.corrections[f'Z_ewcorr'](vpt)
        #     events.weight = events.weight*self.corrections['Z_kfactor'](genpt)*self.corrections[f'Z_ewcorr'](genpt)

        jecfactors = {
            "nominal": events.jecfactor,
            "up": events.jecfactor_up,
            "down": events.jecfactor_down,
        }
        if "data" in dataset.lower():
            jecfactor = jecfactors["nominal"]
        else:
            jecfactor = jecfactors[self._jec]

        # apply top-pt reweighting weight
        events["weight"] = events.weight * events["toppt_weight"]
        if self._variation_weight != "nominal" and "data" not in dataset.lower() and len(events) > 0:
            variation_weights = {
                "toppt_off": 1. / events["toppt_weight"],
                "pu_down": events["weight_pu_down"]/events["weight_pu"],
                "pu_up": events["weight_pu_up"]/events["weight_pu"],
            }
            if len(events["ps_weights"][0]) == 46:
                variation_weights["fsr_down"] = events["ps_weights"][:, 4] / events["ps_weights"][:, 0]
                variation_weights["fsr_up"] = events["ps_weights"][:, 5] / events["ps_weights"][:, 0]
                variation_weights["isr_down"] = events["ps_weights"][:, 26] / events["ps_weights"][:, 0]
                variation_weights["isr_up"] = events["ps_weights"][:, 27] / events["ps_weights"][:, 0]

            else:
                variation_weights["fsr_down"] = 1.0
                variation_weights["fsr_up"] = 1.0
                variation_weights["isr_down"] = 1.0
                variation_weights["isr_up"] = 1.0

            if "WJets" in dataset or "ZJets" in dataset:
                boson = "W" if "W" in dataset else "Z"
                v_qcd_systs = [
                    f"{syst}_{direction}" for direction in ["up", "down"] for syst in ["d1K_NLO", "d2K_NLO", "d3K_NLO"]
                ]
                variation_weights["v_qcd_down"] = self.vjets_syst_weights_envelope(
                    v_qcd_systs,
                    events["V_pt"],
                    edge="lower",
                    boson=boson,
                )
                variation_weights["v_qcd_up"] = self.vjets_syst_weights_envelope(
                    v_qcd_systs,
                    events["V_pt"],
                    edge="upper",
                    boson=boson,
                )

                w_ewk_systs = [
                    f"{syst}_{direction}"
                    for direction in ["up", "down"]
                    for syst in ["d1kappa_EW", "W_d2kappa_EW", "W_d3kappa_EW"]
                ]
                z_ewk_systs = [
                    f"{syst}_{direction}"
                    for direction in ["up", "down"]
                    for syst in ["d1kappa_EW", "Z_d2kappa_EW", "Z_d3kappa_EW"]
                ]
                if boson == "W":
                    variation_weights["w_ewk_down"] = self.vjets_syst_weights_envelope(
                        w_ewk_systs,
                        events["V_pt"],
                        edge="lower",
                        boson=boson,
                    )
                    variation_weights["w_ewk_up"] = self.vjets_syst_weights_envelope(
                        w_ewk_systs,
                        events["V_pt"],
                        edge="upper",
                        boson=boson,
                    )
                    variation_weights["z_ewk_down"] = 1.0
                    variation_weights["z_ewk_up"] = 1.0
                else:
                    variation_weights["w_ewk_down"] = 1.0
                    variation_weights["w_ewk_up"] = 1.0
                    variation_weights["z_ewk_down"] = self.vjets_syst_weights_envelope(
                        z_ewk_systs,
                        events["V_pt"],
                        edge="lower",
                        boson=boson,
                    )
                    variation_weights["z_ewk_up"] = self.vjets_syst_weights_envelope(
                        z_ewk_systs,
                        events["V_pt"],
                        edge="upper",
                        boson=boson,
                    )
            else:
                variation_weights["v_qcd_down"] = 1.0
                variation_weights["v_qcd_up"] = 1.0
                variation_weights["w_ewk_down"] = 1.0
                variation_weights["w_ewk_up"] = 1.0
                variation_weights["z_ewk_down"] = 1.0
                variation_weights["z_ewk_up"] = 1.0

            events["weight"] = events.weight * variation_weights[self._variation_weight]
            event_weight_for_gen = event_weight_for_gen * variation_weights[self._variation_weight]
        pt_raw = events.pt
        pt = pt_raw * jecfactor

        ptgen_ = events.pt_gen_ak8

        mjet_raw = events.mjet
        mjet = mjet_raw * jecfactor

        mPnet_raw = events["ParticleNetMassRegression_mass"]
        mPnet = mPnet_raw * jecfactor

        mJgen_ = events.msd_gen_ak8

        # # not taking leading particle level jet but rather order by either dR(reco, gen) or N2(gen)
        # if "WJetsMatched" in dataset:
        #     if self._gen_sort == "n2":
        #         ptgen_ = ak.where(
        #             events.gentopjet_n2_0 < events.gentopjet_n2_1, events.gentopjet_pt_0, events.gentopjet_pt_1
        #         )
        #         mJgen_ = ak.where(
        #             events.gentopjet_n2_0 < events.gentopjet_n2_1, events.gentopjet_msd_0, events.gentopjet_msd_1
        #         )
        #     elif self._gen_sort == "dR":
        #         ptgen_ = ak.where(
        #             events.gentopjet_dR_reco_0 < events.gentopjet_dR_reco_1,
        #             events.gentopjet_pt_0,
        #             events.gentopjet_pt_1,
        #         )
        #         mJgen_ = ak.where(
        #             events.gentopjet_dR_reco_0 < events.gentopjet_dR_reco_1,
        #             events.gentopjet_msd_0,
        #             events.gentopjet_msd_1,
        #         )

        rho = 2 * np.log(mjet / pt)
        rho_raw = 2 * np.log(mjet_raw / pt_raw)
        rho_corrected_pt = 2 * np.log(mjet_raw / pt)
        rho_corrected_mJ = 2 * np.log(mjet / pt_raw)

        eta_ = events.eta
        phi_ = events.phi

        # apply trigger sf
        if isMC and selection == "vjets":
            trigger_sf_evaluator_450 = self.trigger_scalefactors[
                f"HLT_AK8PFJet450_triggersf_{self._year}"
            ]
            trigger_sf_evaluator_500 = self.trigger_scalefactors[
                f"HLT_AK8PFJet500_triggersf_{self._year}"
            ]
            first_ptbin = (pt < 650.0)
            events["weight"] = events.weight * ak.where(
                first_ptbin,
                trigger_sf_evaluator_450.evaluate(
                    pt, self._trigger_sf_variation
                ),
                trigger_sf_evaluator_500.evaluate(
                    pt, self._trigger_sf_variation
                )
            )

        events["ParticleNetMDDiscriminators_WvsQCD"] = (
            events["ParticleNetMD_probXqq"] + events["ParticleNetMD_probXcc"]
        ) / (events["ParticleNetMD_probXqq"] + events["ParticleNetMD_probXcc"] + events["ParticleNetMD_probQCD"])

        # apply trigger selection of vjets also to control plots
        m_cps = np.ones_like(events.pt, dtype=bool)
        vjets_trigger_mask = ak.where(
            pt < 650.0,
            events["trigger_bits"][:, self._triggerbits.index("HLT_AK8PFJet450_v*")] == 1,
            events["trigger_bits"][:, self._triggerbits.index("HLT_AK8PFJet500_v*")] == 1,
        )

        if selection == "vjets":
            m_cps = m_cps & vjets_trigger_mask

        out["pt"].fill(dataset=dataset, jecAppliedOn="pt", pt=pt[m_cps], weight=events.weight[m_cps])
        out["pt"].fill(dataset=dataset, jecAppliedOn="none", pt=pt_raw[m_cps], weight=events.weight[m_cps])
        out["eta"].fill(dataset=dataset, eta=eta_[m_cps], weight=events.weight[m_cps])
        out["phi"].fill(dataset=dataset, phi=phi_[m_cps], weight=events.weight[m_cps])
        out["mjet"].fill(dataset=dataset, jecAppliedOn="mJ", mJ=mjet[m_cps], weight=events.weight[m_cps])
        out["mjet"].fill(dataset=dataset, jecAppliedOn="none", mJ=mjet_raw[m_cps], weight=events.weight[m_cps])
        out["rho"].fill(dataset=dataset, jecAppliedOn="pt&mJ", rho=rho[m_cps], weight=events.weight[m_cps])
        out["rho"].fill(
            dataset=dataset,
            jecAppliedOn="pt",
            rho=rho_corrected_pt[m_cps],
            weight=events.weight[m_cps],
        )
        out["rho"].fill(
            dataset=dataset,
            jecAppliedOn="mJ",
            rho=rho_corrected_mJ[m_cps],
            weight=events.weight[m_cps],
        )
        out["rho"].fill(dataset=dataset, jecAppliedOn="none", rho=rho_raw[m_cps], weight=events.weight[m_cps])

        out["npv"].fill(dataset=dataset, npv=events.n_pv[m_cps], weight=events.weight[m_cps])
        if isMC:
            out["ntrueint"].fill(dataset=dataset, ntrueint=events.n_trueint[m_cps], weight=events.weight[m_cps])
        out["chf"].fill(dataset=dataset, chf=events.CHF[m_cps], weight=events.weight[m_cps])
        out["nhf"].fill(dataset=dataset, nhf=events.NHF[m_cps], weight=events.weight[m_cps])

        out["n2"].fill(dataset=dataset, n2=events.N2[m_cps], weight=events.weight[m_cps])
        out["tau21"].fill(dataset=dataset, tau21=events.tau21[m_cps], weight=events.weight[m_cps])
        out["tau32"].fill(dataset=dataset, tau32=events.tau32[m_cps], weight=events.weight[m_cps])
        out["pNet_WvsQCD"].fill(
            dataset=dataset, pNet_WvsQCD=events["ParticleNetDiscriminators_WvsQCD"][m_cps], weight=events.weight[m_cps]
        )
        out["pNet_MD_WvsQCD"].fill(
            dataset=dataset,
            pNet_MD_WvsQCD=events["ParticleNetMDDiscriminators_WvsQCD"][m_cps],
            weight=events.weight[m_cps],
        )
        out["pNet_TvsQCD"].fill(
            dataset=dataset, pNet_TvsQCD=events["ParticleNetDiscriminators_TvsQCD"][m_cps], weight=events.weight[m_cps]
        )
        # for jec_applied_on in ['none','pt','pt&mJ']:
        for jec_applied_on in ["pt", "pt&mJ"]:
            selections = PackedSelection()
            pt_ = pt_raw
            if "pt" in jec_applied_on:
                pt_ = pt
            mJ_ = mjet_raw
            mPnet_ = mPnet_raw

            if "mJ" in jec_applied_on:
                mJ_ = mjet
                mPnet_ = mPnet

            rho_ = 2 * np.log(mJ_ / pt_)

            selections.add("pt500cut", (pt_ > 500))
            selections.add("pt200cut", (pt_ > 200))

            selections.add(
                "n2ddt",
                (
                    events.N2 > 0
                )
                # make sure to fail on events with default values for topjets
                # (do we need to remove those from fail region as well?)
                & (
                    self.n2ddt(pt_, rho_, events.N2, corrected=jec_applied_on) < 0
                ),  # actual N2-DDT tagger
            )
            out["n2ddt"].fill(
                dataset=dataset,
                jecAppliedOn=jec_applied_on,
                n2ddt=self.n2ddt(pt_, rho_, events.N2, corrected=jec_applied_on),
                weight=events.weight,
            )

            selections.add(
                "n2",
                (events.N2 > 0) & (events.N2 < 0.2)
            )

            # selections.add("rhocut",
            #                (rho_<-2.1)
            #                &(rho_>-6.0))
            selections.add("rhocut", rho_ < -2.1)
            # selections.add("rhocut", np.ones_like(events.pt, dtype=bool))

            selections.add("tau21", events.tau21 < 0.45)
            selections.add("tau32", events.tau32 < 0.5)

            selections.add("particlenetWvsQCD", events["ParticleNetMDDiscriminators_WvsQCD"] > 0.91)
            # selections.add("particlenetWvsQCD", events["ParticleNetDiscriminators_WvsQCD"] > 0.97)
            selections.add("particlenetTvsQCD", events["ParticleNetDiscriminators_TvsQCD"] > 0.96)
            selections.add("particlenetMDWvsQCD", events["ParticleNetMDDiscriminators_WvsQCD"] > 0.91)
            selections.add(
                "particlenetMDWvsQCD_DDT",
                (events["ParticleNetMDDiscriminators_WvsQCD"] > 0)
                & (
                    self.pNetMDWvsQCDddt(
                        pt_, rho_, events["ParticleNetMDDiscriminators_WvsQCD"], corrected=jec_applied_on
                    )
                    < 0
                ),
            )
            # selections.add(
            #     "gensel_drmatch",
            #     (events.pass_gen_selection == 1) & (events.dR_reco_gen < 0.4) & (mJgen_ > 30.),
            # )

            if "n2_beta1_gen" not in events.fields:
                events.n2_beta1_gen = ak.ones_like(events.pt)

            selections.add(
                "gensel_drmatch",
                (events.pass_gen_selection == 1)
                & (events.dR_reco_gen < 0.4)
                & (mJgen_ > 30.)
                & (events.n2_beta1_gen < 0.2),
            )
            selections.add(
                "recosel",
                (events.pass_reco_selection == 1)
            )

            passing_hem_treatment = (
                ~HEM_affected_event_jets if treat_HEM else
                ak.ones_like(events.pt, dtype=bool)
            )

            selections.add("jetpfid",
                           (
                               (events.jetpfid == 1)
                               & (passing_hem_treatment == 1)
                           )
                           )
            selection = events.metadata["selection"]

            selections.add(
                "trigger",
                vjets_trigger_mask
            )

            for region in self._regions[selection].keys():
                smask_unfolding = selections.require(**self._regions[selection][region], recosel=True, jetpfid=True)
                selection_no_trigger = deepcopy(self._regions[selection][region])
                if "trigger" in selection_no_trigger:
                    selection_no_trigger.pop("trigger")
                smask_unfolding_gen = selections.require(**selection_no_trigger, recosel=True, jetpfid=True)
                corrector = self.mjet_reco_correction[
                    "response_g_" + ("jec" if "mJ" in jec_applied_on else "nojec") + f"_{self._year}"
                ]
                msd_correction = (
                    1.0/corrector.evaluate(pt_[smask_unfolding])  # divide since this is the response msdrec/mgen
                    if selection == "vjets"
                    else 1.0  # ak.ones_like(pt_[smask_unfolding])
                )

                # smask_unfolding_phasespace = selections.require(unfolding=True)
                out[f"{selection}_mjet_unfolding_{region}"].fill(
                    ptreco=pt_[smask_unfolding],
                    mJreco=mJ_[smask_unfolding] * msd_correction,
                    ptgen=ptgen_[smask_unfolding],
                    mJgen=mJgen_[smask_unfolding],
                    dataset=dataset,
                    fakes=selections.require(gensel_drmatch=False)[smask_unfolding]
                    if "WJetsMatched" in dataset
                    else np.zeros_like(smask_unfolding)[smask_unfolding],
                    jecAppliedOn=jec_applied_on,
                    weight=events.weight[smask_unfolding],
                )
                for split_size in [90.0, 80.0, 60.0]:
                    if "WJetsMatched" in dataset:
                        split_dataset = np.where(
                            np.random.rand(len(smask_unfolding[smask_unfolding])) < (split_size / 100.),
                            "vjets_WJetsMatched0p%i" % (split_size/10.),
                            "vjets_WJetsMatched0p%i" % (10 - (split_size/10.)),
                        )
                        out[f"{selection}_mjet_unfolding_{region}"].fill(
                            ptreco=pt_[smask_unfolding],
                            mJreco=mJ_[smask_unfolding] * msd_correction,
                            ptgen=ptgen_[smask_unfolding],
                            mJgen=mJgen_[smask_unfolding],
                            dataset=split_dataset,
                            fakes=selections.require(gensel_drmatch=False)[smask_unfolding],
                            jecAppliedOn=jec_applied_on,
                            weight=events.weight[smask_unfolding],
                        )

                # out[f"{selection}_mjet_unfolding_{region}"].fill(
                #     ptreco=pt_[smask_unfolding & smask_unfolding_phasespace],
                #     mJreco=mJ_[smask_unfolding & smask_unfolding_phasespace] * msd_correction,
                #     ptgen=ptgen_[smask_unfolding & smask_unfolding_phasespace],
                #     mJgen=mJgen_[smask_unfolding & smask_unfolding_phasespace],
                #     dataset=dataset,
                #     fakes="fakes",
                #     jecAppliedOn=jec_applied_on,
                #     weight=events.weight[smask_unfolding],
                # )

                out[f"{selection}_mjetgen_unfolding_{region}"].fill(
                    ptgen=ptgen_[smask_unfolding_gen],
                    mJgen=mJgen_[smask_unfolding_gen],
                    dataset=dataset,
                    weight=event_weight_for_gen[smask_unfolding_gen],
                )

                smask = selections.require(**self._regions[selection][region], jetpfid=True)

                out[f"{selection}_mjet_{region}"].fill(
                    dataset=dataset,
                    jecAppliedOn=jec_applied_on,
                    pt=pt_[smask],
                    abs_eta_regions=np.abs(eta_[smask]),
                    mJ=mJ_[smask],
                    weight=events.weight[smask],
                )

                out[f"{selection}_mPnet_{region}"].fill(
                    dataset=dataset,
                    jecAppliedOn=jec_applied_on,
                    pt=pt_[smask],
                    abs_eta_regions=np.abs(eta_[smask]),
                    mPnet=mPnet_[smask],
                    weight=events.weight[smask],
                )

                out[f"{selection}_pt_{region}"].fill(
                    dataset=dataset,
                    jecAppliedOn=jec_applied_on,
                    pt=pt_[smask],
                    weight=events.weight[smask],
                )

                out[f"{selection}_rho_{region}"].fill(
                    dataset=dataset,
                    jecAppliedOn=jec_applied_on,
                    rho=rho_[smask],
                    weight=events.weight[smask],
                )

                if jec_applied_on == "pt":
                    out[f"{selection}_mjet_v_jecfactor_{region}"].fill(
                        dataset=dataset,
                        mJ=mJ_[smask],
                        abs_eta_regions=np.abs(eta_[smask]),
                        jecfactor=jecfactor[smask],
                        weight=events.weight[smask]
                    )

                if jec_applied_on == "pt":
                    out[f"{selection}_chf_{region}"].fill(
                        dataset=dataset, chf=events.CHF[smask], pt=pt_[smask], weight=events.weight[smask]
                    )
                    out[f"{selection}_nhf_{region}"].fill(
                        dataset=dataset, nhf=events.NHF[smask], pt=pt_[smask], weight=events.weight[smask]
                    )
                    out[f"{selection}_eta_{region}"].fill(dataset=dataset, eta=eta_[smask], weight=events.weight[smask])

                if isMC:
                    for variation in self._variations:
                        mJVar_ = events[f"mjet_{variation}"]

                        if "mJ" in jec_applied_on:
                            mJVar_ = mJVar_ * jecfactor

                        out[f"{selection}_mjet_{variation}_variation_{region}__up"].fill(
                            dataset=dataset,
                            jecAppliedOn=jec_applied_on,
                            pt=pt_[smask],
                            mJ=mJVar_[:, 0][smask],
                            abs_eta_regions=np.abs(eta_[smask]),
                            weight=events.weight[smask],
                        )
                        out[f"{selection}_mjet_{variation}_variation_{region}__down"].fill(
                            dataset=dataset,
                            jecAppliedOn=jec_applied_on,
                            pt=pt_[smask],
                            mJ=mJVar_[:, 1][smask],
                            abs_eta_regions=np.abs(eta_[smask]),
                            weight=events.weight[smask],
                        )
                        # out[f"{selection}_mjet_unfolding_{variation}_variation_{region}__up"].fill(
                        #     ptreco=pt_[smask_unfolding],
                        #     mJreco=mJVar_[:, 0][smask_unfolding] * msd_correction,
                        #     ptgen=ptgen_[smask_unfolding],
                        #     mJgen=mJgen_[smask_unfolding],
                        #     dataset=dataset,
                        #     jecAppliedOn=jec_applied_on,
                        #     weight=events.weight[smask_unfolding],
                        # )
                        # out[f"{selection}_mjet_unfolding_{variation}_variation_{region}__down"].fill(
                        #     ptreco=pt_[smask_unfolding],
                        #     mJreco=mJVar_[:, 1][smask_unfolding] * msd_correction,
                        #     ptgen=ptgen_[smask_unfolding],
                        #     mJgen=mJgen_[smask_unfolding],
                        #     dataset=dataset,
                        #     jecAppliedOn=jec_applied_on,
                        #     weight=events.weight[smask_unfolding],
                        # )

                    out[f"{selection}_mPnet_0_0_all_variation_{region}__up"].fill(
                        dataset=dataset,
                        jecAppliedOn=jec_applied_on,
                        pt=pt_[smask],
                        mPnet=mPnet_[smask]*1.005,
                        abs_eta_regions=np.abs(eta_[smask]),
                        weight=events.weight[smask],
                    )
                    out[f"{selection}_mPnet_0_0_all_variation_{region}__down"].fill(
                        dataset=dataset,
                        jecAppliedOn=jec_applied_on,
                        pt=pt_[smask],
                        mPnet=mPnet_[smask]*0.995,
                        abs_eta_regions=np.abs(eta_[smask]),
                        weight=events.weight[smask],
                    )

        return out


if __name__ == "__main__":
    workflow = CoffeaWorkflow("JMSTemplates")

    workflow.parser.add_argument("--output", "-o", type=str, default="jms_templates.coffea")
    workflow.parser.add_argument("--year", default="UL17")
    workflow.parser.add_argument("--JEC", default="nominal")
    workflow.parser.add_argument("--variation", default="nominal", choices=[
        "isr_up", "isr_down",
        "fsr_up", "fsr_down",
        "pu_up", "pu_down",
        "toppt_off",
        "v_qcd_up", "v_qcd_down",
        "w_ewk_up", "w_ewk_down",
        "z_ewk_up", "z_ewk_down",
    ])
    workflow.parser.add_argument("--triggersf", default="nominal", choices=["nominal", "up", "down"])
    workflow.parser.add_argument("--maxfiles", type=int, default=-1)
    workflow.parser.add_argument(
        "--tagger", default="substructure", choices=["substructure", "particlenet", "particlenetDDT"]
    )

    workflow.parser.add_argument("--VJetsOnly", action="store_true")

    args = workflow.parse_args()

    workflow.processor_instance = JMSTemplates(
        year=args.year,
        jec=args.JEC,
        variation_weight=args.variation,
        trigger_sf_var=args.triggersf,
        tagger=args.tagger
    )
    workflow.processor_schema = BaseSchema

    sample_pattern = (
        "/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/{SELECTION}Trees/workdir_{SELECTION}_{YEAR}/*{SAMPLE}*.root"
    )
    sample_names = {
        "vjets": [
            "Data",
            "WJets",
            "WJetsMatched",
            "WJetsUnmatched",
            "ZJets",
            "ZJetsMatched",
            "ZJetsUnmatched",
            "TTToHadronic",
            "TTToSemiLeptonic",
            "ST_tW_top",
            "ST_tW_antitop",
            "QCD",
        ],
        # "vjets":["Data.JetHT_RunB","Data.JetHT_RunC","Data.JetHT_RunD","Data.JetHT_RunE","Data.JetHT_RunF",],
        "ttbar": [
            "Data",
            "WJets",
            "DYJets",
            "TTToHadronic",
            "TTToSemiLeptonic",
            "TTToSemiLeptonic_mergedTop",
            "TTToSemiLeptonic_mergedW",
            "TTToSemiLeptonic_mergedQB",
            "TTToSemiLeptonic_semiMergedTop",
            "TTToSemiLeptonic_notMerged",
            "TTTo2L2Nu",
            "ST_t",
            "ST_tW",
            "ST_s",
            "QCD",
        ],
    }
    if args.VJetsOnly:
        sample_names = {
            "vjets": [
                "WJets",
                "WJetsMatched",
                "WJetsUnmatched",
                "ZJets",
                "ZJetsMatched",
                "ZJetsUnmatched",
            ]
        }

    files = {}
    # for selection in ["vjets", "ttbar"]:
    for selection in sample_names.keys():
        def parent_samplename(s):
            return min(
                [
                    s.replace(child_samplename_suffix, "")
                    for child_samplename_suffix in [
                        "Matched",
                        "Unmatched",
                        "_mergedTop",
                        "_mergedW",
                        "_mergedQB",
                        "_semiMergedTop",
                        "notMerged",
                    ]
                ],
                key=len,
            )

        samples = {
            sample: glob.glob(
                sample_pattern.format(
                    SELECTION=selection,
                    YEAR=args.year,
                    SAMPLE=parent_samplename(sample),
                )
            )
            for sample in sample_names[selection]
        }

        for k, v in samples.items():
            if k not in files:
                files[f"{selection}_{k}"] = {
                    "files": [],
                    "treename": "AnalysisTree",
                    "metadata": {"selection": selection},
                }
            if args.maxfiles > 0:
                files[f"{selection}_{k}"]["files"] += v[: args.maxfiles]
            else:
                files[f"{selection}_{k}"]["files"] += v
            print(f"{selection}_{k}", len(files[f"{selection}_{k}"]["files"]))

    output_file_path = os.path.join(os.getcwd(), args.output)
    print("changing into /tmp dir")
    os.chdir(os.environ["TMPDIR"])
    if args.scaleout > 0:
        print("init dask client")
        # if(workflow.args.debug):
        #     workflow.init_dask_local_client()
        # else:
        # q = __import__("functools").partial(__import__("os")._exit, 0)  # FIXME
        # __import__("IPython").embed()  # FIXME
        workflow.processor_args["chunksize"] = 500000
        workflow.init_dask_htcondor_client(1, 2, 5)

    print("starting coffea runner")

    output = workflow.run(files)

    if not output_file_path.endswith(".coffea"):
        output_file_path += ".coffea"

    save(output, output_file_path)
