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

jetmass_path = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass"
ddtmaps_path = f"{jetmass_path}/Histograms/ddtmaps.npy"
kfactor_path = f"{jetmass_path}/NLOweights/"


class JMSTemplates(processor.ProcessorABC):
    def __init__(self, year="2017", jec="nominal"):
        self._year = year
        self._jec = jec

        dataset_ax = hist.axis.StrCategory([], name="dataset", growth=True)
        # shift_ax = hist.axis.StrCategory([], name="shift", growth=True)
        jec_applied_ax = hist.axis.StrCategory(
            [], name="jecAppliedOn", label="JEC applied on", growth=True
        )
        # jec_ax = hist.axis.StrCategory(
        #     ["raw", "pt", "pt_up", "pt_down","pt_","pt_mJ_up", "down"], name="JEC", label="JEC"
        # )
        mJ_ax = hist.axis.Regular(50, 0.0, 500.0, name="mJ", label=r"$m_{SD}$ [GeV]")
        pT_ax = hist.axis.Regular(100, 0.0, 3000.0, name="pt", label=r"$p_{T}$ [GeV]")
        eta_ax = hist.axis.Regular(100, -6.5, 6.5, name="eta", label=r"$\eta$")
        eta_regions_ax = hist.axis.Variable([0, 1.3, 2.5], name="abs_eta_regions", label=r"$|\eta|$")
        phi_ax = hist.axis.Regular(100, -4, 4, name="phi", label=r"$\Phi$")

        mJ_fit_ax = hist.axis.Regular(
            500, 0.0, 500.0, name="mJ", label=r"$m_{SD}$ [GeV]"
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
                    np.array([0.0, 70., 80.5, 89.5, np.inf]),
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
        n2ddtmap = np.load(ddtmaps_path, allow_pickle=True).item()
        corrected_str = {
            "none": "",
            "pt": "_corrected_pt",
            "pt&mJ": "_corrected_pt_mass",
        }

        self._n2ddtmaps = {
            jec_applied_on: coffea.lookup_tools.dense_lookup.dense_lookup(
                n2ddtmap[
                    f"n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                ][0],
                dims=(
                    n2ddtmap[
                        f"n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                    ][1],
                    n2ddtmap[
                        f"n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}"
                    ][2],
                ),
            )
            for jec_applied_on in ["none", "pt", "pt&mJ"]
        }

        self.mjet_reco_correction = correctionlib.CorrectionSet.from_file(
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/"
            + "jms_corrections_quadratic_40c365c4ab.json"
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
                "pt": hist.Hist(
                    pT_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()
                ),
                "eta": hist.Hist(eta_ax, dataset_ax, storage=hist.storage.Weight()),
                "phi": hist.Hist(phi_ax, dataset_ax, storage=hist.storage.Weight()),
                "mjet": hist.Hist(
                    mJ_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()
                ),
                "rho": hist.Hist(
                    rho_ax, dataset_ax, jec_applied_ax, storage=hist.storage.Weight()
                ),
                "npv": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(80, 0, 80, name="npv", label=r"$N_{PV}$"),
                    storage=hist.storage.Weight(),
                ),
                "ntrueint": hist.Hist(
                    dataset_ax,
                    hist.axis.Regular(
                        80, 0, 80, name="ntrueint", label=r"$N_{TrueInt}$"
                    ),
                    storage=hist.storage.Weight(),
                ),
            }
        )

        self._regions = {
            "vjets": {
                "inclusive": {"rhocut": True, "pt500cut": True},
                "pass": {"n2ddt": True, "rhocut": True, "pt500cut": True},
                "fail": {"n2ddt": False, "rhocut": True, "pt500cut": True},
                "inclusive_trigger": {
                    "rhocut": True,
                    "pt500cut": True,
                    "HLT_AK8PFJet450": True,
                },
                "pass_trigger": {
                    "n2ddt": True,
                    "rhocut": True,
                    "pt500cut": True,
                    "HLT_AK8PFJet450": True,
                },
                "fail_trigger": {
                    "n2ddt": False,
                    "rhocut": True,
                    "pt500cut": True,
                    "HLT_AK8PFJet450": True,
                },
            },
            "ttbar": {
                "inclusive": {"pt200cut": True},
                "pass": {"tau32": True, "pt200cut": True},
                "passW": {"tau32": False, "tau21": True, "pt200cut": True},
                "fail": {"tau32": False, "tau21": False, "pt200cut": True},
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
                        f"{selection}_mjet_{region}_jec_up": hist.Hist(
                            mJ_fit_ax,
                            self._pT_fit_ax[selection],
                            dataset_ax,
                            jec_applied_ax,
                            eta_regions_ax,
                            storage=hist.storage.Weight(),
                        ),
                        f"{selection}_mjet_{region}_jec_down": hist.Hist(
                            mJ_fit_ax,
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

    def postprocess(self, accumulator):
        return accumulator

    def process(self, events):
        out = self.accumulator()

        dataset = events.metadata["dataset"]
        selection = events.metadata["selection"]

        isMC = "data" not in dataset.lower()

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

        pt_raw = events.pt
        pt = pt_raw * jecfactor

        ptgen_ = events.pt_gen_ak8

        mjet_raw = events.mjet
        mjet = mjet_raw * jecfactor

        mJgen_ = events.msd_gen_ak8

        rho = 2 * np.log(mjet / pt)
        rho_raw = 2 * np.log(mjet_raw / pt_raw)
        rho_corrected_pt = 2 * np.log(mjet_raw / pt)
        rho_corrected_mJ = 2 * np.log(mjet / pt_raw)

        eta_ = events.eta
        phi_ = events.phi

        out["pt"].fill(dataset=dataset, jecAppliedOn="pt", pt=pt, weight=events.weight)
        out["pt"].fill(
            dataset=dataset, jecAppliedOn="none", pt=pt_raw, weight=events.weight
        )
        out["eta"].fill(dataset=dataset, eta=eta_, weight=events.weight)
        out["phi"].fill(dataset=dataset, phi=phi_, weight=events.weight)
        out["mjet"].fill(
            dataset=dataset, jecAppliedOn="mJ", mJ=mjet, weight=events.weight
        )
        out["mjet"].fill(
            dataset=dataset, jecAppliedOn="none", mJ=mjet_raw, weight=events.weight
        )
        out["rho"].fill(
            dataset=dataset, jecAppliedOn="pt&mJ", rho=rho, weight=events.weight
        )
        out["rho"].fill(
            dataset=dataset,
            jecAppliedOn="pt",
            rho=rho_corrected_pt,
            weight=events.weight,
        )
        out["rho"].fill(
            dataset=dataset,
            jecAppliedOn="mJ",
            rho=rho_corrected_mJ,
            weight=events.weight,
        )
        out["rho"].fill(
            dataset=dataset, jecAppliedOn="none", rho=rho_raw, weight=events.weight
        )

        out["npv"].fill(dataset=dataset, npv=events.n_pv, weight=events.weight)
        if isMC:
            out["ntrueint"].fill(
                dataset=dataset, ntrueint=events.n_trueint, weight=events.weight
            )

        # for jec_applied_on in ['none','pt','pt&mJ']:
        for jec_applied_on in ["pt", "pt&mJ"]:
            selections = PackedSelection()
            pt_ = pt_raw
            if "pt" in jec_applied_on:
                pt_ = pt
            mJ_ = mjet_raw
            if "mJ" in jec_applied_on:
                mJ_ = mjet

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

            # selections.add("rhocut",
            #                (rho_<-2.1)
            #                &(rho_>-6.0))
            selections.add("rhocut", rho_ < -2.1)

            selections.add("tau21", events.tau21 < 0.45)
            selections.add("tau32", events.tau32 < 0.5)

            selections.add(
                "unfolding",
                (events.pass_reco_selection == 1) & (events.pass_gen_selection == 1),
            )

            selections.add("jetpfid", events.jetpfid == 1)

            selection = events.metadata["selection"]

            selections.add(
                # "HLT_AK8PFJet450", self.passes_trigger(events, "HLT_AK8PFJet450_v*")
                "HLT_AK8PFJet450", ak.ones_like(events.pt, dtype=bool)
            )

            for region in self._regions[selection].keys():
                smask_unfolding = selections.require(**self._regions[selection][region], unfolding=True, jetpfid=True)
                corrector = self.mjet_reco_correction[
                    "response_g_" + ("jec" if "mJ" in jec_applied_on else "nojec") + f"_{self._year}"
                ]
                msd_correction = (
                    corrector.evaluate(pt_[smask_unfolding])
                    if selection == "vjets"
                    else 1.0  # ak.ones_like(pt_[smask_unfolding])
                )

                out[f"{selection}_mjet_unfolding_{region}"].fill(
                    ptreco=pt_[smask_unfolding],
                    mJreco=mJ_[smask_unfolding] * msd_correction,
                    ptgen=ptgen_[smask_unfolding],
                    mJgen=mJgen_[smask_unfolding],
                    dataset=dataset,
                    jecAppliedOn=jec_applied_on,
                    weight=events.weight[smask_unfolding],
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
        return out


if __name__ == "__main__":
    workflow = CoffeaWorkflow("JMSTemplates")

    workflow.parser.add_argument("--output", "-o", type=str, default="jms_templates.coffea")
    workflow.parser.add_argument("--year", default="UL17")
    workflow.parser.add_argument("--JEC", default="nominal")
    workflow.parser.add_argument("--maxfiles", type=int, default=-1)

    args = workflow.parse_args()

    workflow.processor_instance = JMSTemplates(args.year, args.JEC)
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

    files = {}
    for selection in ["vjets", "ttbar"]:
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
        workflow.init_dask_htcondor_client(1, 8, 5)

    print("starting coffea runner")

    output = workflow.run(files)

    if not output_file_path.endswith(".coffea"):
        output_file_path += ".coffea"

    save(output, output_file_path)
