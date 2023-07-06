#!/usr/bin/env pythonJMS.sh
import uproot
import coffea.lookup_tools
import correctionlib
from coffea.nanoevents.methods import vector
import numpy as np
import awkward as ak
import glob


def load_tree(
        dirname: str,
        fname_pattern: str,
        year: str,
):
    tree_name = "AnalysisTree"

    tree_paths = [
        f"{fpath}:{tree_name}"
        for fpath in glob.glob(dirname + fname_pattern)
        if uproot.open(fpath)[tree_name].num_entries > 0
    ]

    def filter_(x):
        return x in [
            "pt",
            "eta",
            "phi",
            "mass",
            "jecfactor",
            "jecfactor_SD",
            "pt_1",
            "eta_1",
            "phi_1",
            "mass_1",
            "jecfactor_1",
            "mjet",
            "pt_gen_ak8",
            "msd_gen_ak8",
            "weight",
            "N2",
            "IsMergedWZ",
            "trigger_bits",
            "jetpfid",
            "ParticleNetMD_probXbb",
            "ParticleNetMD_probXcc",
            "ParticleNetMD_probXqq",
            "ParticleNetMD_probQCD",
            "ParticleNetMDDiscriminators_XbbvsQCD",
            "ParticleNetMDDiscriminators_XccvsQCD",
            "ParticleNetMDDiscriminators_XqqvsQCD",
        ]

    events = uproot.concatenate(tree_paths, filter_name=filter_)

    trigger_scalefactors = correctionlib.CorrectionSet.from_file(
        "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/notebooks/data/"
        + "HLT_AK8PFJet_MC_trigger_sf_c2e731345f.json"
    )
    trigger_sf_evaluator_450 = trigger_scalefactors[f"HLT_AK8PFJet450_triggersf_{year}"]
    trigger_sf_evaluator_500 = trigger_scalefactors[f"HLT_AK8PFJet500_triggersf_{year}"]

    pt_corr = events.pt * events.jecfactor
    first_ptbin = pt_corr < 650.0

    events["triggersf"] = ak.where(
        first_ptbin,
        trigger_sf_evaluator_450.evaluate(pt_corr, "nominal"),
        trigger_sf_evaluator_500.evaluate(pt_corr, "nominal"),
    )

    events['rho'] = 2*np.log(events.mjet/events.pt)

    events["ParticleNetMD_probXQQ"] = (
        events["ParticleNetMD_probXqq"] + events["ParticleNetMD_probXcc"] + events["ParticleNetMD_probXbb"]
    )
    events["ParticleNetMDDiscriminators_XbbvsQCD"] = events["ParticleNetMD_probXQQ"] / (
        events["ParticleNetMD_probXQQ"] + events["ParticleNetMD_probQCD"]
    )

    leading_jet = ak.zip(
        {
            "pt": events.pt,
            "eta": events.eta,
            "phi": events.phi,
            "mass": events.mass,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )
    subleading_jet = ak.zip(
        {
            "pt": events.pt_1,
            "eta": events.eta_1,
            "phi": events.phi_1,
            "mass": events.mass_1,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )
    events["Jets"] = ak.concatenate((leading_jet[:, np.newaxis], subleading_jet[:, np.newaxis]), axis=1)
    events["jecfactor"] = ak.concatenate((events.jecfactor[:, np.newaxis], events.jecfactor_1[:, np.newaxis]), axis=1)
    events["ht"] = events.pt + events.pt_1

    return events


def apply_selection(
        events: ak.Array,
        selection: str = "vjets",
        year: str = "UL17",
        tagger: str = "n2ddt",
):
    jetmass_path = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass"
    ddtmaps_n2_path = f"{jetmass_path}/Histograms/ddtmaps_n2.npy"
    ddtmaps_particlenet_path = f"{jetmass_path}/Histograms/ddtmaps_particlenet.npy"

    if selection == "vjets":

        n2ddtmap = np.load(ddtmaps_n2_path, allow_pickle=True).item()
        corrected_str = '_corrected_pt_mass'

        n2ddt_LUT = coffea.lookup_tools.dense_lookup.dense_lookup(
            n2ddtmap[f'discddt_map_{year}_smooth_4_0p05{corrected_str}'][0],
            dims=(n2ddtmap[f'discddt_map_{year}_smooth_4_0p05{corrected_str}'][1],
                  n2ddtmap[f'discddt_map_{year}_smooth_4_0p05{corrected_str}'][2]
                  )
        )

        def n2ddt(e):
            return (e.N2 - n2ddt_LUT(e.rho, e.pt)) < 0

        particlenetddtmap = np.load(ddtmaps_particlenet_path, allow_pickle=True).item()
        pNetMDWvsQCDddt_LUT = coffea.lookup_tools.dense_lookup.dense_lookup(
            particlenetddtmap[
                f"discddt_map_{year}_smooth_4_0p975{corrected_str}"
            ][0],
            dims=(
                particlenetddtmap[
                    f"discddt_map_{year}_smooth_4_0p975{corrected_str}"
                ][1],
                particlenetddtmap[
                    f"discddt_map_{year}_smooth_4_0p975{corrected_str}"
                ][2],
            )
        )

        def tagger_ddt(e):
            if tagger == "n2ddt":
                # 5% QCD eff. -> (tagger score - 5th percentile ) < 0 )
                return (e.N2 - n2ddt_LUT(e.rho, e.pt)) < 0
            elif tagger == "pNetddt":
                # for pNet tagger things have to reversed (5% QCD eff. -> ( 95th percentile - tagger score) < 0 )
                return (pNetMDWvsQCDddt_LUT(e.rho, e.pt) - e["ParticleNetMDDiscriminators_XbbvsQCD"]) < 0

        events_sel = events[
            # (events.rho > -6.0)
            (events.rho < -2.1)
            & (tagger_ddt(events))
            & (events.IsMergedWZ == 1)
            & (events.trigger_bits[:, 7] == 1)
            & (events.jetpfid == 1)
        ]

    return events_sel


def dump_tiny_tree(output_name: str, events: ak.Array):
    events_sel_tosave = events[
        [
            x
            for x in ak.fields(events)
            if x
            not in [
                "pt",
                "eta",
                "phi",
                "mass",
                "pt_1",
                "eta_1",
                "phi_1",
                "mass_1",
                "IsMergedWZ",
                "N2",
                "trigger_bits",
                "jetpfid",
                "jecfactor_1" "ParticleNetMDDiscriminators_XbbvsQCD",
                "ParticleNetMDDiscriminators_XccvsQCD",
                "ParticleNetMDDiscriminators_XqqvsQCD",
                "ParticleNetMDDiscriminators_XQQvsQCD",
            ]
        ]
    ]

    ak.to_parquet(events_sel_tosave, output_name)
    return events_sel_tosave


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--year", default="UL17")
    parser.add_argument("--tagger", default="n2ddt", choices=["n2ddt", "pNetddt"])

    args = parser.parse_args()

    preselection_tree = load_tree(
        dirname=f"/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/workdir_vjets_{args.year}/",
        fname_pattern="*WJetsToQQ*.root",
        year=args.year,
    )
    selection_tree = apply_selection(preselection_tree, "vjets", year=args.year, tagger=args.tagger)
    dump_tiny_tree(f"WJetsToQQ_tinyTree_{args.year}_{args.tagger}.parquet", selection_tree)
