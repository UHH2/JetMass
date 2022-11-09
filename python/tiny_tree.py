#!/usr/bin/env coffea_python.sh
import uproot
import coffea.lookup_tools
from coffea.nanoevents.methods import vector
import numpy as np
import awkward as ak
import glob


def load_tree(
        dirname: str,
        fname_pattern: str,
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
):

    if selection == "vjets":
        ddtmaps_path = (
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/"
            "JetMass/Histograms/ddtmaps.npy"
        )

        n2ddtmap = np.load(ddtmaps_path, allow_pickle=True).item()
        year = 'UL17'
        corrected_str = '_corrected_pt_mass'

        n2ddt_LOT = coffea.lookup_tools.dense_lookup.dense_lookup(
            n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str}'][0],
            dims=(n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str}'][1],
                  n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str}'][2]))

        def n2ddt(e):
            return (e.N2 - n2ddt_LOT(e.rho, e.pt)) < 0

        events_sel = events[
            (events.rho > -6.0)
            & (events.rho < -2.1)
            & (n2ddt(events))
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

    preselection_tree = load_tree(
        dirname="/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/workdir_vjets_UL17/",
        fname_pattern="*WJetsToQQ*.root",
    )
    selection_tree = apply_selection(preselection_tree, "vjets")
    dump_tiny_tree("WJetsToQQ_tinyTree.parquet", selection_tree)
