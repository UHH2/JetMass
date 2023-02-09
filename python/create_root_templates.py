#!/usr/bin/env pythonJMS.sh
from coffea.util import load, save
import os


def flatten_templates(
        hists, hist_name, samples, jec_applied_on="pt&mJ", legacy_hist_name_pattern="",
        eta_region='inclusive'
):
    hist_axes = hists[hist_name].axes
    templates = {}
    eta_regions = {"barrel": 0, "endcap": 1, "inclusive": sum}

    for pt_bin in range(len(hist_axes["pt"].edges)):
        pt_inclusive = pt_bin == len(hist_axes["pt"].edges) - 1
        if pt_inclusive:
            pt_bin_name = "inclusive"
        else:
            pt_bin_name = (
                f"{hist_axes['pt'].edges[pt_bin]}to{hist_axes['pt'].edges[pt_bin+1]}"
            )
            pt_bin_name = pt_bin_name.replace(".0", "")
            pt_bin_name = pt_bin_name.replace("inf", "Inf")
        for sample in samples:
            # print(sample,hist_axes['dataset'])
            if sample not in hist_axes["dataset"]:
                continue
            if sample == "Data" and "variation" in hist_name:
                continue
            legacy_hist_name = legacy_hist_name_pattern
            legacy_hist_name = legacy_hist_name.replace("vjets_", f"W_{sample}__")
            legacy_hist_name = legacy_hist_name.replace("vjets_", "")
            legacy_hist_name = legacy_hist_name.replace("ttbar_", f"top_{sample}__")
            legacy_hist_name = legacy_hist_name.replace("ttbar_", "")
            # legacy_hist_name = legacy_hist_name.replace('ttbar_',f'top_{sample}__')
            # legacy_hist_name = legacy_hist_name.replace('variation_',f'{pt_bin_name}_')
            legacy_hist_name = legacy_hist_name.replace("PTBIN", pt_bin_name)
            # print(legacy_hist_name)

            hist_id_dict = {
                "pt": sum if pt_inclusive else pt_bin,
                "dataset": sample,
                "jecAppliedOn": jec_applied_on,
            }
            # TODO remove ones variation hists also have eta bins!
            if "abs_eta_regions" in [ax.name for ax in hist_axes]:
                hist_id_dict["abs_eta_regions"] = eta_regions[eta_region]
            h_ = hists[hist_name][hist_id_dict]

            templates[legacy_hist_name] = h_
    return templates


def make_unfolding_templates(
    hists,
    hist_name,
    samples,
    signal_samples,
    jec_applied_on="pt&mJ",
    legacy_hist_name_pattern="",
):
    templates = {}
    hist_axes = hists[hist_name].axes

    pt_gen_edges = hist_axes["ptgen"].edges
    msd_gen_edges = hist_axes["mJgen"].edges

    # templates for fit
    for pt_bin in range(len(hist_axes["ptreco"].edges)):
        pt_inclusive = pt_bin == len(hist_axes["ptreco"].edges) - 1
        if pt_inclusive:
            pt_bin_name = "inclusive"
        else:
            pt_bin_name = f"{hist_axes['ptreco'].edges[pt_bin]}to{hist_axes['ptreco'].edges[pt_bin+1]}"
            pt_bin_name = pt_bin_name.replace(".0", "")
            pt_bin_name = pt_bin_name.replace("inf", "Inf")

        for sample in samples:
            # print(sample,hist_axes['dataset'])
            if sample not in hist_axes["dataset"]:
                continue

            import re

            sample_info_pattern = re.compile(
                "(?P<selection>vjets|ttbar)_(?P<sample>[A-Za-z]*)"
            )
            sample_name = sample_info_pattern.match(sample).groupdict()["sample"]

            legacy_hist_name = legacy_hist_name_pattern
            legacy_hist_name = legacy_hist_name.replace("vjets_", f"W_{sample}__")
            legacy_hist_name = legacy_hist_name.replace("vjets_", "")
            legacy_hist_name = legacy_hist_name.replace("ttbar_", f"top_{sample}__")
            legacy_hist_name = legacy_hist_name.replace("ttbar_", "")
            legacy_hist_name = legacy_hist_name.replace("PTBIN", pt_bin_name)

            hist_id_dict = {
                "ptreco": sum if pt_inclusive else pt_bin,
                "dataset": sample,
                "jecAppliedOn": jec_applied_on,
                "mJgen": sum,
                "ptgen": sum,
            }

            h_ = hists[hist_name][hist_id_dict]

            templates[legacy_hist_name] = h_

            if sample in signal_samples:
                templates[
                    legacy_hist_name.replace(sample_name, f"{sample_name}_onegenbin")
                ] = h_

            for iptgen in range(len(pt_gen_edges) - 1):
                for imsdgen in range(len(msd_gen_edges) - 1):
                    if sample in signal_samples:

                        gen_bin_name = f"ptgen{iptgen}_msdgen{imsdgen}"

                        legacy_hist_name_unfolding = legacy_hist_name_pattern
                        legacy_hist_name_unfolding = legacy_hist_name_unfolding.replace(
                            "vjets_", f"W_{sample}_{gen_bin_name}__"
                        )
                        legacy_hist_name_unfolding = legacy_hist_name_unfolding.replace(
                            "vjets_", ""
                        )
                        legacy_hist_name_unfolding = legacy_hist_name_unfolding.replace(
                            "ttbar_", f"top_{sample}_{gen_bin_name}__"
                        )
                        legacy_hist_name_unfolding = legacy_hist_name_unfolding.replace(
                            "ttbar_", ""
                        )
                        legacy_hist_name_unfolding = legacy_hist_name_unfolding.replace(
                            "PTBIN", pt_bin_name
                        )

                        hist_id_dict = {
                            "ptreco": sum if pt_inclusive else pt_bin,
                            "dataset": sample,
                            "jecAppliedOn": jec_applied_on,
                            "mJgen": imsdgen,
                            "ptgen": iptgen,
                        }
                        h_ = hists[hist_name][hist_id_dict]

                        templates[legacy_hist_name_unfolding] = h_

        # truth templates
        for iptgen in range(len(pt_gen_edges) - 1):
            for sample in samples:
                if sample in signal_samples:
                    hist_id_dict = {
                        "ptreco": sum,
                        "mJreco": sum,
                        "dataset": sample,
                        "jecAppliedOn": jec_applied_on,
                        "ptgen": iptgen,
                    }

                    h_ = hists[hist_name][hist_id_dict]
                    truth_hist_name = hist_name
                    truth_hist_name = truth_hist_name.replace(
                        "unfolding", f"ptgen{iptgen}"
                    )
                    templates[truth_hist_name] = h_

    return templates


if __name__ == "__main__":
    import argparse
    from utils import hist_to_th1
    import uproot

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", default="templates_2017.coffea")
    parser.add_argument("--output", "-o", default="templates_1d_2017")
    parser.add_argument("--JEC", default="pt&mJ", choices=["none", "pt", "pt&mJ"])
    parser.add_argument("--eta", default="inclusive", choices=["inclusive", "barrel", "endcap"])
    parser.add_argument("--unfolding", action="store_true")

    args = parser.parse_args()

    hists = load(args.input)

    selections = {
        "vjets": {
            "regions": ["inclusive", "pass", "fail"],
            "samples": [
                "Data",
                "QCD",
                "WJets",
                "WJetsUnmatched",
                "WJetsMatched",
                "ZJets",
                "ZJetsUnmatched",
                "ZJetsMatched",
                "TTToHadronic",
                "TTToSemiLeptonic",
                "ST_tW_top",
                "ST_tW_antitop",
            ],
            # 'signal_samples':["WJets","WJetsMatched","ZJets"]
            "signal_samples": ["WJetsMatched"],
        },
        "ttbar": {
            "regions": ["inclusive", "pass", "passW", "fail"],
            "samples": [
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
            "signal_samples": ["TTToSemiLeptonic_mergedW"],
        },
    }
    print(f"flattening templates and saving to ROOT from {args.input}")
    templates = {}
    for selection in selections.keys():
        samples = [f"{selection}_{s}" for s in selections[selection]["samples"]]
        signal_samples = [
            f"{selection}_{s}" for s in selections[selection]["signal_samples"]
        ]
        for region in selections[selection]["regions"]:

            if args.unfolding:
                templates.update(
                    make_unfolding_templates(
                        hists,
                        f"{selection}_mjet_unfolding_{region}",
                        samples=samples,
                        signal_samples=signal_samples,
                        jec_applied_on=args.JEC,
                        legacy_hist_name_pattern=f"{selection}_mjet_PTBIN_{region}",
                    )
                )
            else:
                templates.update(
                    flatten_templates(
                        hists,
                        f"{selection}_mjet_{region}",
                        samples=samples,
                        jec_applied_on=args.JEC,
                        legacy_hist_name_pattern=f"{selection}_mjet_PTBIN_{region}",
                        eta_region=args.eta
                    )
                )
                for variation in ["all", "chargedH", "neutralH", "gamma", "other"]:
                    templates.update(
                        flatten_templates(
                            hists,
                            f"{selection}_mjet_0_0_{variation}_variation_{region}__up",
                            samples=samples,
                            jec_applied_on=args.JEC,
                            legacy_hist_name_pattern=f"{selection}_mjet_0_0_{variation}_PTBIN_{region}__up",
                            eta_region=args.eta
                        )
                    )
                    templates.update(
                        flatten_templates(
                            hists,
                            f"{selection}_mjet_0_0_{variation}_variation_{region}__down",
                            samples=samples,
                            jec_applied_on=args.JEC,
                            legacy_hist_name_pattern=f"{selection}_mjet_0_0_{variation}_PTBIN_{region}__down",
                            eta_region=args.eta
                        )
                    )

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    save(templates, f"{args.output}.coffea")
    fout = uproot.recreate(f"{args.output}.root")
    for hist_name, H in templates.items():
        fout[hist_name] = hist_to_th1(H, hist_name)
    fout.close()
