#!/usr/bin/env python3
from __future__ import print_function
import ROOT
import numpy as np
import os
import json
import glob
import tqdm
from CombineUtils import import_config
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
import sys
import gc
mpl.use("Agg")
hep.set_style(hep.style.CMS)
sys.path.append("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/rhalphalib/")


def correlation_coefficient(fname, pois):
    file_ = ROOT.TFile.Open(fname, "READ")
    fit_s = None
    corr_coeff = np.array([-1])
    try:
        fit_s = file_.Get("fit_s")
        n_pars = fit_s.floatParsFinal().getSize()

        poi_ind = [fit_s.floatParsFinal().index(p) for p in pois]

        cov = np.array([[fit_s.correlationMatrix()[row][col] for col in range(n_pars)] for row in range(n_pars)])
        cov_inv = np.linalg.inv(cov)
        corr_coeff = np.sqrt(1 - 1 / (cov_inv[poi_ind, poi_ind] * cov[poi_ind, poi_ind]))
        file_.Close()
        del fit_s

        del cov, cov_inv
    except BaseException as e:
        print(e)
    return corr_coeff


if __name__ == "__main__":
    import argparse
    # from FitSubmitter import FitSubmitter
    from GenericJobSubmitter import JobSubmitter
    rhalph_path = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/"

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", default="")
    parser.add_argument("--year", help="year to be used in scan as basis.", required=True)
    parser.add_argument("-s", "--submit", action="store_true", help="submit jobs for scan.")
    parser.add_argument("--extract", action="store_true", help="extract correlation from scan.")
    parser.add_argument(
        "--plots", action="store_true", help="plot extracted average correlation coefficient from scans."
    )
    parser.add_argument("--debug", action="store_true", help="perform a dry-run without submission.")
    parser.add_argument("--tagger", choices=["substructure", "particlenetDDT"], default="substructure")

    args = parser.parse_args()

    # Uniform, NoWeight, NoWeightLargeCutoff, WidthWeight, WidthWeightLargeCutoff
    # for strategy in ["Uniform", "NoWeight", "NoWeightLargeCutoff", "WidthWeight", "WidthWeightLargeCutoff"]:
    stepsizes = {
        "substructure": {
            "Uniform": 0.05,
            "NoWeight": 0.05,
            "WidthWeight": 0.000005,
            "NoWeightLargeCutoff": 0.05,
            "WidthWeightLargeCutoff": 0.000005,
        },
        "particlenetDDT": {
            "Uniform": 0.05,
            "NoWeight": 0.05,
            "WidthWeight": 0.000005,
            "NoWeightLargeCutoff": 0.05,
            "WidthWeightLargeCutoff": 0.000005,
        },
    }
    starts = {
        "substructure": {
            "Uniform": 0.05,
            "NoWeight": 0.05,
            "WidthWeight": 0.0001,
            "NoWeightLargeCutoff": 0.2,
            "WidthWeightLargeCutoff": 0.0001,
        },
        "particlenetDDT": {
            "Uniform": 0.05,
            "NoWeight": 0.05,
            "WidthWeight": 0.0001,
            "NoWeightLargeCutoff": 0.2,
            "WidthWeightLargeCutoff": 0.0001,
        },
    }
    for strategy in ["Uniform", "NoWeight"]:#, "WidthWeight"]:
    # for strategy in ["NoWeightLargeCutoff", "WidthWeightLargeCutoff"]:

        n_steps = 250
        # stepsize = 0.025
        stepsize = stepsizes[args.tagger][strategy]
        # stepsize = 0.01
        # if "WidthWeight" in strategy:
        #     stepsize = 0.00002
        #     # stepsize = 0.00001
        start = starts[args.tagger][strategy]
        end = stepsize * n_steps + start
        scan_dict = {
            "regularizationStrength": [
                ("delta{}".format("%.8f" % delta).replace(".", "p"), delta)
                for delta in np.arange(0, n_steps) * stepsize + start
            ]
        }

        # prefixe = ["MSD", "PT", "All"]
        prefixe = ["PT", "All"]
        # if strategy == "Uniform" or strategy == "NoWeight":
        #     prefixe = ["PT"]
        # if strategy == "WidthWeight":
        #     prefixe = ["All"]
        base_dir = "DUST/RegStrScans_{}/".format(args.tagger)
        workdir_pattern = base_dir+"{}_RegStrScan/"
        strat_name = strategy+args.year+args.name
        names = [prefix+strat_name for prefix in prefixe]
        for name in names:
            workdir = workdir_pattern.format(name)

            regularization_dict = {
                "MSD"+strat_name: ["msd"],
                "PT"+strat_name: ["pt"],
                "All"+strat_name: ["msd", "pt"],
            }

            if args.submit:
                # config = import_config("configs/unfolding/WJets_{}.py".format(args.year))
                config = {}
                config["regularization"] = regularization_dict[name]
                if "LargeCutoff" in strategy:
                    config["pt_cutoff"] = 3000
                    config["msd_cutoff"] = 300
                else:
                    config["pt_cutoff"] = 1400
                    config["msd_cutoff"] = 100

                if "Uniform" in strategy:
                    config["uniformGenbins"] = "True"
                    config["scaleGenBinWidth"] = "False"
                elif "NoWeight" in strategy:
                    config["uniformGenbins"] = "False"
                    config["scaleGenBinWidth"] = "False"
                elif "WidthWeight" in strategy:
                    config["uniformGenbins"] = "False"
                    config["scaleGenBinWidth"] = "True"

                if name == "":
                    name = os.path.splitext(os.path.basename(args.config))[0]

            if args.submit:
                # submitter = FitSubmitter(
                #     config=config,
                #     workdir=workdir,
                #     dry_run=args.debug,
                # )
                # submitter.extra_jetmass_options += " -M unfolding --skipTemplatePlots "
                # if args.tagger != "substructure":
                #     submitter.extra_jetmass_options += " --tagger {} ".format(args.tagger)
                # submitter.scan(scan_dict, name+"_RegStrScan")
                common_wrapper_lines = [
                    "#!/bin/bash",
                    "source /cvmfs/cms.cern.ch/cmsset_default.sh",
                    "cd {}/CMSSW_11_3_4/src".format(rhalph_path),
                    "eval `scramv1 runtime -sh`",
                    "cd {}".format(rhalph_path),
                ]
                jetmass_arguments = ""
                jetmass_exe = "jetmass.py"
                if args.year == "RunII":
                    jetmass_arguments += " --configs configs/unfolding/WJets_UL*.py "
                    jetmass_arguments += " --name WJetsRunIIUnfolding${1} "
                    jetmass_arguments += ' --extra-options "'
                    jetmass_exe = "jetmass_combination.py"
                else:
                    jetmass_arguments += " configs/unfolding/WJets_{}.py".format(args.year)
                    jetmass_arguments += "  --skipTemplatePlots "
                if args.tagger != "substructure":
                    jetmass_arguments += " --tagger {} ".format(args.tagger)
                # jetmass_arguments += ' --config-update \"{\\"regularizationStrength\\":$1}\" --job_index $2'
                config["regularizationStrength"] = "$2"
                if args.year == "RunII":
                    jetmass_arguments += '" '
                jetmass_arguments += ' --config-update \"{}\"'.format(str(config))
                jetmass_arguments += " --job_index $1 "
                wrapper_lines = common_wrapper_lines + [
                    "python jetmass.py -M unfolding {} --workdir {}".format(jetmass_arguments, workdir)
                ]
                submitter = JobSubmitter(base_dir)
                submitter.add_job_array(
                    name + "RegStrScan",
                    os.path.basename(workdir),
                    common_wrapper_lines
                    + ["python {} -M unfolding {} --workdir {} --correlationCoeffcients".format(jetmass_exe, jetmass_arguments, workdir)],
                    ["$(MyJobIndex)", "$(RegStrPoint)"],
                    "MyJobIndex, RegStrPoint from ("
                    + ("\n".join("%s %f" % point for point in scan_dict["regularizationStrength"])) + "\n)",
                )
                submitter.submit_jobs(debug=args.debug)

            if args.extract:
                all_missing_npys = []
                for fitD_path in glob.glob("{}/*{}*/fitDiagnostics.root".format(os.path.abspath(workdir), args.year)):
                    fit_dir = os.path.dirname(fitD_path)
                    if not os.path.isfile("{}/corr_coeff_new.npy"):
                        all_missing_npys.append(fitD_path)

                parallel = 10
                for chunK_missing_npys in tqdm.tqdm(
                        [all_missing_npys[i:i+parallel] for i in range(0, len(all_missing_npys), parallel)],
                        desc="calculating correlation coefficient in chunks",
                ):
                    recover_cmd = 'printf "%s\\n" {}'.format(" ".join(chunK_missing_npys))
                    recover_cmd += "| xargs -I{} -n1 -P 10 python dump_correlation_coefficients.py -i {}"
                    os.system(recover_cmd)

                corr_coeff_dict = {}
                for npy_path in tqdm.tqdm(
                        glob.glob("{}/*{}*/corr_coeff.npy".format(os.path.abspath(workdir), args.year)),
                        desc="Extracting correlation coefficients"
                        ):
                    fit_dir = os.path.dirname(npy_path)
                    with open("{}/config.json".format(fit_dir), "r") as config_file:
                        config = json.load(config_file)
                        corr_coeff = np.load(npy_path, allow_pickle=True)
                        delta_string = config["ModelName"].split("delta")[-1]
                        delta = float(delta_string.replace("p", "."))
                        corr_coeff_dict[delta] = corr_coeff

                np.save("{}_{}_corr_coeff.npy".format(name, args.tagger), corr_coeff_dict)

        if args.plots:
            plot_dir = "{}_{}_plots".format(strat_name, args.tagger)
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
            y_min = 0.85
            colors = [
                "#1b9e77",
                "#d95f02",
                "#7570b3"
            ]

            coeffs = {
                name: np.load("{}_{}_corr_coeff.npy".format(name, args.tagger), allow_pickle=True, encoding="bytes")[()]
                for name in names
            }

            average_corr_graphs = {}
            max_corr_graphs = {}
            min_deltas = {}
            min_metrics = []
            for name in names:
                deltas = sorted(coeffs[name].keys())
                average_corr = [coeffs[name][delta].sum()/len(coeffs[name][delta]) for delta in deltas]
                max_corr = [coeffs[name][delta].max() for delta in deltas]
                x = np.array(deltas)
                average_corr_graphs[name] = (x, np.array(average_corr))
                max_corr_graphs[name] = (x, np.array(max_corr))

            f, ax = plt.subplots(figsize=(10, 7))

            for iname, name in enumerate(names):
                delta, metric = average_corr_graphs[name]
                # delta, metric = max_corr_graphs[name]
                m = metric > 0
                delta = delta[m]
                metric = metric[m]
                if len(metric) == 0:
                    continue
                min_delta = delta[np.argmin(metric)]
                min_deltas[name] = min_delta
                min_metric = metric.min()
                min_metrics.append(min_metric)
                ax.plot(
                    delta,
                    metric,
                    "-",
                    label=r"%s $\delta_{\mathrm{min}}=%.2e$" % (name.replace(strat_name, ""), min_delta),
                    color=colors[iname],
                )
                ax.plot([min_delta] * 2, [0, min_metric], "--", color=colors[iname])
            if min(min_metrics) < y_min:
                y_min = 0.8 * min(min_metrics)
            ax.set_ylim(0.7, 1)
            hep.label._exp_label(llabel="Private work (CMS data/simulation)", fontsize=20, ax=ax, year=args.year)
            ax.set_xlabel(r"$\delta$")
            ax.set_ylabel("average correlation coeff.")
            ax.legend(fontsize=18)
            f.savefig("{}/{}_scan.pdf".format(plot_dir, strat_name), bbox_inches="tight")

            for name in names:
                workdir = workdir_pattern.format(name)
                fit_d_paths = glob.glob("{}/*/fitDiagnostics.root".format(workdir))
                if name not in min_deltas:
                    continue
                delta_str = ("%.8f" % (min_deltas[name])).replace(".", "p")
                best_fit = [d for d in fit_d_paths if delta_str in d][0]

                cmd = "./plot_pois.py -f {} -r -1.0 5.0".format(best_fit)
                os.system(cmd)
                cmd = "cp {} {}/{}_delta{}_pois.pdf".format(
                    best_fit.replace(".root", "_pois.pdf"), plot_dir, name, delta_str
                )
                os.system(cmd)
