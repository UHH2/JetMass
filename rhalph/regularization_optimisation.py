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
mpl.use("Agg")
import mplhep as hep
hep.set_style(hep.style.CMS)
import sys
sys.path.append("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/rhalphalib/")


def correlation_coefficient(fname, pois):
    file_ = ROOT.TFile.Open(fname, "READ")
    fit_s = file_.Get("fit_s")

    try:
        n_pars = fit_s.floatParsFinal().getSize()

        poi_ind = [fit_s.floatParsFinal().index(p) for p in pois]

        cov = np.array([[fit_s.correlationMatrix()[row][col] for col in range(n_pars)] for row in range(n_pars)])
        cov_inv = np.linalg.inv(cov)
        return np.sqrt(1 - 1 / (cov_inv[poi_ind, poi_ind] * cov[poi_ind, poi_ind]))
    except BaseException:
        return np.array([-1])


if __name__ == "__main__":
    import argparse
    from FitSubmitter import FitSubmitter

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fitDiagnostics", help="fitDiagnostics file containing fit results.")
    parser.add_argument("-n", "--name", default="")
    parser.add_argument("-c", "--config", help="config to be used in scan as basis.")
    parser.add_argument("-s", "--submit", action="store_true", help="submit jobs for scan.")
    parser.add_argument("--extract", action="store_true", help="extract correlation from scan.")
    parser.add_argument(
        "--plots", action="store_true", help="plot extracted average correlation coefficient from scans."
    )
    parser.add_argument("--debug", action="store_true", help="perform a dry-run without submission.")

    args = parser.parse_args()

    n_steps = 75
    stepsize = 0.025
    # stepsize = 0.00002
    start = 0.0
    end = stepsize * n_steps + start
    scan_dict = {
        "regularizationStrength": [
            ("delta{}".format("%.5f" % delta).replace(".", "p"), delta)
            for delta in np.arange(0, n_steps) * stepsize + start
        ]
    }

    prefixe = ["MSD", "PT", "All"]
    workdir_pattern = "DUST/{}_RegStrScan/"
    names = [prefix+args.name for prefix in prefixe]
    for name in names:
        workdir = workdir_pattern.format(name)

        regularization_dict = {
            "MSD"+args.name: ["msd"],
            "PT"+args.name: ["pt"],
            "All"+args.name: ["msd", "pt"],
        }

        if args.submit:
            if args.config == "":
                raise AttributeError("You did not provide a config")
            config = import_config(args.config)
            config["regularization"] = regularization_dict[name]
            if name == "":
                name = os.path.splitext(os.path.basename(args.config))[0]

        workdir = "DUST/{}_RegStrScan/".format(name)

        if args.submit:
            submitter = FitSubmitter(
                config=config,
                workdir=workdir,
                dry_run=args.debug,
            )
            submitter.extra_jetmass_options += " -M unfolding "
            submitter.scan(scan_dict, name+"_RegStrScan")

        if args.extract:
            corr_coeff_dict = {}
            for config_path in tqdm.tqdm(
                    glob.glob("{}/*.json".format(os.path.abspath(workdir))),
                    desc="Extracting correlation coefficients"
                    ):
                fitD_path = os.path.abspath(config_path.replace(".json", "/fitDiagnostics.root"))
                fit_dir = os.path.dirname(fitD_path)
                config = json.load(open("{}/config.json".format(fit_dir), "r"))
                pois = ["r_{}".format(genbin) for genbin in config.get("genbins", [])]
                corr_coeff = correlation_coefficient(fitD_path, pois)
                delta_string = config["ModelName"].split("delta")[-1]
                delta = float(delta_string.replace("p", "."))
                corr_coeff_dict[delta] = corr_coeff

            np.save("{}_corr_coeff.npy".format(name), corr_coeff_dict)

    if args.plots:
        plot_dir = "{}_plots".format(args.name)
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        y_min = 0.85
        colors = [
            "#1b9e77",
            "#d95f02",
            "#7570b3"
        ]
        coeffs = {
            name: np.load("{}_corr_coeff.npy".format(name), allow_pickle=True, encoding="bytes")[()] for name in names
        }
        average_corr_graphs = {}
        max_corr_graphs = {}
        min_deltas = {}
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
            m = metric > 0
            delta = delta[m]
            metric = metric[m]
            min_delta = delta[np.argmin(metric)]
            min_deltas[name] = min_delta
            min_metric = metric.min()
            ax.plot(
                delta,
                metric,
                "-",
                label=r"%s $\delta_{\mathrm{min}}=%.2e$" % (name.replace(args.name,""), min_delta),
                color=colors[iname],
            )
            ax.plot([min_delta] * 2, [y_min, min_metric], "--", color=colors[iname])

        ax.set_ylim(y_min, 1)
        hep.label._exp_label(llabel="Private work (CMS data/simulation)", fontsize=20, ax=ax, year="UL17")
        ax.set_xlabel(r"$\delta$")
        ax.set_ylabel("average correlation coeff.")
        ax.legend(fontsize=18)
        f.savefig("{}/{}_scan.pdf".format(plot_dir, args.name), bbox_inches="tight")

        for name in names:
            workdir = workdir_pattern.format(name)
            fit_d_paths = glob.glob("{}/*/fitDiagnostics.root".format(workdir))
            delta_str = ("%.5f" % (min_deltas[name])).replace(".", "p")
            best_fit = [d for d in fit_d_paths if delta_str in d][0]

            cmd = "./plot_pois.py -f {}".format(best_fit)
            os.system(cmd)
            cmd = "cp {} {}/{}_delta{}_pois.pdf".format(
                best_fit.replace(".root", "_pois.pdf"), plot_dir, name, delta_str
            )
            os.system(cmd)
