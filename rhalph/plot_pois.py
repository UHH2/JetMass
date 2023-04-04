#!/usr/bin/env python3
from __future__ import print_function
import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
hep.set_style(hep.style.CMS)


def extract_pois(fname, pois, fit_result_name="fit_s"):
    file_ = ROOT.TFile.Open(fname, "READ")
    try:
        fit_result = file_.Get(fit_result_name)
    except BaseException as e:
        print(e)
        print("could not find fit result ({}) in provided file.".format(fit_result))
        exit(-1)

    poi_result_tmp = {}
    for p in fit_result.floatParsFinal():
        if p.GetName() in pois or len(pois) == 0:
            poi_result_tmp[p.GetName()] = [p.getVal(), p.getErrorHi(), p.getErrorLo()]

    return [np.array([poi_result_tmp[name][i] for name in pois]) for i in [0, 1, 2]]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fitDiagnostics", help="fitDiagnostics file containing fitResults.")
    parser.add_argument("-P", "--pois", nargs="+", default=[])
    parser.add_argument("-e", "--expectSignal", default=1.0, type=float)

    args = parser.parse_args()

    fit_dir = os.path.dirname(os.path.abspath(args.fitDiagnostics))

    prefit_pois = extract_pois(args.fitDiagnostics, args.pois, fit_result="nuisances_prefit_res")
    postfit_pois = extract_pois(args.fitDiagnostics, args.pois, fit_result="fit_s")

    f, ax = plt.subplots(figsize=(20, 10))

    xtick_ids = np.arange(len(args.pois))
    xtick_labels = np.array(args.pois)
    ax.plot(xtick_ids, prefit_pois[0], "k--", alpha=0.7, label="prefit")
    ax.errorbar(
        xtick_ids,
        postfit_pois[0],
        yerr=[postfit_pois[0] - postfit_pois[2], postfit_pois[0] + postfit_pois[1]],
        label="postfit",
        fmt=".",
    )

    plt.xticks(xtick_ids, xtick_labels, rotation=20)
    ax.legend()
    f.savefig(f"{args.fitDiagnostics.replace('.root','')}_pois.pdf", bbox_inches="tight")
