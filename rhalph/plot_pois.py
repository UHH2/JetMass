#!/usr/bin/env python3
from __future__ import print_function
import ROOT
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
hep.set_style(hep.style.CMS)


def extract_pois(fname, pois_pattern, fit_result_name="fit_s"):
    file_ = ROOT.TFile.Open(fname, "READ")
    try:
        fit_result = file_.Get(fit_result_name)
    except BaseException as e:
        print(e)
        print("could not find fit result ({}) in provided file.".format(fit_result))
        exit(-1)

    poi_result_tmp = {}
    pois = []
    for p in fit_result.floatParsFinal():
        if pois_pattern in p.GetName():
            poi_result_tmp[p.GetName()] = [p.getVal(), p.getErrorHi(), p.getErrorLo()]
            pois.append(p.GetName())
    pois.sort()
    return [np.array([poi_result_tmp[name][i] for name in pois]) for i in [0, 1, 2]], pois


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fitDiagnostics", help="fitDiagnostics file containing fitResults.")
    parser.add_argument("-P", "--pois", default="r_")
    parser.add_argument("-e", "--expectSignal", default=1.0, type=float)

    args = parser.parse_args()

    fit_dir = os.path.dirname(os.path.abspath(args.fitDiagnostics))

    # prefit_pois = extract_pois(args.fitDiagnostics, args.pois, fit_result_name="nuisances_prefit_res")
    postfit_pois, pois = extract_pois(args.fitDiagnostics, args.pois, fit_result_name="fit_s")

    f, ax = plt.subplots(figsize=(20, 10))

    xtick_ids = np.arange(len(pois))
    xtick_labels = np.array(pois)
    # ax.plot(xtick_ids, prefit_pois[0], "k--", alpha=0.7, label="prefit")
    # ax.plot(xtick_ids, prefit_pois[0], "k--", alpha=0.7, label="prefit")
    ax.errorbar(
        xtick_ids,
        postfit_pois[0],
        # yerr=[postfit_pois[0] - postfit_pois[2], postfit_pois[0] + postfit_pois[1]],
        yerr=[-postfit_pois[2], postfit_pois[1]],
        label="postfit",
        fmt=".",
    )
    ax.plot(ax.get_xlim(), [args.expectSignal]*2, "k--", label="expectSignal", alpha=0.7)

    plt.xticks(xtick_ids, xtick_labels, rotation=20)
    ax.legend()
    ax.set_ylim(-11,11)
    f.savefig(f"{args.fitDiagnostics.replace('.root','')}_pois.pdf", bbox_inches="tight")
