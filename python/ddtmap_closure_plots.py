#!/usr/bin/env pythonJMS.sh
from coffea.util import load
import matplotlib.pyplot as plt
import hist
import numpy as np
import os
from hist.intervals import ratio_uncertainty
from hist.plot import plot_ratio_array
import mplhep as hep
hep.style.use("CMS")

dirname = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/coffea_hists/"


def plot_closure_pt_bin(
    ipt, qcd_pass, qcd_fail, rebin_factor=1, denisty=True, out_dir="ddt_closure", year="", logY=False, sample="QCD"
):
    pt_ax = qcd_pass.axes[1]
    pt_strs = tuple(map(str, pt_ax[ipt]))
    pt_bin_str = f"{pt_strs[0]}to{pt_strs[1]}" if ipt >= 0 else "inclusive"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    pt_tex = r"$%s~\mathrm{GeV} \leq p_{T} < %s~\mathrm{GeV}$" % pt_strs if ipt >= 0 else "inclusive"
    pt_tex = pt_tex.replace(r" < inf~\mathrm{GeV}", "")
    ipt = ipt if ipt >= 0 else sum

    fig = plt.figure(figsize=(11, 9))
    grid = fig.add_gridspec(2, 1, hspace=0, height_ratios=[1, 1 - 9 / 12])
    ax = fig.add_subplot(grid[0])
    rax = fig.add_subplot(grid[1], sharex=ax)
    plt.setp(ax.get_xticklabels(), visible=False)

    phist = qcd_pass[:, ipt][hist.rebin(rebin_factor)]
    fhist = qcd_fail[:, ipt][hist.rebin(rebin_factor)]

    hep.histplot(phist, ax=ax, label="%s pass-region" % sample, color="#1b9e77", density=denisty)
    hep.histplot(fhist, ax=ax, label="%s fail-region" % sample, color="#d95f02", density=denisty)
    ax.set_ylabel(r"$\Delta N / N$")
    y_half = ax.get_ylim()[1] / 2
    if logY:
        ax.set_yscale("log")
        ymin, ymax = np.log10(ax.get_ylim())
        y_half = 10 ** ((ymin + ymax) / 2)

    ax.text(ax.get_xlim()[1] * 0.5, y_half, pt_tex, fontsize=20)
    # hep.cms.label(label=', Work in Progress',year=year,ax=ax,fontsize=20)
    hep.label.exp_label(
        llabel="Private Work (CMS %s)" % ("data/simulation" if sample == "Data" else "simulation"),
        year=year,
        ax=ax,
        fontsize=20,
    )
    ax.legend()
    # ax.plot([50,50],ax.get_ylim(),'k--')
    sumw_num = phist.values()
    sumw_denom = fhist.values()
    rsumw = np.nan_to_num(sumw_num / sumw_denom)
    rsumw_err = np.nan_to_num(ratio_uncertainty(sumw_num, sumw_denom, uncertainty_type="efficiency"))

    plot_ratio_array(phist, rsumw, rsumw_err, rax)

    rax.set_ylim(0.035, 0.065)
    rax.set_ylabel("pass/fail", loc="center")
    fig.align_ylabels()
    # plotratio(phist,fhist,rax)
    # rax.step(centers,rsumw,'k')
    xlim = rax.get_xlim()
    rax.plot(xlim, [0.05, 0.05], "k--", alpha=0.8)
    rax.plot(xlim, [0.045, 0.045], "k--", alpha=0.6)
    rax.plot(xlim, [0.055, 0.055], "k--", alpha=0.6)
    fig.savefig(f"{out_dir}/ddt_closure_{pt_bin_str}_{year}.pdf", bbox_inches="tight")

if __name__ == "__main__":
    for suffix in ["", "_particlenetDDT", "_particlenet"]:
        for sample in ["QCD", "Data"]:
            for year in ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]:
                fname = f"templates_{year}{suffix}.coffea"
                hists = load(dirname + ("DPNote_06-07-23/" if suffix == "_particlenet" else "") + fname)
                hist_name_prefix = "vjets_mjet_"
                # hist_name_prefix = 'vjets_mjet_'
                hists_pass = hists[hist_name_prefix + "pass"]
                hists_fail = hists[hist_name_prefix + "fail"]
                bin_identifiers = {
                    # 'mJgen':sum,'ptgen':sum,
                    "mJ": slice(30, 300),
                    "dataset": f"vjets_{sample}",
                    "jecAppliedOn": "pt&mJ",
                    "abs_eta_regions": sum,
                }
                qcd_pass = hists_pass[bin_identifiers]
                qcd_fail = hists_fail[bin_identifiers]

                for ipt in range(-1, len(hists_pass.axes[1])):
                    plot_closure_pt_bin(ipt, qcd_pass, qcd_fail, rebin_factor=5, year=year, out_dir=f"ddt_closure_{sample}{suffix}", sample=sample)
