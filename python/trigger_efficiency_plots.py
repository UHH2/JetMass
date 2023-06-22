#!/usr/bin/env pythonJMS.sh

import uproot
import hist
from hist.intervals import ratio_uncertainty
import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import mplhep as hep
import os
from copy import deepcopy
hep.style.use("CMS")


def setup_axes():
    f_eff, ax_eff = plt.subplots(figsize=(9, 9))
    fig, ax = plt.subplots(figsize=(9, 9))
    return (f_eff, ax_eff, fig, ax)


def finalize_axes(axes, sample, year, outdir, suffix_str="", xlim=(200, 1000), data=True):
    f_eff, ax_eff, fig, ax = axes
    xlabel = "$p_{T,AK8}$"
    ax_eff.axhline(1.0, color="black", linestyle="dashed", linewidth=1.0)
    fig.align_ylabels()
    # hep.label.exp_label(
    #     llabel="Private Work (CMS %s)" % ("data"),
    #     year=year,
    #     ax=ax,
    #     fontsize=20,
    # )
    ax_eff.legend(loc="upper right", fontsize=14, ncol=2)
    ax_eff.set_ylabel(r"$\varepsilon$", loc="center")
    ax_eff.set_ylim(-0.01, 1.5)

    ax.legend(loc="upper right", fontsize=14)
    ax.set_ylabel("Events")
    ax.set_yscale("log")

    cms_label = "Private Work (CMS data/simulation)" if data else "Private Work (CMS simulation)"

    for ax_ in [ax, ax_eff]:
        hep.label.exp_label(llabel=cms_label, year=year, ax=ax_, fontsize=20)
        ax_.set_xlabel(xlabel)
        ax_.set_xlim(*xlim)

    f_eff.savefig(f"{outdir}/trigger_curve_{sample}_{year}{suffix_str}.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/pt_trigger_comparison_{sample}_{year}{suffix_str}.pdf", bbox_inches="tight")


def get_eff_hists(h, varname, ref_trig, probe_trig, **kwargs):
    rebin_factor = kwargs.pop("rebin_factor", 1.0)
    lookup_str = kwargs.pop("lookup_str", "")
    prescale_str = kwargs.pop("prescale_str", "")
    variation = kwargs.pop("variation", "nominal")
    h_ref = h[f"HLTEffHists/{varname}_{ref_trig}_{probe_trig}{lookup_str}_denom"].to_hist()[hist.rebin(rebin_factor)]
    h_probe = h[f"HLTEffHists/{varname}_{ref_trig}_{probe_trig}{lookup_str}{prescale_str}_num"].to_hist()[
        hist.rebin(rebin_factor)
    ]
    if variation != "nominal":
        direction = -1.0 if variation == "down" else 1.0
        h_ref = h_ref + direction * np.sqrt(h_ref.variances())
        h_probe = h_probe + direction * np.sqrt(h_probe.variances())

    sumw_num = h_probe.values()
    sumw_denom = h_ref.values()
    rsumw = np.nan_to_num(np.divide(sumw_num, sumw_denom, out=np.zeros_like(sumw_num), where=sumw_denom != 0))
    try:
        rsumw_err = np.nan_to_num(ratio_uncertainty(sumw_num, sumw_denom, "efficiency"))
    except ValueError as e:
        rsumw_err = 0.0
        print("denominator is larger than numerator.")
        print("Probably the denominator was filled with prescales as weight")
        print("setting errors to zero.")
        print(e)

    return h_ref, h_probe, rsumw, rsumw_err


def efficiency_scalefactor(hists, year, ref_trig, probe_trig, **kwargs):
    h_data = hists["Data"][year]
    h_mc = hists["QCD"][year]
    f, ax = plt.subplots(figsize=(9, 9))
    varname = "AK8_PT"
    xlabel = "$p_{T,AK8}$"
    rebin_factor = kwargs.pop("rebin", 1)
    outdir = kwargs.pop("outdir", ".")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    lookup_str = "_lookup" if kwargs.pop("lookup_trigger", False) else ""
    prescale_str = "_prescale" if kwargs.pop("use_prescales", False) else ""

    def fit_func(z, a, b, c, d):
        return a * (b + 0.5 * (1 - b) * (1 + erf((z - c) / d)))

    def fit_erf(x_, y_, fit_range=(None, None)):
        x = deepcopy(x_)
        y = deepcopy(y_)
        p0 = (1.0, 0.0, 500.0, 40.0)
        if fit_range[0] is not None:
            y = y[x >= fit_range[0]]
            x = x[x >= fit_range[0]]
        if fit_range[1] is not None:
            y = y[x < fit_range[1]]
            x = x[x < fit_range[1]]
        popt, pcov = curve_fit(fit_func, x, y, maxfev=8000, p0=p0)  # , bounds=([0.99,0,300,0],[1,0.04,900,200]))

    x_min = 400
    x_max = 750

    variations = ["nominal", "up", "down"]
    eff_data = {
        var: get_eff_hists(
            h_data,
            varname,
            ref_trig,
            probe_trig,
            rebin_factor=rebin_factor,
            lookup_str=lookup_str,
            prescale_str=prescale_str,
            variation=var,
        )
        for var in variations
    }
    eff_mc = {
        var: get_eff_hists(
            h_mc,
            varname,
            ref_trig,
            probe_trig,
            rebin_factor=rebin_factor,
            lookup_str=lookup_str,
            prescale_str=prescale_str,
            variation=var,
        )
        for var in variations
    }

    eff_fits = {k: {} for k in variations}
    sf_curves = {}
    x_plot = np.linspace(x_min, x_max, 1000)

    for var in variations:
        alpha = 1.0 if var == "nominal" else 0.2
        ax.errorbar(
            eff_data[var][0].axes[0].centers,
            eff_data[var][2],
            yerr=eff_data[var][3],
            label=f"data ({var})",
            color="tab:red",
            linestyle="",
            marker=".",
            alpha=alpha,
        )
        ax.errorbar(
            eff_mc[var][0].axes[0].centers,
            eff_mc[var][2],
            yerr=eff_mc[var][3],
            label=f"mc ({var})",
            color="tab:blue",
            linestyle="",
            marker="x",
            alpha=alpha,
        )

        popt_data, _ = fit_erf(eff_data[var][0].axes[0].centers, eff_data[var][2], fit_range=(100, 1000))
        popt_mc, _ = fit_erf(eff_mc[var][0].axes[0].centers, eff_mc[var][2], fit_range=(100, 1000))
        fit_data = fit_func(x_plot, *popt_data)
        fit_mc = fit_func(x_plot, *popt_mc)

        eff_fits[var]["data"] = fit_data
        eff_fits[var]["mc"] = fit_mc
        sf_curves[var] = fit_data / fit_mc

    ax.fill_between(x_plot, eff_fits["down"]["data"], eff_fits["up"]["data"], alpha=0.4, color="tab:red")
    ax.fill_between(x_plot, eff_fits["down"]["mc"], eff_fits["up"]["mc"], alpha=0.4, color="tab:blue")
    ax.plot(x_plot, eff_fits["nominal"]["data"], label="data fit", color="tab:red", alpha=0.9)
    ax.plot(x_plot, eff_fits["nominal"]["mc"], label="mc fit", color="tab:blue", alpha=0.9)

    ax.fill_between(
        x_plot[x_plot < 500], sf_curves["down"][x_plot < 500], sf_curves["up"][x_plot < 500], alpha=0.1, color="k"
    )
    ax.fill_between(
        x_plot[x_plot >= 500], sf_curves["down"][x_plot >= 500], sf_curves["up"][x_plot >= 500], alpha=0.4, color="k"
    )
    ax.plot(x_plot[x_plot < 500], sf_curves["nominal"][x_plot < 500], alpha=0.2, lw=2, ls="-.", color="k")
    ax.plot(
        x_plot[x_plot >= 500],
        sf_curves["nominal"][x_plot >= 500],
        label="Data/MC SF",
        alpha=0.9,
        lw=2,
        ls="-.",
        color="k",
    )

    ax.axhline(1.0, color="black", linestyle="dashed", linewidth=1.0)
    ax.legend(ncols=1, fontsize=15)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\varepsilon$", loc="center")
    f.savefig(f"{outdir}/data_mc_scalefactor_{probe_trig}_{year}.pdf", bbox_inches="tight")


def efficiency_curve(hists, year, sample, triggers, use_prescales=True, **kwargs):
    colors = kwargs.pop(
        "colors",
        [
            "#66c2a5",
            "#fc8d62",
            "#8da0cb",
            "#e78ac3",
            "#a6d854",
        ],
    )
    axes = kwargs.pop("axes", None)
    if axes is None:
        f_eff, ax_eff, fig, ax = setup_axes()
    else:
        f_eff, ax_eff, fig, ax = axes
    icolor = 0
    varname = "AK8_PT"
    rebin_factor = kwargs.pop("rebin", 1)
    outdir = kwargs.pop("outdir", ".")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    lookup_str = ("_lookup" if kwargs.pop("lookup_trigger", False) else "")
    prescale_str = "_prescale" if use_prescales else ""
    for ref_trig, probe_triggers in triggers.items():
        alpha = 0.8

        for probe_trig in probe_triggers:
            h_ref, h_probe, rsumw, rsumw_err = get_eff_hists(
                hists,
                varname,
                ref_trig,
                probe_trig,
                rebin_factor=rebin_factor,
                lookup_str=lookup_str,
                prescale_str=prescale_str,
            )
            label_suffix = kwargs.get("label_suffix", "")
            hep.histplot(h_ref, ax=ax, color="k", label=ref_trig, linestyle="--")
            hep.histplot(h_probe, ax=ax, color=colors[icolor], alpha=alpha, label=f"{probe_trig}&&{ref_trig}")
            if isinstance(rsumw_err, float) and rsumw_err == 0:
                label_suffix += " (errors=0)/"

            ax_eff.errorbar(
                h_probe.axes[0].centers,
                rsumw,
                yerr=rsumw_err,
                color=colors[icolor],
                alpha=alpha,
                label=r"$\frac{%s&&%s}{%s}$" % (probe_trig, ref_trig, ref_trig) + " " + label_suffix,
            )
            icolor += 1

    if axes is None:
        finalize_axes(
            (f_eff, ax_eff, fig, ax),
            sample,
            year,
            outdir,
            suffix_str=prescale_str + lookup_str,
            data=(sample == "Data"),
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default="/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/ttbarTrees/ForTriggerEff")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
    # years = ["UL16postVFP", "UL17"]
    samples = ["QCD", "Data"]
    triggers_year_comparison = {
        "IsoMuonReference": [
            "AK8PFJet450",
        ],
    }
    triggers = {
        "IsoMuonReference": [
            # "AK8PFJET400",
            "AK8PFJet450",
            "AK8PFJet500",
            "AK8PFJet550",
            # "AK8PFHT900",
            "PFHT1050",
        ],
        # "AK8PFJET320":[
        # "AK8PFJET450",
        # "AK8PFJET500",
        # "AK8PFJET550"
        # ],
        # "AK8PFJET450":["AK8PFJET500"],
        # "AK8PFJET500":["AK8PFJET550"],
        # "o450":["PFJET450"],
        # "o500":["PFJET500"],
        # "PFJET320":["PFJET450","PFJET500","PFJET550"],
        # "PFJET450":["PFJET500"],
        # "PFJET500":["PFJET550"],
    }
    hists = {
        sample: {year: uproot.open(f"{args.workdir}/{sample}_{year}.root") for year in years} for sample in samples
    }

    for sample in samples:
        axes = setup_axes()
        axes_prescale = setup_axes()
        colors = [
            "#66c2a5",
            "#fc8d62",
            "#8da0cb",
            "#e78ac3",
            "#a6d854",
        ]
        for iyear, year in enumerate(years):
            colors_ = colors[iyear:]
            for ref_trig, probe_triggers in triggers.items():
                for probe_trig in probe_triggers:
                    efficiency_scalefactor(
                        hists, year, ref_trig, probe_trig, outdir=args.outdir, rebin=5, lookup_trigger=True
                    )
            for lookup_trigger in [True, False]:
                lookup_trigger_str = ("_lookup" if lookup_trigger else "")
                if sample == "Data":
                    efficiency_curve(
                        hists[sample][year],
                        year,
                        sample,
                        triggers,
                        use_prescales=True,
                        outdir=args.outdir,
                        rebin=5,
                        lookup_trigger=lookup_trigger,
                    )
                    efficiency_curve(
                        hists[sample][year],
                        year,
                        sample,
                        triggers_year_comparison,
                        use_prescales=True,
                        outdir=args.outdir,
                        rebin=5,
                        axes=axes_prescale,
                        label_suffix=year,
                        colors=colors_,
                        lookup_trigger=lookup_trigger,
                    )
                efficiency_curve(
                    hists[sample][year],
                    year,
                    sample,
                    triggers,
                    use_prescales=False,
                    outdir=args.outdir,
                    rebin=5,
                    lookup_trigger=lookup_trigger,
                )
                efficiency_curve(
                    hists[sample][year],
                    year,
                    sample,
                    triggers_year_comparison,
                    use_prescales=False,
                    outdir=args.outdir,
                    rebin=5,
                    axes=axes,
                    label_suffix=year,
                    colors=colors_,
                    lookup_trigger=lookup_trigger,
                )

            if sample == "Data":
                finalize_axes(
                    axes_prescale,
                    sample,
                    year="RunII",
                    outdir=args.outdir,
                    suffix_str="_prescale" + lookup_trigger_str,
                    data=True,
                )
            finalize_axes(
                axes, sample, year="RunII", outdir=args.outdir, suffix_str=lookup_trigger_str, data=(sample == "Data")
            )
