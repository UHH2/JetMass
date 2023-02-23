#!/usr/bin/env pythonJMS.sh

import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
from matplotlib.ticker import FormatStrFormatter
import re

plt.style.use(hep.style.CMS)
font_size = 20
mpl.rcParams['axes.labelsize'] = font_size

colors_list = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:cyan"]


def plot_JMS(
    JMS,
    names,
    outname="JMSSF.pdf",
    plot_wtagging=False,
    ylim=None,
    massScales_to_plot=["W", "Z", "top"],
    samples_to_plot=[],
):
    linewidth = 2
    markersize = 12
    f, ax = plt.subplots(figsize=(10, 7))
    if ylim is None:
        vals = np.array([par["vals"] for fit in JMS.keys() for par in JMS[fit]["jms"].values()])
        # hi = (vals[:, 0] + vals[:, 1]).max()
        # low = (vals[:, 0] + vals[:, 2]).min()
        hi = (vals[:, 0]).max()
        low = (vals[:, 0]).min()
        magnitude = max(np.log10(abs(hi - 1)), np.log10(abs(low - 1)))
        if magnitude < 0:
            magnitude = np.floor(magnitude)
        else:
            magnitude = np.ceil(magnitude)
        # ylim = (
        #     np.divide(np.ceil(hi/10**(-magnitude) + 1), 10**(-magnitude)),
        #     np.divide(np.floor(low/10**(-magnitude) - 1), 10**(-magnitude))
        # )
        ylim = (
            low*(1-10**magnitude),
            hi*(1+10**magnitude),
        )
        # make symmetric
        ylim_abs = max(abs(ylim[0] - 1), abs(ylim[1] - 1))
        ylim = (1 - ylim_abs, 1 + ylim_abs)
    # hep.label._exp_text("",text="private work",ax=ax, fontsize=font_size)

    def get_year(name):
        year_pattern = re.compile("(?P<year>UL16preVFP|UL16postVFP|UL17|UL18)")
        year = None
        year_match = year_pattern.search(name)
        if year_match:
            year = year_match.groupdict().get("year", "UL17")
        return year

    years = [get_year(name) for name in names]
    if len(np.unique(years)) == 1:
        hep.cms.label(label=", Work in Progress", year=years[0], ax=ax, fontsize=font_size)
    else:
        hep.cms.text(",Work in progress", ax=ax, fontsize=font_size)

    colors = {
        "Combined": "tab:green",
        "VJets": "tab:blue",
        "TTBar": "tab:orange",
    }

    # get pT min and max
    pt_edges = set()
    for fit, results in JMS.items():
        for i in results["jms"].keys():
            pt_edges.add(results["jms"][i]["edges"][0].lower())
            pt_edges.add(results["jms"][i]["edges"][1].lower())
    include_inf = "inf" in pt_edges
    if include_inf:
        pt_edges.remove("inf")
    pt_edges = np.array([float(p) for p in pt_edges])
    pt_edges.sort()
    # pt_min = min(pt_edges)
    pt_max = max(pt_edges)

    if plot_wtagging:
        eb = ax.errorbar(
            [200.0],
            [1.014],
            yerr=0.003,
            xerr=1.0,
            fmt="vk",
            label="W-tagging 2017 ($\\tau_{21}<0.35$)",
            linewidth=linewidth,
            markersize=markersize,
        )
        eb = ax.errorbar(
            [200.0],
            [0.982],
            yerr=0.002,
            xerr=1.0,
            fmt="^k",
            label="W-tagging 2017 ($\\tau_{21}<0.45$)",
            linewidth=linewidth,
            markersize=markersize,
        )

    def fit_match(n):
        fit_pattern = re.compile("(?P<fit>(Combined|VJets|TTBar){1})")
        fit_search = fit_pattern.search(n)
        _fit_match = None
        if fit_search:
            _fit_match = fit_search.groupdict().get("fit", None)
        if _fit_match is None:
            raise RuntimeError("Could not extract fit type from model name {}".format(n))
        return _fit_match

    fit_matches = [fit_match(n) for n in names]
    sample_comp = all(np.unique(fit_matches, return_counts=True)[1] <= 2)
    # iplot = 1
    for iname, name in enumerate(names):
        scales_ = np.array(list(JMS[name]["jms"].keys()))
        edges_ = np.array([JMS[name]["jms"][n]["edges"] for n in JMS[name]["jms"].keys()], dtype=float)
        edges_[edges_ == np.inf] = 1.1 * pt_max

        if any(["W" in scale for scale in scales_]):
            split_prefixe = massScales_to_plot  # ["W","Z","top"]
        else:
            split_prefixe = [""]

        for isplit_prefix, split_prefix in enumerate(split_prefixe):
            marker_style = "x" if isplit_prefix > 0 else "."
            scales = []
            edges = []
            for iscale, scale in enumerate(scales_):
                if split_prefix in scale:
                    if scale == "massScale_pt0_eta0_all_top_PT200":
                        continue
                    scales.append(scale)
                    edges.append(edges_[iscale])
            if len(scales) == 0:
                continue
            edges = np.array(edges)
            # marker_separation = 1 / (len(names) + 2)
            pt_central = 0.5 * (edges[:, 0] + edges[:, 1])
            # pt_central_a = 0.5 * (edges[:, 0] + edges[:, 1])
            if not sample_comp:
                if "UL16preVFP" in name:
                    pt_central = edges[:, 0] + (0.05+0.112+0.225*0) * (edges[:, 1] - edges[:, 0])
                elif "UL16postVFP" in name:
                    pt_central = edges[:, 0] + (0.05+0.112+0.225*1) * (edges[:, 1] - edges[:, 0])
                elif "UL17" in name:
                    pt_central = edges[:, 0] + (0.05+0.112+0.225*2) * (edges[:, 1] - edges[:, 0])
                elif "UL18" in name:
                    pt_central = edges[:, 0] + (0.05+0.112+0.225*3) * (edges[:, 1] - edges[:, 0])

            # pt_central = edges[:,0]+iplot*marker_separation*(edges[:,0]-edges[:,1])
            # iplot +=1

            pt_low = pt_central - edges[:, 0]
            pt_high = edges[:, 1] - pt_central

            scales.sort()
            central = np.array([JMS[name]["jms"][scale]["vals"][0] for scale in scales])
            error_up = np.array([abs(JMS[name]["jms"][scale]["vals"][1]) for scale in scales])
            error_down = np.array([abs(JMS[name]["jms"][scale]["vals"][2]) for scale in scales])

            # jec up/down fits
            jec_unc = False
            jec_up, jec_down = None, None

            def jec_var_name(name, var="UP"):
                eta_region_match = re.search("(BARREL|ENDCAP)", name)
                eta_region = eta_region_match.group() if eta_region_match else ""
                name_ = name.replace(eta_region, "")
                name_ += "" if var == "" else "JEC" + var
                name_ += eta_region
                return name_

            if (jec_var_name(name, "UP") in JMS.keys() and jec_var_name(name, "DOWN") in JMS.keys()) and (
                jec_var_name(name, "UP") not in names and jec_var_name(name, "DOWN") not in names
            ):
                jec_up = central - np.array(
                    [JMS[jec_var_name(name, "UP")]["jms"][scale]["vals"][0] for scale in scales]
                )
                jec_up = np.where(jec_up < 0, 0, jec_up)
                jec_down = (
                    np.array([JMS[jec_var_name(name, "DOWN")]["jms"][scale]["vals"][0] for scale in scales]) - central
                )
                jec_down = np.where(jec_down < 0, 0, jec_down)
                jec_symmetric = np.maximum(jec_up, jec_down)
                jec_unc = True

            legend_label = name + " " + split_prefix

            # legend_label = legend_label.replace("NOJEC", " (no JEC on msd)")
            if "JEC" in name and "UP" not in name and "DOWN" not in name:
                legend_label = legend_label.replace("NOJEC", "")
            else:
                legend_label += r" (JEC on $m_\mathrm{SD}$)"
                legend_label = legend_label.replace("JECUP", " (jec-up)")
                legend_label = legend_label.replace("JECDOWN", " (jec-down)")
            legend_label = legend_label.replace("Combined", "combined ")
            legend_label = legend_label.replace("VJets", "W(qq)+jets ")
            legend_label = legend_label.replace("TTBar", r"$t\bar{t}$ ")
            legend_label = legend_label.replace("BARREL", "")
            legend_label = legend_label.replace("ENDCAP", "")
            color = colors[fit_matches[iname]] if sample_comp else colors_list[iname]
            if "JEC" not in name:
                marker_style = "v"
            eb_jec = None
            if jec_unc:
                # eb_sum =
                eb_jec = ax.errorbar(
                    pt_central,
                    central,
                    yerr=[error_down+jec_symmetric, error_up+jec_symmetric],
                    xerr=[pt_low, pt_high],
                    fmt=marker_style,
                    color=color,
                    linewidth=linewidth,
                    markersize=markersize,
                    alpha=0.7,
                )
                 # eb_sum[2][1].set_linestyle("-.")

            eb = ax.errorbar(
                pt_central,
                central,
                yerr=[error_down, error_up],
                xerr=None if jec_unc else [pt_low, pt_high],
                fmt=marker_style,
                color=color,
                label=legend_label,
                linewidth=linewidth,
                markersize=markersize,
                alpha=0.7,
                capsize=10 if jec_unc else 0.0,
            )

            # line_type = "-" if "VJets" in name else "-"
            line_type = "--" if "JEC" not in name else "-"            
            if eb_jec:
                # eb[1][0].set_linestyle(line_type)
                eb_jec[-1][0].set_linestyle(line_type)
                eb_jec[-1][1].set_linestyle(line_type)
            else:
                eb[-1][0].set_linestyle(line_type)

        # iplot = +1
    ax.plot(ax.get_xbound(), [1.0, 1.0], "k--")
    if ylim:
        ax.set_ylim(*ylim)
    # ax.set_ylim(0.97,1.03)
    ax.set_ylabel("JMS-SF", loc="center")
    ax.set_xlabel("$p_T$ [GeV]")

    if all(["ENDCAP" in name for name in names]):
        ax.text(ax.get_xlim()[1]*0.7, 1.04, r"$1.3 < |\eta| < 2.4$", fontsize=font_size - 2)
    if all(["BARREL" in name for name in names]):
        ax.text(ax.get_xlim()[1]*0.7, 1.04, r"$|\eta| < 1.3$", fontsize=font_size - 2)

    # pT_labels = np.array(['inclusive']+pT_ if plot_wtagging else pT_)
    # plt.xticks(
    #     rotation=45, ticks=(range(-1, len(pT_labels) - 1) if plot_wtagging else range(len(pT_labels))),
    #     labels=pT_labels
    # )
    plt.xticks(fontsize=font_size - 2)
    plt.yticks(fontsize=font_size - 2)
    # ax.legend(fontsize=font_size - 2, loc="upper left", bbox_to_anchor=(1, 1))
    ax.legend(fontsize=font_size - 6, loc="upper left")
    # ax.legend(fontsize = font_size-5,loc='upper left', ncol=1,frameon=True)
    for ext in ["pdf", "png"]:
        f.savefig(outname.replace(outname.split(".")[-1], ext), bbox_inches="tight")


def plot_mass(
    JMS,
    names,
    outname="mass",
    plotJMS=False,
    plot_wtagging=False,
    ylim=None,
    regions_to_plot=["pass", "passW"],
    prefit_postfit="postfit",
    samples_to_plot=["ZJetsMatched"],
    skipfilter=None,
    skip_colors=0,
    yerr=True,
):
    linewidth = 2
    colors_list = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    ]
    markersize = 12
    mpl.rcParams["xtick.major.pad"] = "8"

    f, ax = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)

    # hep.label._exp_text("",text="private work",ax=ax, fontsize=font_size)
    hep.cms.text(",work in progress", ax=ax[0], fontsize=font_size)

    # get pT min and max
    pt_edges = set()
    for fit, results in JMS.items():
        for channel in results["peak_positions"].values():
            for region in channel.values():
                for sample in region.values():
                    pt_edges.add(sample["prefit"]["pt_edges"][0])
                    pt_edges.add(sample["prefit"]["pt_edges"][1])
    include_inf = "inf" in pt_edges
    if include_inf:
        pt_edges.remove("inf")
    pt_edges = np.array([float(p) for p in pt_edges])
    pt_edges.sort()
    # pt_min = min(pt_edges)
    # pt_max = max(pt_edges)

    iplot = 0 + skip_colors
    for iname, name in enumerate(names):

        samples = set()
        regions = set()
        channels = set()
        for ch_name, ch in JMS[name]["peak_positions"].items():
            channels.add(ch_name)
            for region_name, region in ch.items():
                regions.add(region_name)
                for sample_name, sample in region.items():
                    # print(sample_name,sample)
                    samples.add(sample_name)
        channels = list(channels)
        channels.sort()
        for region in regions:
            if region not in regions_to_plot:
                continue
            for sample in samples:
                if sample not in samples_to_plot:
                    continue
                if skipfilter is not None:
                    if skipfilter(region + "_" + sample):
                        continue

                results = JMS[name]["peak_positions"]
                # print(results)
                prefit_mu = np.array([results[ch][region][sample]["prefit"]["values"][0] for ch in channels])
                prefit_sigma = np.array([abs(results[ch][region][sample]["prefit"]["values"][1]) for ch in channels])

                postfit_mu = np.array([results[ch][region][sample]["postfit"]["values"][0] for ch in channels])
                postfit_sigma = np.array([abs(results[ch][region][sample]["postfit"]["values"][1]) for ch in channels])

                edges = np.array([results[ch][region][sample]["prefit"]["pt_edges"] for ch in channels])

                pt_central = 0.5 * (edges[:, 0] + edges[:, 1])
                # print(pt_central)
                pt_low = pt_central - edges[:, 0]
                pt_high = edges[:, 1] - pt_central

                legend_label = name + " " + region + " " + sample

                if prefit_postfit == "prefit" or prefit_postfit == "both":
                    eb_prefit = ax[0].errorbar(
                        pt_central,
                        prefit_mu,
                        yerr=prefit_sigma if yerr else 0.0,
                        xerr=[pt_low, pt_high],
                        fmt=".",
                        color=colors_list[iplot],
                        label=legend_label + " prefit",
                        linewidth=linewidth,
                        markersize=markersize,
                        alpha=0.7,
                    )
                    eb_prefit[-1][0].set_linestyle("--")
                    ax[1].plot(
                        pt_central,
                        prefit_sigma / prefit_mu,
                        color=colors_list[iplot],
                        linewidth=linewidth,
                        linestyle="--",
                        markersize=markersize,
                        alpha=0.7,
                        marker=".",
                    )
                if prefit_postfit == "postfit" or prefit_postfit == "both":
                    eb_prefit = ax[0].errorbar(
                        pt_central,
                        postfit_mu,
                        yerr=postfit_sigma if yerr else 0.0,
                        xerr=[pt_low, pt_high],
                        fmt="x",
                        color=colors_list[iplot],
                        label=legend_label + " postfit",
                        linewidth=linewidth,
                        markersize=markersize,
                        alpha=0.7,
                    )
                    eb_prefit[-1][0].set_linestyle("-")
                    ax[1].plot(
                        pt_central,
                        postfit_sigma / postfit_mu,
                        color=colors_list[iplot],
                        linewidth=linewidth,
                        linestyle="-",
                        markersize=markersize,
                        alpha=0.7,
                        marker="x",
                    )

                iplot += 1

    if ylim:
        ax[0].set_ylim(*ylim)

    msd_bounds = ax[0].get_ybound()
    pt_bounds = ax[0].get_xbound()
    if msd_bounds[1] > 80:
        ax[0].plot(pt_bounds, [80.379, 80.379], "k--")  # W-Mass
    if msd_bounds[1] > 90:
        ax[0].plot(pt_bounds, [91.1876, 91.1876], "k--")  # Z-Mass
    if msd_bounds[1] > 170:
        ax[0].plot(pt_bounds, [172.76, 172.76], "k--")  # Top-Mass

    ratio_bounds = ax[1].get_ybound()
    ax[1].set_ylim(ratio_bounds[0] * 0.9, ratio_bounds[1] * 1.1)

    ax[0].set_ylabel(r"$\mu$ [GeV]", loc="center")
    ax[1].set_ylabel(r"$\sigma /\ \mu$", loc="center")
    ax[1].set_xlabel(r"$p_T$ [GeV]")

    ax[0].yaxis.set_major_formatter(FormatStrFormatter("%7.0f"))

    plt.xticks(fontsize=font_size - 2)
    plt.yticks(fontsize=font_size - 2)
    ax[0].tick_params(axis="y", labelsize=font_size - 2)
    ax[0].legend(fontsize=font_size - 2, loc="upper left", bbox_to_anchor=(1, 1))
    for ext in ["pdf", "png"]:
        f.savefig(outname.replace(outname.split(".")[-1], ext), bbox_inches="tight")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fitResults",
        "-f",
        default="fitResults.json",
        help="File containing FitResults in JSON-Format (generated by extractMasscales.py)",
    )
    parser.add_argument(
        "--fits",
        nargs="+",
        default=[],
        help="List of fits (separated by spaces) to be included in the plot.",
        required=True,
    )
    parser.add_argument(
        "--outName",
        "-o",
        default="JMSSF.pdf",
        help="Name of output-file(s). Regardless of provided file-type, plot will be saved both as pdf and png.",
    )
    parser.add_argument(
        "--plotCMSWTag", action="store_true", help="Add a marker with the JMS-SF from CMS W-tagging studies."
    )
    parser.add_argument("--plot", "-p", default="JMS", choices=["JMS", "mass"], help="Method to use for plotting.")
    parser.add_argument(
        "--stage",
        "-s",
        default="postfit",
        help="For plotting the mass: take mass from >prefit< or >postfit< templates.",
        choices=["prefit", "postfit", "both"],
    )
    parser.add_argument("--samples", default=["WJetsMatched", "ZJetsMatched"])

    args = parser.parse_args()

    fit_results = json.load(open(args.fitResults))
    if any(fit not in fit_results.keys() for fit in args.fits):
        print("provided fit names:", args.fits)
        raise AttributeError(
            "At least one of the provided names of fits is not included in the provided fitResults JSON! "
            f"Available fits: {fit_results.keys()}"
        )
    for fit in args.fits:
        print(fit)

    kwargs = {}
    if args.plot == "mass":
        kwargs["prefit_postfit"] = args.stage

    plotting_functions = {f_: globals()[f_] for f_ in globals() if "plot_" in f_}
    plotting_function_name = f"plot_{args.plot}"
    if plotting_function_name not in plotting_functions:
        raise NotImplementedError(
            f"The plotting function you asked for is not implemented. Available functions: {plotting_functions.keys()}"
        )

    ylim = (0.9478006051858293, 1.0521993948141706)

    plotting_functions[plotting_function_name](
        fit_results,
        args.fits,
        args.outName,
        plot_wtagging=args.plotCMSWTag,
        samples_to_plot=args.samples,
        **kwargs, ylim=ylim
    )
