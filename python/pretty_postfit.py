#!/usr/bin/env pythonJMS.sh
import uproot
import os
import hist
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import numpy as np
import itertools
import logging
import unfolding_plotting
from coffea.util import load

hep.style.use("CMS")

diverging_colors = [
    "#a50026",
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",

    "#e0f3f8",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#313695",

    "#40004b",
    "#762a83",
    "#9970ab",
    "#c2a5cf",
    "#e7d4e8",

    "#d9f0d3",
    "#a6dba0",
    "#5aae61",
    "#1b7837",
    "#00441b",
]

sequential_colors = list(reversed([
    '#fee6ce',
    '#fdd0a2',
    '#fdae6b',
    '#fd8d3c',
    '#f16913',
    '#d94801',
    '#a63603',
    '#7f2704'
]))+list(reversed([
    '#efedf5',
    '#dadaeb',
    '#bcbddc',
    '#9e9ac8',
    '#807dba',
    '#6a51a3',
    '#54278f',
    '#3f007d',
]))+list(reversed([

    '#e5f5e0',
    '#c7e9c0',
    '#a1d99b',
    '#74c476',
    '#41ab5d',
    '#238b45',
    '#006d2c',
    '#00441b',
]))

lumis = {
    "UL16preVFP": 19301.591954787407 / 1000.0,
    "UL16postVFP": 16626.734093195286 / 1000.0,
    "UL17": 41479.68052876168 / 1000.0,
    "UL18": 59832.47533908866 / 1000.0,
}
lumis["RunII"] = sum(lumis[year] for year in ["UL16preVFP", "UL16postVFP", "UL17", "UL18"])


mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=diverging_colors)


def cms_label(exp_label, year, fs, ax, data):
    if "private" in exp_label.lower():
        hep.label.exp_label(
            exp="",
            llabel=exp_label,
            year=year,
            ax=ax,
            fontsize=fs,
            data=data
        )
    else:
        hep.cms.label(exp_label, year=year, fontsize=fs, ax=ax, data=data)


def plot_templates(
        fit_dir: str,
        fs: int = 20, max_plots: int = 1000,
        logy: bool = False,
        year: str = "UL17",
        exp_label: str = "Private work (CMS simulation)",
        data: bool = False,
        region: str = "pass",
):
    logger = logging.getLogger()
    fit_shapes = uproot.open(os.path.join(args.fit_dir, "fit_shapes.root"))

    configs = json.loads(open(os.path.join(args.fit_dir, "config.json"), "r").read())
    ptgen_edges = configs["unfolding_bins"]["ptgen"]
    msdgen_edges = configs["unfolding_bins"]["msdgen"]

    # getting templates
    pt_reco_edges = configs["pt_edges"]
    if "Inf" in pt_reco_edges:
        pt_reco_edges[pt_reco_edges.index("Inf")] = np.inf
    genbin_colors = {
        f"ptgen{ipt}_msdgen{imsd}": diverging_colors[imsd+ipt*len(msdgen_edges)]
        for ipt in range(len(ptgen_edges) - 1)
        for imsd in range(len(msdgen_edges) - 1)
    }
    for ptreco_low, ptreco_hi in hist.axis.Variable(pt_reco_edges):
        templates = {
            state: {
                f"ptgen{ipt}_msdgen{imsd}": {
                    "pt_low": ptgen_edges[ipt],
                    "pt_hi": ptgen_edges[ipt + 1],
                    "msd_low": msdgen_edges[imsd],
                    "msd_hi": msdgen_edges[imsd + 1],
                }
                for ipt in range(len(ptgen_edges) - 1)
                for imsd in range(len(msdgen_edges) - 1)
            }
            for state in ["prefit", "postfit"]
        }

        # getting actual histograms
        for state, bins in templates.items():
            for genbin in bins.keys():
                hist_name = f"WJetsPt{ptreco_low:.0f}{region}_{state}/WJetsMatched_{genbin}"
                if hist_name not in fit_shapes:
                    logger.info(f"{state} template in genbin {genbin} does not exist.")
                    continue
                templates[state][genbin]["hist"] = fit_shapes[hist_name]

        for state in ["prefit", "postfit"]:
            if max_plots < 0:
                break
            max_plots -= 1

            f, ax = plt.subplots(figsize=(9, 9))
            # hists = []
            # colors = []
            labels = ["fakes"]
            hists = [fit_shapes[f"WJetsPt{ptreco_low:.0f}{region}_{state}/WJetsMatched_fakes"]]
            colors = ["tab:grey"]

            for genbin_name, genbin in templates[state].items():
                if "hist" in genbin:
                    hists.append(genbin["hist"])
                    colors.append(genbin_colors[genbin_name])
                    pt_gen_bin_tex = r"$p_{T,gen} \in" + f"[{genbin['pt_low']},{genbin['pt_hi']})$"
                    msd_gen_bin_tex = r"$m_{SD,gen} \in" + f"[{genbin['msd_low']},{genbin['msd_hi']})$"
                    labels.append(pt_gen_bin_tex + " " + msd_gen_bin_tex)
            hep.histplot(hists, label=labels, stack=True, histtype="fill", color=colors)

            cms_label(exp_label=exp_label, year=year, fs=fs, ax=ax, data=data and state != "prefit")

            ax.text(
                np.diff(ax.get_xlim()) * 0.4,
                ax.get_ylim()[1] * 0.2,
                r"$\it{" + state + r"}~W(qq)\mathrm{+jets}$",
                fontsize=fs + 2,
            )
            ax.text(
                np.diff(ax.get_xlim()) * 0.4,
                ax.get_ylim()[1] * 0.15,
                r"$%.0f~\mathrm{GeV} \leq p_{T,\mathrm{reco}} < %.0f~\mathrm{GeV}$" % (ptreco_low, ptreco_hi),
                fontsize=fs + 2,
            )

            ax.legend(fontsize=fs - 6)
            ax.set_ylabel("Events")
            ax.set_xlabel(r"$m_{SD,\mathrm{reco}}$ [GeV]")
            # if logy:
            #     ax.set_yscale("log")
            #     ax.set_ylim(10**(-3), 10**(3))
            # f.tight_layout()

            plot_dir = os.path.join(fit_dir, "plots/pretty_signal")
            if not os.path.isdir(plot_dir):
                os.makedirs(plot_dir)
            for do_logY in [False]:
                if do_logY:
                    ax.set_yscale("log")
                    ax.set_ylim(10**(-3), 10**(3))

                plt.savefig(
                    f"{plot_dir}/WJetsMatchedSignal_{ptreco_low:.0f}_{region}_{state}{'_logY' if do_logY else ''}.pdf",
                    bbox_inches="tight",
                    # transparent=True,
                    pad_inches=0.01,
                )
        # os.system(
        #     "convert -delay 100 "
        #     + f"{plot_dir}/WJetsMatchedSignal_{ptreco_low:.0f}_prefit.pdf "
        #     + f"{plot_dir}/WJetsMatchedSignal_{ptreco_low:.0f}_postfit.pdf "
        #     + f"{plot_dir}/WJetsMatchedSignal_{ptreco_low:.0f}.gif"
        # )


def finite_edges(axis):
    edges = axis.edges
    if edges[-1] == np.inf:
        edges[-1] = edges[-2] + axis.widths[-2] * 2
    labels = [f"{e:.0f}".lower().replace("inf", r"$\infty$") for e in axis.edges]
    return edges, labels


def unrolled_edge_positions(edges, n_bins_super):
    bin_widths = (edges[1:]-edges[:-1])
    edge_positions = np.concatenate(([0], bin_widths / bin_widths[0]))
    unrolled_positions = np.cumsum(
        np.concatenate((np.array([0]), np.array([edge_positions[1:]] * n_bins_super).flatten()))
    )
    return unrolled_positions


def plot_migration_matrix1(migmat: hist.Hist, outname: str, r_gen=1.0):
    f = plt.figure(figsize=(11, 12))
    grid = f.add_gridspec(1, 2, width_ratios=[1, 0.02])
    ax = f.add_subplot(grid[0])
    cax = f.add_subplot(grid[1])

    pt_reco_edges, pt_reco_labels = finite_edges(migmat.axes["pt_reco"])
    nbins_pt_reco = len(pt_reco_edges) - 1
    pt_gen_edges, pt_gen_labels = finite_edges(migmat.axes["pt_gen"])
    nbins_pt_gen = len(pt_gen_edges) - 1

    msd_reco_edges, msd_reco_labels = finite_edges(migmat.axes["msd_reco"])
    nbins_msd_reco = len(msd_reco_edges) - 1
    msd_gen_edges, msd_gen_labels = finite_edges(migmat.axes["msd_gen"])
    nbins_msd_gen = len(msd_gen_edges) - 1

    if isinstance(r_gen, np.ndarray):
        if len(r_gen) != nbins_msd_gen*nbins_pt_gen:
            raise RuntimeError("r_gen has to have size of nbind_msd_gen*nbins_pt_gen!!")
        r_gen = np.multiply(
            np.ones((
                nbins_pt_gen * nbins_msd_gen,
                nbins_pt_reco * nbins_msd_reco,
            )),
            r_gen[:, np.newaxis]
        )

    gen_positions = unrolled_edge_positions(msd_gen_edges, nbins_pt_gen)
    reco_positions = unrolled_edge_positions(msd_reco_edges, nbins_pt_reco)

    pt_reco_labels_positions = [reco_positions[ipt * nbins_msd_reco] for ipt in range(nbins_pt_reco + 1)]
    pt_gen_labels_positions = [gen_positions[ipt * nbins_msd_gen] for ipt in range(nbins_pt_gen + 1)]
    content = (
        migmat.values().reshape(
            nbins_pt_gen * nbins_msd_gen,
            nbins_pt_reco * nbins_msd_reco,
        )
        * r_gen
    )

    c = ax.pcolormesh(
        *np.meshgrid(
            gen_positions,
            reco_positions,
        ),
        content.T,
        cmap="magma",
        norm=LogNorm(),
    )

    ax.set_xticks(pt_gen_labels_positions, pt_gen_labels)
    for iptgen in range(1, nbins_pt_gen):
        ax.plot([pt_gen_labels_positions[iptgen]] * 2, ax.get_ylim(), "k--", alpha=0.6)
    ax.set_xlabel("generator bin")

    ax.set_yticks(pt_reco_labels_positions, pt_reco_labels)
    for iptreco in range(1, nbins_pt_reco):
        ax.plot(ax.get_xlim(), [pt_reco_labels_positions[iptreco]] * 2, "k--", alpha=0.6)
    ax.set_ylabel("detector bin")

    msd_gen_subax_length = gen_positions[nbins_msd_gen]
    msd_reco_subax_length = reco_positions[nbins_msd_reco]

    # pt_gen_underflow= patches.Rectangle((0,0),msd_gen_subax_length,msd_reco_subax_length*nbins_pt_reco,
    #                                   alpha=0.1,facecolor="green",label="under-&overflow bins")
    # ax.add_patch(pt_gen_underflow)

    # pt_gen_underflow= patches.Rectangle((0,0),msd_gen_subax_length,msd_reco_subax_length*nbins_pt_reco,
    #                                   alpha=0.1,facecolor="green",label="under-&overflow bins")
    # ax.add_patch(pt_gen_underflow)

    for iptgen, iptreco in itertools.product(range(nbins_pt_gen), range(nbins_pt_reco)):
        min_gen = msd_gen_edges.min()
        max_gen = msd_gen_edges.max()
        min_reco = msd_reco_edges.min()
        max_reco = msd_reco_edges.max()
        min_ = max(min_gen, min_reco)
        max_ = min(max_gen, max_reco)
        gen_pos = (
            (min_ - min_gen) * msd_gen_subax_length / (max_gen - min_gen) + iptgen * msd_gen_subax_length,
            (max_ - min_gen) * msd_gen_subax_length / (max_gen - min_gen) + iptgen * msd_gen_subax_length,
        )
        reco_pos = (
            (min_ - min_reco) * msd_reco_subax_length / (max_reco - min_reco) + iptreco * msd_reco_subax_length,
            (max_ - min_reco) * msd_reco_subax_length / (max_reco - min_reco) + iptreco * msd_reco_subax_length,
        )
        ax.plot(gen_pos, reco_pos, "k--", alpha=0.6)

    f.colorbar(c, cax=cax, orientation="vertical", label="Events")
    f.savefig(
        outname,
        bbox_inches="tight"
    )


def plot_unfolded_mass(
    fit_dir,
    data: bool = False,
    exp_label: str = "Private work (CMS simulation)",
    plot_truth: bool = False,
    year: str = "UL17",
    yaxis: str = "au",
    inclusive_tagger: bool = True,
    plot_no_matching: bool = False,
    plot_matching_comp: bool = False,
    n2cut: str = "",
):
    import json
    from copy import deepcopy
    configs = json.load(open(f"{fit_dir}/config.json", "r"))
    unfolding_fit_results = json.load(open(f'{fit_dir}/{configs["ModelName"]}fitResult.json', "r"))

    # f_ = uproot.open(configs["histLocation"])
    hist_location = configs["histLocation"]
    # import re
    # year_match = re.search("(?P<year>UL1[678]{1}(preVFP|postVFP)*)", hist_location)
    # year_hist_location = year
    # if year_match:
    #     year_hist_location = year_match.groupdict()["year"]
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"] if year == "RunII" else [year]
    tagger = "_particlenetDDT" if "particlenetddt" in hist_location.lower() else ""
    print("loading files")
    files = [
        load(f"/nfs/dust/cms/user/albrechs/JetMassFits/coffea_hists/templates_{year}{tagger}_mctruth.coffea")
        for year in years
    ]
    n2cut_str = "no_n2" if n2cut == "" else n2cut
    # n2_0p124

    theory_systs = ["v_qcd", "w_ewk"]
    theory_vars = [f"{syst}_{direc}" for syst in theory_systs for direc in ["up", "down"]]

    data_label = ("unfolded data" if data else "unfolded pseudo data") + r"(stat. $\bigoplus$ syst. unc.)"
    out_dir = f"{fit_dir}/plots/pretty_unfold_{yaxis}{n2cut}/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    def get_edges(axis):
        edges = deepcopy(configs["unfolding_bins"][axis])
        if edges[-1] == "Inf":
            edges[-1] = 2 * edges[-2] + edges[-3]
        return edges
    msd_edges = get_edges("msdgen")
    pt_edges = get_edges("ptgen")

    region = "inclusive" if inclusive_tagger else "pass"
    region_str = "inclusive" if inclusive_tagger else ""

    matchings = ["matching"]
    if plot_no_matching:
        matchings = ["no matching"]
    if plot_matching_comp:
        matchings = ["matching", "no matching"]

    matching_str = {
        "matching": " (with $W$ matching)",
        "no matching": " (without $W$ matching)",
    }
    if not plot_no_matching and not plot_matching_comp:
        matching_str = {
            "no matching": "",
            "matching": "",
        }

    if n2cut == "":
        matching_kwargs = {"matching": dict(color="tab:blue"), "no matching": dict(color="tab:red")}
    else:
        matching_kwargs = {"matching": dict(color="tab:red"), "no matching": dict(color="tab:red")}

    acceptance_file_matching_str = {
        "matching": "",
        "no matching": "_nomatching",
    }
    acceptance_files = {
        matching: [
            load(
                f"acceptance_efficiency_plots_{n2cut_str}{acceptance_file_matching_str[matching]}/" +
                f"misses_acceptance{year}.coffea"
            )
            for year in years
        ]
        for matching in matchings
    }

    acceptance_files_no_n2 = {
        matching: [
            load(
                f"acceptance_efficiency_plots_no_n2{acceptance_file_matching_str[matching]}/" +
                f"misses_acceptance{year}.coffea"
            )
            for year in years
        ]
        for matching in matchings
    }

    print("getting hists")

    acceptance_years = {
        matching: [
            [
                acceptance_files[matching][ifile]["acceptance"][iptgen][0]
                * (
                    acceptance_files_no_n2[matching][ifile]["total"][{"ptgen": iptgen}].values()
                    / acceptance_files[matching][ifile]["total"][{"ptgen": iptgen}].values()
                )
                for ifile in range(len(acceptance_files[matching]))
            ]
            for iptgen in range(len(pt_edges) - 1)
        ]
        for matching in matchings
    }

    # acceptance_years_no_n2 = [
    #     [
    #         acceptance_file["acceptance"][iptgen][0]
    #         for acceptance_file in acceptance_files_no_n2
    #     ]
    #     for iptgen in range(len(pt_edges) - 1)
    # ]

    efficiency_years = {
        matching: [
            [
                # acceptance_file["sf_efficiency"]["pass_gen_pass_reco_pass_dR"][iptgen][0]
                acceptance_file["pass_gen_pass_reco_pass_dR_withSFHEM"][{"ptgen": iptgen}].values()
                / acceptance_file["pass_gen_pass_reco_pass_dR"][{"ptgen": iptgen}].values()
                for acceptance_file in acceptance_files[matching]
            ]
            for iptgen in range(len(pt_edges) - 1)
        ]
        for matching in matchings
    }

    truth_mc_noreco_years = {
        matching: [
            [
                acceptance_file["total"][{"ptgen": iptgen}]
                for acceptance_file in acceptance_files[matching]
            ]
            for iptgen in range(len(pt_edges) - 1)
        ]
        for matching in matchings
    }

    truth_mc_years = {
        "matching": [
            [
                deepcopy(f_[f"vjets_mjet_unfolding_{region}"][{"ptgen": iptgen, "dataset": "vjets_WJetsMatched"}])
                for f_ in files
            ]
            for iptgen in range(len(pt_edges) - 1)
        ]
    }
    if "no matching" in matchings:
        truth_mc_years["no matching"] = [
            [
                deepcopy(f_["pass_gen_pass_reco_pass_dR_withSFHEM"][{"ptgen": iptgen}])
                for f_ in acceptance_files["no matching"]
            ]
            for iptgen in range(len(pt_edges) - 1)
        ]

    del files

    var_files = {
        var: [
            load(f"/nfs/dust/cms/user/albrechs/JetMassFits/coffea_hists/templates_{year}{tagger}_{var}_mctruth.coffea")
            for year in years
        ]
        for var in theory_vars
    }
    theory_vars_hists_years = {
        matching: {
            var: [
                [
                    deepcopy(
                        f_var[f"vjets_mjet_unfolding_{region}"][
                            {
                                "ptgen": iptgen,
                                "dataset": "vjets_WJetsMatched" if matching == "matching" else "vjets_WJets",
                            }
                        ]
                    )
                    for f_var in var_files[var]
                ]
                for iptgen in range(len(pt_edges) - 1)
            ]
            for var in theory_vars
        }
        for matching in matchings
    }

    del var_files

    print("adding hists")

    truth_mc = {
        matching: [deepcopy(hists[0]) for hists in truth_mc_years[matching]]
        for matching in matchings
    }

    truth_mc_noreco = {
        matching: [
            deepcopy(hists[0])
            for hists in truth_mc_noreco_years[matching]
        ]
        for matching in matchings
    }

    theory_vars_hists = {
        matching: {
            var: [deepcopy(hists[0]) for hists in theory_vars_hists_years[matching][var]]
            for var in theory_vars
        }
        for matching in matchings
    }
    for matching in matchings:
        for iptgen in range(len(pt_edges) - 1):
            truth_mc[matching][iptgen] /= acceptance_years[matching][iptgen][0] * efficiency_years[matching][iptgen][0]
            for var in theory_vars:
                theory_vars_hists[matching][var][iptgen] /= (
                    acceptance_years[matching][iptgen][0] * efficiency_years[matching][iptgen][0]
                )
    for iyear in range(1, len(years)):
        for iptgen in range(len(pt_edges) - 1):
            for matching in matchings:
                truth_mc_noreco[matching][iptgen] += deepcopy(truth_mc_noreco_years[matching][iptgen][iyear])
                truth_mc[matching][iptgen] += deepcopy(truth_mc_years[matching][iptgen][iyear]) / (
                    acceptance_years[matching][iptgen][iyear] * efficiency_years[matching][iptgen][iyear]
                )

                for var in theory_vars:
                    theory_vars_hists[matching][var][iptgen] += deepcopy(
                        theory_vars_hists_years[matching][var][iptgen][iyear]
                    ) / (acceptance_years[matching][iptgen][iyear] * efficiency_years[matching][iptgen][iyear])

    mc_truth_noreco_sum = None
    mc_truth_noreco_variance_sum = None
    mc_truth_sum = None
    mc_truth_variance_sum = None
    # capsize = 5
    fontsize = 14
    legend_fontsize = fontsize + 4
    unfolding_sum = None
    unfolding_variance_sum = None
    f_all, ax_all = plt.subplots(figsize=(10, 10))
    theory_upper_sum = None
    theory_lower_sum = None
    ymaxs = [100.0, 23.0, 7.0, 0.6]
    ymax_sum = 122.0
    if n2cut != "":
        ymaxs = [50., 10., 3., 0.3]
        ymax_sum = 50.

    # markers = ["o", "v", "^"]
    # linestyle = ["-", "--", ":"]
    marker_kwargs = {
        "no matching": dict(fillstyle="none", markeredgewidth=2),
        "matching": dict(),
    }

    for ipt in range(0, len(pt_edges) - 1):
        # plot_kwargs = {
        #     "mc": dict(ls=linestyle[ipt]),
        #     "data": dict(fmt=markers[ipt], markersize=6)
        # }
        f, ax = plt.subplots(figsize=(10, 10))
        pt_bin_tex = r" ($%s \leq p_{T,\mathrm{truth}} < %s $)" % (
            str(configs["unfolding_bins"]["ptgen"][ipt]),
            str(configs["unfolding_bins"]["ptgen"][ipt + 1]),
        )

        flat_scale = {match: 1.0 for match in matchings}
        y_label = "Events"
        binwnorm = None if yaxis == "au" else 1.0

        if yaxis == "au":
            for matching in matchings:
                flat_scale[matching] = 1.0 / truth_mc[matching][ipt].values().sum()
            y_label = "a.u."
        elif yaxis == "sigma":
            for matching in matchings:
                flat_scale[matching] = 1.0 / lumis[year]
            y_label = r"$\frac{d\sigma}{d m_\mathrm{SD}}~[\frac{fb}{\mathrm{GeV}}]$"

        x_label = r"$m_{\mathrm{SD, gen}} [GeV]$"
        msd_max = 260.
        msd_min = 30.
        msd_edges_ = truth_mc["matching"][ipt].axes[0].edges.copy()
        msd_edges_[0] = msd_min
        msd_edges_[-1] = msd_max
        msd_centers = (msd_edges_[:-1]+msd_edges_[1:])/2
        msd_xerr = [msd_centers - msd_edges_[:-1], msd_edges_[1:] - msd_centers]
        msd_edges_for_binwidth = deepcopy(msd_edges)
        # setting last bin edge to 1050
        # this number is the upper bound of mtruth values in the signal region
        # (estimated from tinyTrees (~1042))
        msd_edges_for_binwidth[-1] = 1050.
        scale = {matching: 1.0 for matching in matchings}
        if binwnorm:
            scale = {matching: binwnorm / np.diff(msd_edges_for_binwidth) for matching in matchings}

        truth_noreco_values = {
            matching: truth_mc_noreco[matching][ipt].values() * scale[matching] * flat_scale[matching]
            for matching in matchings
        }
        truth_noreco_variances = {
            matching: (
                truth_mc_noreco[matching][ipt].variances()
                * scale[matching]
                * scale[matching]
                * flat_scale[matching]
                * flat_scale[matching]
            )
            for matching in matchings
        }

        truth_values = {
            matching: truth_mc[matching][ipt].values() * flat_scale[matching] * scale[matching]
            for matching in matchings
        }
        truth_variances = {
            matching: truth_mc[matching][ipt].variances()
            * flat_scale[matching]
            * flat_scale[matching]
            * scale[matching]
            * scale[matching]
            for matching in matchings
        }

        theory_var_values = {
            matching: {
                var: theory_vars_hists[matching][var][ipt].values() * flat_scale[matching] * scale[matching]
                for var in theory_vars
            }
            for matching in matchings
        }

        theory_lower = {
            matching: {
                syst: np.abs(
                    np.min(
                        np.array([theory_var_values[matching][f"{syst}_{direc}"] for direc in ["up", "down"]]), axis=0
                    )
                    - truth_values[matching]
                )
                for syst in theory_systs
            }
            for matching in matchings
        }
        theory_upper = {
            matching: {
                syst: np.abs(
                    np.max(
                        np.array([theory_var_values[matching][f"{syst}_{direc}"] for direc in ["up", "down"]]), axis=0
                    )
                    - truth_values[matching]
                )
                for syst in theory_systs
            }
            for matching in matchings
        }

        theory_band_low = {
            matching: theory_lower[matching]["v_qcd"] ** 2 + theory_lower[matching]["w_ewk"] ** 2
            for matching in matchings
        }
        theory_band_hi = {
            matching: theory_upper[matching]["v_qcd"] ** 2 + theory_upper[matching]["w_ewk"] ** 2
            for matching in matchings
        }

        unfolding_signal_strengths = np.array(
            [
                unfolding_fit_results.get(f"r_ptgen{ipt}_msdgen{imsdgen}", [1, 0, 0])[0]
                for imsdgen in range(len(msd_edges) - 1)
            ]
        )
        unfolding_uncertainties = np.array(
            [
                [
                    abs(unfolding_fit_results.get(f"r_ptgen{ipt}_msdgen{imsdgen}", [1, 0, 0])[2]),
                    abs(unfolding_fit_results.get(f"r_ptgen{ipt}_msdgen{imsdgen}", [1, 0, 0])[1])
                ]
                for imsdgen in range(len(msd_edges) - 1)
            ]
        ).T

        unfolding_values = {m: unfolding_signal_strengths*deepcopy(truth_values[m]) for m in matchings}
        unfolding_variances = {m: (deepcopy(truth_values[m]) * unfolding_uncertainties)**2 for m in matchings}
        # unfolding_variances = (unfolding_signal_strengths**2 * truth_variances) + (
        #     truth_values * unfolding_uncertainties
        # ) ** 2
        if True:  # ipt == 1 or ipt == 2:
            if mc_truth_sum is None:
                mc_truth_noreco_sum = {m: truth_noreco_values[m] for m in matchings}
                mc_truth_noreco_variance_sum = {m: truth_noreco_variances[m] for m in matchings}
                mc_truth_sum = {m: truth_values[m] for m in matchings}
                mc_truth_variance_sum = {m: truth_variances[m] for m in matchings}
                theory_lower_sum = {m: theory_band_low[m] for m in matchings}
                theory_upper_sum = {m: theory_band_hi[m] for m in matchings}
                unfolding_sum = {m: unfolding_values[m] for m in matchings}
                unfolding_variance_sum = {m: unfolding_variances[m] for m in matchings}
            else:
                for matching in matchings:
                    mc_truth_noreco_sum[matching] += truth_noreco_values[matching]
                    mc_truth_noreco_variance_sum[matching] += truth_noreco_variances[matching]
                    mc_truth_sum[matching] += truth_values[matching]
                    mc_truth_variance_sum[matching] += truth_variances[matching]
                    theory_lower_sum[matching] += theory_band_low[matching]
                    theory_upper_sum[matching] += theory_band_hi[matching]
                    unfolding_sum[matching] += unfolding_values[matching]
                    unfolding_variance_sum[matching] += unfolding_variances[matching]

        for ax_ in [ax, ax_all]:
            alpha = 1.0-0.2*ipt if ax_ == ax_all else 1.0
            if plot_truth:

                for matching in matchings:
                    # hep.histplot(
                    #     truth_noreco_values[matching],
                    #     msd_edges_,
                    #     yerr=np.sqrt(truth_noreco_variances[matching]),
                    #     label="Simulation (no reco-sel){}".format(pt_bin_tex),
                    #     ax=ax_,
                    #     alpha=alpha,
                    #     color="tab:grey",
                    #     ls="-" if matching == "matching" else "--"
                    # )
                    hep.histplot(
                        truth_values[matching],
                        msd_edges_,
                        yerr=np.sqrt(truth_variances[matching]),
                        label="Simulation{}{}".format(matching_str[matching], pt_bin_tex),
                        ax=ax_,
                        ls="-",
                        alpha=alpha,
                        **matching_kwargs[matching]
                    )

                    for ibin in range(len(msd_edges_)-1):
                        ax_.fill_between(
                            msd_edges_[ibin:ibin+2],
                            (truth_values[matching]-np.sqrt(theory_band_low[matching]))[ibin],
                            (truth_values[matching]+np.sqrt(theory_band_hi[matching]))[ibin],
                            alpha=0.2,
                            label="Theory uncertainties{}".format(matching_str[matching]) if ibin == 0 else None,
                            **matching_kwargs[matching]
                        )
            for matching in matchings:
                ax_.errorbar(
                    msd_centers,
                    unfolding_values[matching],
                    yerr=np.sqrt(unfolding_variances[matching]),
                    xerr=msd_xerr,
                    label=data_label + " " + pt_bin_tex,
                    color="k",
                    alpha=alpha,
                    fmt="o",
                    markersize=6,
                    **marker_kwargs[matching]
                )
        cms_label(exp_label=exp_label, year=year, ax=ax, fs=20, data=data)

        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        ax.set_xlim(msd_min, msd_max)
        ax.set_xticks([30.0, 50.0, 100.0, 150.0, 200.0, 250.0])
        ax.set_ylim(0, ymaxs[ipt])
        ax.legend(fontsize=legend_fontsize)
        f.savefig(f"{out_dir}/m_unfold_pt{ipt}{region_str}.pdf", bbox_inches="tight")
        del f, ax

    ax_all.set_yscale("log")
    ax_all.set_ylabel(y_label)
    ax_all.set_xlabel(x_label)
    ax_all.set_xlim(msd_min, msd_max)
    ax_all.set_xticks([30.0, 50.0, 100.0, 150.0, 200.0, 250.0])
    ax_all.legend(fontsize=legend_fontsize)
    cms_label(exp_label=exp_label, year=year, fs=20, ax=ax_all, data=data)
    f_all.savefig(f"{out_dir}/m_unfold_pt_all{region_str}.pdf", bbox_inches="tight")
    plt.close()
    del ax_all, f_all
    f, ax = plt.subplots(figsize=(9, 9))
    if plot_truth:
        # hep.histplot(
        #     mc_truth_noreco_sum,
        #     msd_edges_,
        #     yerr=np.sqrt(mc_truth_noreco_variance_sum),
        #     ax=ax,
        #     label="Simulation (no reco-sel)",
        #     alpha=0.8,
        #     color="tab:grey",
        # )
        for matching in matchings:
            hep.histplot(
                mc_truth_sum[matching],
                msd_edges_,
                yerr=np.sqrt(mc_truth_variance_sum[matching]),
                ax=ax,
                label="Simulation {}".format(matching_str[matching]),
                alpha=0.8,
                ls="-",
                **matching_kwargs[matching]
            )
            for ibin in range(len(msd_edges_)-1):
                ax.fill_between(
                    msd_edges_[ibin:ibin+2],
                    (mc_truth_sum[matching]-np.sqrt(theory_lower_sum[matching]))[ibin],
                    (mc_truth_sum[matching]+np.sqrt(theory_upper_sum[matching]))[ibin],
                    alpha=0.2,
                    label="Theory uncertainties{}".format(matching_str[matching]) if ibin == 0 else None,
                    **matching_kwargs[matching]
                )
    for matching in matchings:
        ax.errorbar(
            msd_centers,
            unfolding_sum[matching],
            yerr=np.sqrt(unfolding_variance_sum[matching]),
            xerr=msd_xerr,
            label=data_label,
            color="k",
            fmt="o",
            markersize=6,
            **marker_kwargs[matching]
        )

    cms_label(exp_label=exp_label, year=year, fs=20, ax=ax, data=data)
    ax.legend(fontsize=legend_fontsize)
    ax.text(
        ax.get_xlim()[0]+0.5*np.diff(ax.get_xlim()),
        ax.get_ylim()[1]*0.6,
        # r"$650 \leq p_{T,\mathrm{truth}} < 1200~$GeV ",
        r"$p_{T,\mathrm{truth}} > 500~$GeV ",
        fontsize=20
    )

    ax.set_xlim(msd_min, msd_max)
    ax.set_xticks([30.0, 50.0, 100.0, 150.0, 200.0, 250.0])

    # ax.set_yscale("log")
    # ax.set_ylim(-1, 90.0)
    ax.set_ylim(0., ymax_sum)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    f.savefig(f"{out_dir}/m_unfold_sum{region_str}.pdf", bbox_inches="tight")
    del f, ax
    plt.close()


if __name__ == "__main__":
    import argparse
    import correctionlib
    import awkward as ak
    import json
    import glob
    from utils import jms_correction_files

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("fit_dir", help="Fit-directory, holding at least fit_shapes.root and config.json")
    parser.add_argument("--year", choices=["UL16preVFP", "UL16postVFP", "UL17", "UL18", "RunII"], default="RunII")
    parser.add_argument("--mctruth", action="store_true", help="plot mc truth against unfolded result")
    parser.add_argument("--data", action="store_true", help="whether or not the fit includes data or pseudo-data.")
    parser.add_argument("--tagger", default="substructure", choices=["substructure", "particlenetDDT"])
    parser.add_argument("--migmat", action="store_true")
    parser.add_argument("--matching-comp", action="store_true")
    parser.add_argument("--no-matching", action="store_true")
    parser.add_argument("--skip-templates", action="store_true")
    parser.add_argument("--n2cut", default="", choices=["", "n2_0p147", "n2_0p17", "n2_0p25", "n2_0p2"])
    parser.add_argument("--skipmunfold", action="store_true")

    args = parser.parse_args()
    exp_label = "Work in progress"
    if not args.skip_templates:
        plot_templates(args.fit_dir, year=args.year, exp_label=exp_label, data=args.data)
        plot_templates(args.fit_dir, year=args.year, exp_label=exp_label, data=args.data, region="fail")
    if not args.skipmunfold:
        # for yaxis in ["nevents", "au", "sigma"]:
        for yaxis in ["sigma"]:
            # print(f"plotting munfold ({yaxis})")
            # plot_unfolded_mass(
            #     args.fit_dir,
            #     data=args.data,
            #     plot_truth=args.mctruth,
            #     year=args.year,
            #     exp_label=exp_label,
            #     yaxis=yaxis,
            #     inclusive_tagger=False,
            #     plot_no_matching=args.matching_comp,
            # )
            print(f"plotting munfold ({yaxis}) inclusive tagger region")
            plot_unfolded_mass(
                args.fit_dir,
                data=args.data,
                plot_truth=args.mctruth,
                year=args.year,
                exp_label=exp_label,
                yaxis=yaxis,
                inclusive_tagger=True,
                plot_no_matching=args.no_matching,
                plot_matching_comp=args.matching_comp,
                n2cut=args.n2cut,
            )
    if not args.migmat:
        exit(0)
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"] if args.year == "RunII" else [args.year]

    pt_reco_edges = np.array([575.0, 650.0, 725.0, 800.0, 1000.0, 1200.0])
    pt_gen_edges = np.array([500.0, 650.0, 800.0, 1200.0, np.inf])
    msd_reco_edges = np.arange(50.0, 300.0, 5)
    msd_gen_edges = np.array([30.0, 70., 80., 90., np.inf])

    migmat = hist.Hist(
        hist.axis.Variable(pt_gen_edges, name="pt_gen"),
        hist.axis.Variable(msd_gen_edges, name="msd_gen"),
        hist.axis.Variable(pt_reco_edges, name="pt_reco"),
        hist.axis.Variable(msd_reco_edges, name="msd_reco"),
        storage=hist.storage.Weight(),
    )
    polynomial_msd_correction_set = correctionlib.CorrectionSet.from_file(
        jms_correction_files["notagger"]
    )
    tagger_str = "n2ddt" if args.tagger == "substructure" else "pNetddt"
    for year in years:
        msd_corr = polynomial_msd_correction_set[f"response_g_jec_{year}"]
        sel_events_tree_fname = f"WJetsToQQ_tinyTree_{year}_notagger.parquet"
        events = ak.from_parquet(sel_events_tree_fname)
        events["pt_raw"] = events.Jets.pt[:, 0]
        events["pt"] = events.pt_raw * events.jecfactor[:, 0]
        events["ptgen"] = events.pt_gen_ak8
        events["mjet_raw"] = events.mjet
        events["mjet"] = events.mjet_raw * events.jecfactor[:, 0]
        events["mjetgen"] = events.msd_gen_ak8

        migmat.fill(
            pt_reco=events.pt,
            pt_gen=events.ptgen,
            msd_reco=events.mjet/msd_corr.evaluate(events.pt),
            msd_gen=events.mjetgen,
            weight=events.weight
        )
        del sel_events_tree_fname
    n_reco = (msd_reco_edges.shape[0] - 1) * (pt_reco_edges.shape[0] - 1)
    n_gen = (msd_gen_edges.shape[0] - 1) * (pt_gen_edges.shape[0] - 1)
    print("Condition number:", np.linalg.cond(migmat.values().reshape((n_gen, n_reco))))
    unfolding_plotting.plot_migration_matrix(
        migmat,
        f"{args.fit_dir}/migration_matrix_prefit.pdf"
    )

    fit_results = json.load(open(glob.glob(args.fit_dir + "/*fitResult.json")[0], "r"))
    r_gen = np.array(
        [v[0] for p, v in fit_results.items() if "r_ptgen" in p.lower()]
    )
    unfolding_plotting.plot_migration_matrix(migmat, f"{args.fit_dir}/migration_matrix_postfit.pdf", r_gen=r_gen)

    os.system(
        "convert -delay 100 "
        + f"{args.fit_dir}/migration_matrix_prefit.pdf "
        + f"{args.fit_dir}/migration_matrix_postfit.pdf "
        + f"{args.fit_dir}/migration_matrix.gif"
    )
