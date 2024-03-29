#!/usr/bin/env pythonJMS.sh
import hist
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import awkward as ak
import mplhep as hep
from collections.abc import Callable
from typing import Union
from copy import deepcopy
import itertools
from utils import year_alias

hep.style.use("CMS")

stability_color = "#1b9e77"
purity_color = "#d95f02"

label_tex_dict = {
    "pt": "p_{T}",
    "mjet": "m_{SD}",
}


lumi = 41.47968052876168  # fb^-1
global year
year = "UL17"
common_plot_kwargs = dict(flow="none")


def cms_label(ax, fs=20):
    hep.cms.label("Work in Progress", year=year_alias.get(year, year), ax=ax, fontsize=fs, data=False)
    # try:
    #     hep.label.exp_label(llabel="Private work (CMS simulation)", year=year, ax=ax, fontsize=fs)
    # except BaseException as e:
    #     print(e)
    #     hep._exp_label(llabel="Private work (CMS simulation)", year=year, ax=ax, fontsize=fs)


def fax(w=9, h=9):
    return plt.subplots(figsize=(w, h))


def safe_div(a, b):
    return np.divide(a, b, out=np.zeros_like(a), where=b != 0)
    # return np.divide(a, b, where=b != 0)


def migration_metric(h: hist.Hist, axis_name: str = "pt_reco", flow: bool = False):
    axes = [a.name for a in h.axes]
    if axis_name not in axes:
        raise NameError(f"Did not find {axis_name} among axes of hist!")
    axis_index = axes.index(axis_name)
    mat = h.to_numpy(flow=flow)[0]

    # we assume relevant axis is 1, so if hist-axis is not index 0 transpose matrix
    # in order to renorm along correct axis
    if axis_index != 0:
        mat = mat.T
    metric_bin_contents = []
    metric_bin_edges = h.to_numpy(flow=flow)[axis_index + 1]
    other_dim_edges = h.to_numpy(flow=flow)[(1 - axis_index) + 1]
    main_dim_max_bin = len(metric_bin_edges) - 1 - (2 if flow else 0)  # -1 for last edge and -2 for uflow oflow
    second_dim_max_bin = len(other_dim_edges) - 1 - (2 if flow else 0)
    renormed_mat = safe_div(mat, np.sum(mat, axis=1)[:, None])

    for bin_ind in range(main_dim_max_bin):
        bin_ind_x = bin_ind
        bin_ind_y = bin_ind
        bin_ind_y = bin_ind
        if flow:
            bin_ind_x = bin_ind + 1
            bin_ind_y = min(bin_ind + 1, second_dim_max_bin + 1)  # set to overflow of secondary dim if bin_ind exceeds
            if bin_ind_y != bin_ind + 1:
                print("Warning! reached end of binning in secondary dimension! taking counts from overflow for renorm")
        metric_bin_contents.append(renormed_mat[bin_ind_x, bin_ind_y])

    # return calculated metric and bin edges without flow
    return (
        np.array(metric_bin_contents),
        metric_bin_edges[1:-1] if flow else metric_bin_edges,
    )


class Unfolding1DPlotter(object):
    def __init__(self, variable: str = "mjet", reco_correction: Union[float, Callable] = 1.0, flow: bool = False):
        self.reco_correction = reco_correction
        self.variable = variable
        self.flow = flow

    def reco_corr_eval(self, events: ak.Array):
        if isinstance(self.reco_correction, Callable):
            return self.reco_correction(events["pt"])
        elif isinstance(self.reco_correction, float):
            return self.reco_correction
        else:
            raise TypeError("reco mass correction factors has to be either callable or float!")

    def build_migration_matrix(
        self,
        reco_bins: hist.axis,
        gen_bins: hist.axis,
        events: ak.Array,
    ) -> hist.Hist:
        h = hist.Hist(reco_bins, gen_bins, storage=hist.storage.Weight())
        h.fill(
            **{
                reco_bins.name: events[self.variable] * self.reco_corr_eval(events),
                gen_bins.name: events[f"{self.variable}gen"],
                "weight": events.weight,
            }
        )
        return h

    def plot_migration_metric(
        self, binning: np.ndarray, events: ak.Array, w: int = 9, h: int = 9, threshold: float = None
    ):
        f, ax = fax(w, h)

        migmat = self.build_migration_matrix(
            hist.axis.Variable(binning, name=f"{self.variable}_reco", overflow=self.flow),
            hist.axis.Variable(binning, name=f"{self.variable}_gen", overflow=self.flow),
            events,
        )

        stability, metric_edges = migration_metric(migmat, f"{self.variable}_gen", flow=self.flow)

        def regular_tick(e):
            if e > 1e4:
                return r"$\infty$"
            else:
                return str(int(e))
        metric_ticks = [regular_tick(e) for e in metric_edges]
        if metric_edges[-1] > 1e4:
            metric_edges[-1] = metric_edges[-2]+np.max(np.diff(metric_edges[:-1]))

        hep.histplot(
            (stability, metric_edges),
            label="stability",
            ax=ax,
            **{"color": stability_color},
            **common_plot_kwargs,
        )
        hep.histplot(
            (migration_metric(migmat, f"{self.variable}_reco", flow=self.flow)[0], metric_edges),
            label="purity",
            ax=ax,
            **{"color": purity_color},
            **common_plot_kwargs,
        )
        if threshold is not None:
            ax.plot(ax.get_xlim(), [threshold, threshold], "k--", alpha=0.6)

        ax.legend(loc="lower center")
        ax.set_xticks(metric_edges, metric_ticks)
        ax.set_xlabel(r"$%s~$[GeV]" % label_tex_dict[self.variable])
        ax.set_ylabel("binning metric")
        cms_label(ax, fs=18)

        return f, ax

    def plot_migration_metric_nbins_comparison(
        self,
        xmin: float,
        xmax: float,
        events: ak.Array,
        nbins_list: list[int] = [100, 50],
        w: int = 12,
        h: int = 12,
    ):
        f, ax = fax(w, h)
        for nbins in nbins_list:
            migmat = self.build_migration_matrix(
                hist.axis.Regular(nbins, xmin, xmax, name=f"{self.variable}_reco", overflow=self.flow),
                hist.axis.Regular(nbins, xmin, xmax, name=f"{self.variable}_gen", overflow=self.flow),
                events,
            )
            hep.histplot(
                migration_metric(migmat, f"{self.variable}_gen", flow=self.flow),
                label=f"stability nbins={nbins}",
                ax=ax,
                **common_plot_kwargs,
            )
            hep.histplot(
                migration_metric(migmat, f"{self.variable}_reco", flow=self.flow), label=f"purity nbins={nbins}", ax=ax,
                **common_plot_kwargs,
            )
        ax.legend()
        cms_label(ax, fs=25)
        ax.set_xlabel(r"$%s~$[GeV]" % label_tex_dict[self.variable])
        ax.set_ylabel("binning metric")

        return f, ax

    def plot_migration_matrix(
        self,
        reco_binning: np.ndarray,
        gen_binning: np.ndarray,
        events: ak.Array,
        w: int = 9,
        h: int = 9,
    ):
        migmat = self.build_migration_matrix(
            hist.axis.Variable(reco_binning, name=f"{self.variable}_reco", overflow=self.flow),
            hist.axis.Variable(gen_binning, name=f"{self.variable}_gen", overflow=self.flow),
            events,
        )

        f, ax = fax()
        hep.hist2dplot(migmat.to_numpy(), ax=ax)

        cms_label(ax)
        ax.set_xlabel(r"$%s^\mathrm{reco}}~$[GeV]" % label_tex_dict[self.variable])
        ax.set_ylabel(r"$%s^\mathrm{gen}}~$[GeV]" % label_tex_dict[self.variable])
        return f, ax

    def plot_distributions(
        self,
        reco_binning: np.ndarray,
        gen_binning: np.ndarray,
        events: ak.Array,
        w: int = 9,
        h: int = 9,
    ):
        f, ax = fax(w, h)

        migmat = self.build_migration_matrix(
            hist.axis.Variable(reco_binning, name=f"{self.variable}_reco", overflow=self.flow),
            hist.axis.Variable(gen_binning, name=f"{self.variable}_gen", overflow=self.flow),
            events,
        )

        hep.histplot(
            migmat[::sum, :], ax=ax, label=r"$%s^\mathrm{ptcl}}$" % label_tex_dict[self.variable],
            **common_plot_kwargs,
        )
        hep.histplot(
            migmat[:, ::sum], ax=ax, label=r"$%s^\mathrm{reco}}$" % label_tex_dict[self.variable],
            **common_plot_kwargs,
        )
        ax.legend()
        ax.set_xlabel(r"$%s~$[GeV]" % label_tex_dict[self.variable])
        cms_label(ax)
        return f, ax


def setup_edges(x: np.ndarray, y: np.ndarray, nth_subtick: int = 1):
    edges = x
    edges = np.array(
        list(map(lambda x: round(x, 2), list(edges[:-1:nth_subtick]) * len(y[:-2]) + list(edges[::nth_subtick]))),
        # -2: -1 to remove last y edge and -1 to remove last another one to add complete row of x edges
        dtype=str,
    )
    return edges


def plot_migration_matrix_old(
        h2d_hist: hist.Hist,
        zlog: bool = False,
        output_name: str = "migration_matrix.pdf",
        extratext: str = "",
):
    h2d = h2d_hist.to_numpy()
    gen_pt_edges = h2d[1]
    gen_pt_edges[-1] = 2 * gen_pt_edges[-2] - gen_pt_edges[-3]
    gen_msd_edges = h2d[2]
    gen_msd_edges[-1] = 2 * gen_msd_edges[-2] - gen_msd_edges[-3]

    reco_pt_edges = h2d[3]
    reco_pt_edges[-1] = 2 * reco_pt_edges[-2] - reco_pt_edges[-3]
    reco_msd_edges = h2d[4]

    edges_generator = setup_edges(gen_msd_edges, gen_pt_edges)
    edges_detector = setup_edges(reco_msd_edges, reco_pt_edges)

    x_shape = h2d[0].shape[0] * h2d[0].shape[1]
    y_shape = h2d[0].shape[2] * h2d[0].shape[3]
    mat = h2d[0].reshape((x_shape, y_shape)).T  # gen-bins on x-axis

    fig, (ax, cax) = plt.subplots(ncols=2, figsize=(12, 12), gridspec_kw={"width_ratios": [1, 0.02]})

    ax_pt_generator = ax.twiny()
    ax_pt_detector = ax.twinx()

    # colorbar
    n_tick_msd = 3
    every_nth_tick = int(len(gen_msd_edges) / n_tick_msd)
    generator_labels = setup_edges(gen_msd_edges, gen_pt_edges, every_nth_tick)
    msd_tick_position_generator = np.linspace(0, x_shape, len(generator_labels))
    n_tick_msd = 3
    every_nth_tick = int(len(reco_msd_edges) / n_tick_msd)
    detector_labels = setup_edges(reco_msd_edges, reco_pt_edges, every_nth_tick)
    msd_tick_position_detector = np.linspace(0, y_shape, len(detector_labels))

    sns.heatmap(mat, ax=ax, cmap="magma", norm=LogNorm() if zlog else None, cbar=False)

    fig.colorbar(ax.get_children()[0], cax=cax, orientation="vertical", label="Events")

    ax.invert_yaxis()
    ax_pt_detector.invert_yaxis()

    n_reco_bins = len(reco_msd_edges) - 1
    n_gen_bins = len(gen_msd_edges) - 1

    def diagonal_line(ax, starting_bin, ending_bin):
        ax.plot(
            [
                starting_bin[0] * n_gen_bins,
                (ending_bin[0] + 1) * n_gen_bins,
            ],
            [
                starting_bin[1] * n_reco_bins,
                (ending_bin[1] + 1) * n_reco_bins,
            ],
            "k--",
            alpha=0.6,
        )

    for i in range(1, len(reco_pt_edges) - 1):
        min_ = 0
        max_ = len(edges_generator)
        span = len(reco_msd_edges) - 1
        span_arr = [span * i, span * i]
        ax.plot([min_, max_], span_arr, "k--", alpha=0.6)
        # if draw_diagonal:
        #     i_0 = 0
        #     j_0 = i
        #     max_diff_i = len(gen_pt_edges) - 2 - i_0
        #     max_diff_j = len(reco_pt_edges) - 2 - j_0
        #     if max_diff_i < max_diff_j:
        #         i_f = len(gen_pt_edges) - 2
        #         j_f = j_0 - i_0 + i_f
        #     else:
        #         j_f = len(reco_pt_edges) - 2
        #         i_f = i_0 - j_0 + j_f
        #     diagonal_line(ax, (i_0, j_0), (i_f, j_f))

    for i in range(0, len(gen_pt_edges) - 1):
        # horizontal lines
        if i > 0:
            min_ = 0
            max_ = len(edges_detector)
            span = len(gen_msd_edges) - 1
            span_arr = [span * i, span * i]
            ax.plot(span_arr, [min_, max_], "k--", alpha=0.6)

        # if draw_diagonal:
        #     # diagonal lines lower half
        #     i_0 = i
        #     j_0 = 0
        #     max_diff_i = len(gen_pt_edges) - 2 - i_0
        #     max_diff_j = len(reco_pt_edges) - 2 - j_0
        #     if max_diff_i < max_diff_j:
        #         i_f = len(gen_pt_edges) - 2
        #         j_f = j_0 - i_0 + i_f
        #     else:
        #         j_f = len(reco_pt_edges) - 2
        #         i_f = i_0 - j_0 + j_f
        #     diagonal_line(ax, (i_0, j_0), (i_f, j_f))

    generator_labels[-1] = r"$\infty$"
    ax.set_xticks(msd_tick_position_generator)
    ax.set_xticklabels(generator_labels)

    ax.set_yticks(msd_tick_position_detector)
    ax.set_yticklabels(detector_labels)

    pt_tick_position_generator = np.linspace(0, x_shape, len(gen_pt_edges))
    ax_pt_generator.set_xlim(ax.get_xlim())
    ax_pt_generator.set_xticks(pt_tick_position_generator)
    gen_pt_labels = gen_pt_edges.astype(dtype=str)
    gen_pt_labels[-1] = r"$\infty$"
    ax_pt_generator.set_xticklabels(gen_pt_labels)
    ax_pt_generator.set_xlabel("$p_{T,gen}$ [GeV]")

    pt_tick_position_detector = np.linspace(0, y_shape, len(reco_pt_edges))
    ax_pt_detector.set_ylim(ax.get_ylim())
    ax_pt_detector.set_yticks(pt_tick_position_detector)
    reco_pt_labels = reco_pt_edges.astype(dtype=str)
    reco_pt_labels[-1] = r"$\infty$"
    ax_pt_detector.set_yticklabels(reco_pt_labels)
    ax_pt_detector.set_ylabel("$p_{T,reco}$ [GeV]")

    ax.set_xlabel("$m_{SD,gen}$ [GeV]")
    ax.set_ylabel("$m_{SD,reco}$ [GeV]")

    hep.cms.text("Simulation, Work in progress", ax=ax, fontsize=23, pad=0.03)
    ax.text(0.05 * x_shape, 0.9 * y_shape, r"$W(q\bar{q})$+jets")
    ax.text(0.05 * x_shape, 0.85 * y_shape, extratext)
    fig.tight_layout()
    plt.savefig(
        output_name,
        bbox_inches="tight",
        pad_inches=0.01,
    )


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


def plot_migration_matrix(
    migmat: hist.Hist,
    outname: str,
    r_gen=1.0,
    extratext: str = "",
):
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
    ax.text(0.05 * ax.get_xlim()[1], 0.9 * ax.get_ylim()[1], r"$W(q\bar{q})$+jets")
    ax.text(0.05 * ax.get_xlim()[1], 0.85 * ax.get_ylim()[1], extratext)

    f.colorbar(c, cax=cax, orientation="vertical", label="Events")
    f.savefig(
        outname,
        bbox_inches="tight"
    )


def plot_efficiencies(h, outdir):
    f, ax_ = plt.subplots(2, 3, figsize=(3 * 9, 2 * 9))
    pt_gen_tex = [
        r"500 < $p_{T,\mathrm{gen}}$ < 650 GeV",
        r"650 < $p_{T,\mathrm{gen}}$ < 800 GeV",
        r"800 < $p_{T,\mathrm{gen}}$ < 1200 GeV",
        r"$p_{T,\mathrm{gen}}$ > 1200 GeV",
    ]

    region = "inclusive"
    for i in [0, 1, 2, 3]:
        ax = ax_.flatten()[i]
        msd_mat = deepcopy(h[{"pt_gen": i}].integrate("pt_reco", 575j))
        h1d = msd_mat.integrate("msd_gen", 10j, 70j)
        h1d_region_0 = h1d[{"msd_reco": slice(0j, 50j)}]
        h1d_region = h1d[{"msd_reco": slice(50j, 70j)}]
        h1d_region_2 = h1d[{"msd_reco": slice(70j, None)}]
        h1d_sum = h1d.sum().value
        h1d_region_sum = h1d_region.sum().value
        h1d_region_2_sum = h1d_region_2.sum().value
        h1d_region_0_sum = h1d_region_0.sum().value
        hep.histplot(h1d, ax=ax, label=r"$10<m_\mathrm{SD, gen} < 70$ [GeV] -> total", flow="none")
        hep.histplot(
            h1d_region_0,
            ax=ax,
            histtype="fill",
            color="tab:red",
            alpha=0.3,
            label=r"$m_\mathrm{SD, reco} < 50~\mathrm{[GeV]}$"
            + f" -> {100.*(h1d_region_0_sum/h1d_sum):.2f}% of total",
            flow="none",
        )
        hep.histplot(
            h1d_region,
            ax=ax,
            histtype="fill",
            color="k",
            alpha=0.3,
            label=r"$50<m_\mathrm{SD, reco} < 70~\mathrm{[GeV]}$"
            + f" -> {100.*(h1d_region_sum/h1d_sum):.2f}% of total",
            flow="none",
        )
        hep.histplot(
            h1d_region_2,
            ax=ax,
            histtype="fill",
            color="tab:orange",
            alpha=0.3,
            label=r"$m_\mathrm{SD, reco} > 70~\mathrm{[GeV]}$"
            + f" -> {100.*(h1d_region_2_sum/h1d_sum):.2f}% of total",
            flow="none",
        )
        hep.cms.label("Work in Progress", fontsize=16, ax=ax, year=year)
        ax.legend(fontsize=15)
        ax.set_xlim(0, 150)
        # ax.set_yscale("log")
        ymax = np.max(h1d.values()) * (1.0 / 0.6)
        ymin = ax.get_ylim()[0]

        # ax.text(50, 10**((np.ceil(np.log10(ymin))+np.ceil(np.log10(ymax)))*0.5), pt_gen_tex[i], fontsize=16)
        ax.text(50, ymax*0.65, pt_gen_tex[i], fontsize=16)

        ax.set_ylim(None, ymax)
    f.savefig(f"{outdir}/msd_bin_gen_comparison_ptgen_{year}_{region}.pdf")

    f, ax_ = plt.subplots(2, 3, figsize=(3 * 9, 2 * 9))
    f1, ax1_ = plt.subplots(2, 3, figsize=(3 * 9, 2 * 9))
    pt_reco_tex = [
        r"575 < $p_{T,\mathrm{reco}}$ < 650 GeV",
        r"650 < $p_{T,\mathrm{reco}}$ < 725 GeV",
        r"725 < $p_{T,\mathrm{reco}}$ < 800 GeV",
        r"800 < $p_{T,\mathrm{reco}}$ < 1000 GeV",
        r"1000 < $p_{T,\mathrm{reco}}$ < 1200 GeV",
        r"$p_{T,\mathrm{reco}}$ > 1200 GeV",
    ]
    for i in [0, 1, 2, 3, 4, 5]:
        ax = ax_.flatten()[i]
        ax1 = ax1_.flatten()[i]
        msd_mat = h[{"pt_reco": i}].integrate("pt_gen")
        h1d = deepcopy(msd_mat[{"msd_reco": slice(50j, 300j)}]).integrate("msd_gen", 10j, 300j)
        h1d_sum = h1d.sum().value
        h1d_region_0 = deepcopy(msd_mat[{"msd_reco": slice(50j, 300j)}]).integrate("msd_gen", 10j, 50j)
        h1d_region_0_sum = h1d_region_0.sum().value
        h1d_region_1 = deepcopy(msd_mat[{"msd_reco": slice(50j, 300j)}]).integrate("msd_gen", 50j, None)
        h1d_region_1_sum = h1d_region_1.sum().value

        hep.histplot(h1d, ax=ax, label=r"$50<m_\mathrm{SD, reco} < 300$ [GeV] -> 100% of total", flow="none")
        hists = [h1d_region_0, h1d_region_1]
        labels = [
            r"$10< m_\mathrm{SD, gen} < 50~\mathrm{[GeV]}$" + f" -> {100.*(h1d_region_0_sum/h1d_sum):.2f}% of total",
            r"$m_\mathrm{SD, gen} > 50~\mathrm{[GeV]}$" + f" -> {100.*(h1d_region_1_sum/h1d_sum):.2f}% of total",
        ]
        colors = ["k", "tab:green"]
        hep.histplot(hists, ax=ax, histtype="fill", stack=True, color=colors, alpha=0.3, label=labels, flow="none")
        hep.cms.label("Work in Progress", fontsize=16, ax=ax, year=year)
        ax.set_yscale("log")
        ymax = ax.get_ylim()[1]
        ymin = ax.get_ylim()[0]
        ax.text(50, 10**((np.ceil(np.log10(ymin))+np.ceil(np.log10(ymax)))*0.5), pt_reco_tex[i], fontsize=16)

        # ax.text(150, ax.get_ylim()[1] * 0.75, pt_reco_tex[i], fontsize=16)

        ax.legend(fontsize=15)

        h1d = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 10j, 300j)
        h1d_sum = h1d.sum().value
        h1d_region_0 = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 10j, 50j)
        h1d_region_0_sum = h1d_region_0.sum().value
        h1d_region_1 = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 50j, 300j)
        h1d_region_1_sum = h1d_region_1.sum().value
        labels_ratios = [
            r"fraction of events with $10< m_\mathrm{SD, gen} < 50~\mathrm{[GeV]}$"
            + f" -> {100.*(h1d_region_0_sum/h1d_sum):.2f}% of total",
            r"fraction of events with $m_\mathrm{SD, gen} > 50~\mathrm{[GeV]}$"
            + f" -> {100.*(h1d_region_1_sum/h1d_sum):.2f}% of total",
        ]

        hists_ratios = [(h.values() / h1d.values(), h1d.axes[0].edges) for h in [h1d_region_0, h1d_region_1]]
        hep.histplot(
            hists_ratios, ax=ax1, histtype="fill", stack=True, color=colors, alpha=0.3, label=labels_ratios, flow="none"
        )
        hep.cms.label("Work in Progress", fontsize=16, ax=ax1, year=year)
        ax1.text(150, 0.75, pt_reco_tex[i], fontsize=14)
        ax1.legend(fontsize=15)
        ax1.set_ylim(0, 1.1)

    f.savefig(f"{outdir}/msd_bin_reco_comparison_ptreco_{year}_{region}.pdf")
    f1.savefig(f"{outdir}/msd_bin_reco_fractions_ptreco_{year}_{region}.pdf")

    # same but with gen <-> reco swapped

    f, ax_ = plt.subplots(2, 3, figsize=(3 * 9, 2 * 9))
    f1, ax1_ = plt.subplots(2, 3, figsize=(3 * 9, 2 * 9))
    pt_reco_tex = [
        r"575 < $p_{T,\mathrm{reco}}$ < 650 GeV",
        r"650 < $p_{T,\mathrm{reco}}$ < 725 GeV",
        r"725 < $p_{T,\mathrm{reco}}$ < 800 GeV",
        r"800 < $p_{T,\mathrm{reco}}$ < 1000 GeV",
        r"1000 < $p_{T,\mathrm{reco}}$ < 1200 GeV",
        r"$p_{T,\mathrm{reco}}$ > 1200 GeV",
    ]
    for i in [0, 1, 2, 3]:
        ax = ax_.flatten()[i]
        # ax1 = ax1_.flatten()[i]
        msd_mat = h[{"pt_gen": i}].integrate("pt_reco")
        h1d = deepcopy(msd_mat[{"msd_gen": slice(10j, 300j)}]).integrate("msd_reco", 50j, 300j)
        h1d_sum = h1d.sum().value
        h1d_region_0 = deepcopy(msd_mat[{"msd_gen": slice(10j, 300j)}]).integrate("msd_reco", 50j, 70j)
        h1d_region_0_sum = h1d_region_0.sum().value
        h1d_region_1 = deepcopy(msd_mat[{"msd_gen": slice(10j, 300j)}]).integrate("msd_reco", 70j, 300j)
        h1d_region_1_sum = h1d_region_1.sum().value

        hep.histplot(h1d, ax=ax, label=r"$50<m_\mathrm{SD, gen} < 300$ [GeV] -> 100% of total", flow="none")
        hists = [h1d_region_0, h1d_region_1]
        labels = [
            r"$50< m_\mathrm{SD, reco} < 70~\mathrm{[GeV]}$" + f" -> {100.*(h1d_region_0_sum/h1d_sum):.2f}% of total",
            r"$70 < m_\mathrm{SD, reco} < 300~\mathrm{[GeV]}$" + f" -> {100.*(h1d_region_1_sum/h1d_sum):.2f}% of total",
        ]
        colors = ["k", "tab:green"]
        hep.histplot(hists, ax=ax, histtype="fill", stack=True, color=colors, alpha=0.3, label=labels, flow="none")
        hep.cms.label("Work in Progress", fontsize=16, ax=ax, year=year)
        ax.set_yscale("log")
        ymax = ax.get_ylim()[1]
        ymin = ax.get_ylim()[0]
        ax.text(50, 10**((np.ceil(np.log10(ymin))+np.ceil(np.log10(ymax)))*0.5), pt_reco_tex[i], fontsize=16)

        # ax.text(150, ax.get_ylim()[1] * 0.75, pt_reco_tex[i], fontsize=16)

        ax.legend(fontsize=15)

        # h1d = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 10j, 300j)
        # h1d_sum = h1d.sum().value
        # h1d_region_0 = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 10j, 50j)
        # h1d_region_0_sum = h1d_region_0.sum().value
        # h1d_region_1 = deepcopy(msd_mat[{"msd_reco": slice(0j, 300j)}]).integrate("msd_gen", 50j, 300j)
        # h1d_region_1_sum = h1d_region_1.sum().value
        # labels_ratios = [
        #     r"fraction of events with $10< m_\mathrm{SD, gen} < 50~\mathrm{[GeV]}$"
        #     + f" -> {100.*(h1d_region_0_sum/h1d_sum):.2f}% of total",
        #     r"fraction of events with $m_\mathrm{SD, gen} > 50~\mathrm{[GeV]}$"
        #     + f" -> {100.*(h1d_region_1_sum/h1d_sum):.2f}% of total",
        # ]

        # hists_ratios = [(h.values() / h1d.values(), h1d.axes[0].edges) for h in [h1d_region_0, h1d_region_1]]
        # hep.histplot(
        #     hists_ratios, ax=ax1, histtype="fill", stack=True,
        #     color=colors, alpha=0.3, label=labels_ratios, flow="none"
        # )
        # hep.cms.label("Work in Progress", fontsize=16, ax=ax1, year=year)
        # ax1.text(150, 0.75, pt_reco_tex[i], fontsize=14)
        # ax1.legend(fontsize=15)
        # ax1.set_ylim(0, 1.1)

    f.savefig(f"{outdir}/msd_bin_gen_reco_comparison_ptgen_{year}_{region}.pdf")
    # f1.savefig(f"{outdir}/msd_bin_gen_fractions_ptgen_{year}_{region}.pdf")


if __name__ == "__main__":
    import correctionlib
    import os
    from utils import jms_correction_files
    import argparse
    # import re

    parser = argparse.ArgumentParser()

    # parser.add_argument("--input", "-i", type=str, default="WJetsToQQ_tinyTree_UL17_n2ddt.parquet")
    parser.add_argument("--year", default="RunII")
    parser.add_argument("--tagger", default="notagger", choices=["n2ddt", "pNetddt", "notagger"])
    args = parser.parse_args()

    triggersf = True

    # year = "UL17"
    # year_re = re.search("(UL)*(20)*1[678]{1}(preVFP|postVFP)*", args.input)
    # if year_re:
    #     year = year_re.group()

    workdir = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python"
    outdir = f"unfolding_plots_{args.year}_{args.tagger}"
    if not os.path.exists(f"{workdir}/{outdir}"):
        os.makedirs(f"{workdir}/{outdir}")

    def load_tree(year):
        events_sel = ak.from_parquet(f"{workdir}/WJetsToQQ_tinyTree_{year}_notagger.parquet")

        if triggersf:
            events_sel["weight"] = events_sel["weight"]*events_sel["triggersf"]

        events_sel["pt_raw"] = events_sel.Jets.pt[:, 0]
        events_sel["pt"] = events_sel.pt_raw * events_sel.jecfactor[:, 0]
        events_sel["ptgen"] = events_sel.pt_gen_ak8
        events_sel["mjet_raw"] = events_sel.mjet
        events_sel["mjet"] = events_sel.mjet_raw * events_sel.jecfactor[:, 0]
        events_sel["mjetgen"] = events_sel.msd_gen_ak8
        events_sel["rho"] = 2 * np.log(events_sel.mjet / events_sel.pt)
        events_sel["rhogen"] = 2 * np.log(events_sel.mjetgen / events_sel.ptgen)
        events_sel = events_sel[
            (events_sel.pass_reco_selection == 1)
            & (events_sel.pass_gen_selection == 1)
            & (events_sel.dR_reco_gen < 0.4)
        ]

        return events_sel

    polynomial_msd_correction_set = correctionlib.CorrectionSet.from_file(
        # f"{workdir}/{jms_correction_files[args.tagger]}"
        f"{workdir}/"+jms_correction_files["notagger"]
    )

    pt_reco_ax = hist.axis.Variable(
        np.array([575, 650, 725, 800, 1000, 1200, np.inf]), name="pt_reco", label=r"$p_{T,\mathrm{reco}}$ [GeV]"
    )
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"] if args.year == "RunII" else [args.year]

    # with correction on msd reco
    for correction in [True, False]:

        h_2d_fine = hist.Hist(
            hist.axis.Variable(
                np.array([500, 650, 800, 1200, np.inf]), name="pt_gen", label=r"$p_{T,\mathrm{gen}}$ [GeV]"
            ),
            hist.axis.Regular(54, 30.0, 300.0, name="msd_gen", label=r"$m_{SD,\mathrm{gen}}$ [GeV]"),
            pt_reco_ax,
            hist.axis.Regular(54, 30.0, 300.0, name="msd_reco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"),
            storage=hist.storage.Weight(),
        )

        h_2d_fine_large = hist.Hist(
            hist.axis.Variable(
                np.array([500, 650, 800, 1200, np.inf]), name="pt_gen", label=r"$p_{T,\mathrm{gen}}$ [GeV]"
            ),
            hist.axis.Regular(60, 0.0, 300.0, name="msd_gen", label=r"$m_{SD,\mathrm{gen}}$ [GeV]"),
            pt_reco_ax,
            hist.axis.Regular(60, 0.0, 300.0, name="msd_reco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"),
            storage=hist.storage.Weight(),
        )

        for year in years:
            events_sel = load_tree(year)

            def msd_corr(pt):
                if correction:
                    return 1.0 / polynomial_msd_correction_set[f"response_g_jec_{year}"].evaluate(pt)
                else:
                    return 1.0

            correction_str = "_msdcorrected" if correction else ""

            h_2d_fine.fill(
                **{
                    "pt_gen": events_sel.ptgen,
                    "pt_reco": events_sel.pt,
                    "msd_gen": events_sel.mjetgen,
                    "msd_reco": events_sel.mjet * msd_corr(events_sel.pt),
                    "weight": events_sel.weight,
                }
            )
            h_2d_fine_large.fill(
                **{
                    "pt_gen": events_sel.ptgen,
                    "pt_reco": events_sel.pt,
                    "msd_gen": events_sel.mjetgen,
                    "msd_reco": events_sel.mjet * msd_corr(events_sel.pt),
                    "weight": events_sel.weight,
                }
            )

            del events_sel

            plot_migration_matrix(
                h_2d_fine,
                # zlog=True,
                extratext="msd corrected" if correction else "",
                outname=f"{workdir}/{outdir}/migration_matrix_fine_binning_{args.year}{correction_str}.pdf",
            )

            plot_efficiencies(h_2d_fine_large,
                              outdir=f"{workdir}/{outdir}")
    cmd = (
        "convert -delay 100 "
        + f"{workdir}/{outdir}/migration_matrix_fine_binning_{args.year}.pdf "
        + f"{workdir}/{outdir}/migration_matrix_fine_binning_{args.year}_msdcorrected.pdf "
        + f"{workdir}/{outdir}/migration_matrix_fine_binning_{args.year}.gif"
    )
    os.system(cmd)

    # final choice of binning
    # old binning
    # mjetgen_binning = np.array([0, 48.0, 67.0, 81.0, 92.5, 134.5, np.inf])
    # old bininng from optimization on fail region
    # optimization on pass region in 575 < pT < 650 (response_g_jec JMS factor)
    # mjetgen_binning = np.array([0.0, 44.5, 68.0, 80.5, 92.0, 132.5, 300])
    # mjetgen_binning = np.array([0., 70., 80.5, 89.5, 300.0])
    mjetgen_binning = np.array([30., 70., 80., 90., np.inf])
    h_2d = hist.Hist(
        hist.axis.Variable(np.array([500, 650, 800, 1200, np.inf]), name="pt_gen", label=r"$p_{T,\mathrm{gen}}$ [GeV]"),
        hist.axis.Variable(mjetgen_binning, name="msd_gen", label=r"$m_{SD,\mathrm{gen}}$ [GeV]"),
        pt_reco_ax,
        hist.axis.Regular(54, 50.0, 300.0, name="msd_reco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"),
        storage=hist.storage.Weight(),
    )
    for year in years:
        events_sel = load_tree(year)

        def msd_corr(pt):
            return 1.0 / polynomial_msd_correction_set[f"response_g_jec_{year}"].evaluate(pt)

        h_2d.fill(
            **{
                "pt_gen": events_sel.ptgen,
                "pt_reco": events_sel.pt,
                "msd_gen": events_sel.mjetgen,
                "msd_reco": events_sel.mjet * msd_corr(events_sel.pt),
                "weight": events_sel.weight,
            }
        )
        plot_migration_matrix(
            h_2d,
            # zlog=True,
            outname=f"{workdir}/{outdir}/migration_matrix_final_binning_{args.year}.pdf"
        )

        unfolding_1d_plotter = Unfolding1DPlotter("mjet", msd_corr)
        for pt_low, pt_high in pt_reco_ax:
            pt_low_str, pt_high_str = map(lambda b: str(b).replace(".0", ""), [pt_low, pt_high])
            f, ax = unfolding_1d_plotter.plot_migration_metric(
                mjetgen_binning,
                events_sel[(events_sel.pt >= pt_low) & (events_sel.pt < pt_high)],
                # threshold=0.5,
            )
            # ax.text(
            #     1 / 3 * ax.get_xlim()[1],
            #     2 / 3 * ax.get_ylim()[1],
            #     r"$%s \leq p_T^\mathrm{reco} < %s$" % (pt_low_str, pt_high_str),
            #     fontsize=18,
            # )
            ax.plot(
                [],
                [],
                " ",
                label=r"$%s \leq p_T^\mathrm{reco} < %s$" % (pt_low_str, pt_high_str),
            )
            ax.legend(loc="lower center")
            ax.set_ylim(0.0, 1.0)
            f.savefig(
                f"{workdir}/{outdir}/"
                + f"migration_metric_final_binning_{year}_pt"
                + f"{pt_low_str}To{pt_high_str}.pdf",
                bbox_inches="tight",
            )
        del events_sel
