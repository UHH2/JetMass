import hist
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import awkward as ak
import mplhep as hep
from collections.abc import Callable
from typing import Union

hep.style.use("CMS")

stability_color = "#1b9e77"
purity_color = "#d95f02"

label_tex_dict = {
    "pt": "p_{T}",
    "mjet": "m_{SD}",
}


lumi = 41.47968052876168  # fb^-1


def cms_label(ax, fs=20):
    hep.cms.label(label=", Work in Progress", year=2017, ax=ax, fontsize=fs)


def fax(w=9, h=9):
    return plt.subplots(figsize=(w, h))


def migration_metric(h: hist.Hist, axis_name: str = "pt_reco", flow: bool = False):
    axes = [a.name for a in h.axes]
    if axis_name not in axes:
        raise NameError(f"Did not find {axis_name} among axes of hist!")
    axis_index = axes.index(axis_name)
    mat = h.to_numpy(flow=flow)[0]
    if axis_index != 0:
        mat = mat.T
    metric_bin_contents = []
    metric_bin_edges = h.to_numpy(flow=flow)[axis_index + 1]
    other_dim_edges = h.to_numpy(flow=flow)[(1 - axis_index) + 1]
    main_dim_max_bin = len(metric_bin_edges) - 1 - (2 if flow else 0)  # -1 for last edge and -2 for uflow oflow
    second_dim_max_bin = len(other_dim_edges) - 1 - (2 if flow else 0)

    renormed_mat = mat / np.sum(mat, axis=1)[:, None]

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

        hep.histplot(
            migration_metric(migmat, f"{self.variable}_reco", flow=self.flow),
            label="stability",
            ax=ax,
            **{"color": stability_color},
        )
        hep.histplot(
            migration_metric(migmat, f"{self.variable}_gen", flow=self.flow),
            label="purity",
            ax=ax,
            **{"color": purity_color},
        )
        if threshold is not None:
            ax.plot(ax.get_xlim(), [threshold, threshold], "k--", alpha=0.6)

        ax.legend()
        ax.set_xlabel(r"$%s~$[GeV]" % label_tex_dict[self.variable])
        ax.set_ylabel("binning metric")
        cms_label(ax)

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
                migration_metric(migmat, f"{self.variable}_reco", flow=self.flow),
                label=f"stability nbins={nbins}",
                ax=ax,
            )
            hep.histplot(
                migration_metric(migmat, f"{self.variable}_gen", flow=self.flow), label=f"purity nbins={nbins}", ax=ax
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

        hep.histplot(migmat[::sum, :], ax=ax, label=r"$%s^\mathrm{gen}}$" % label_tex_dict[self.variable])
        hep.histplot(migmat[:, ::sum], ax=ax, label=r"$%s^\mathrm{reco}}$" % label_tex_dict[self.variable])
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


def plot_migration_matrix(
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
    draw_diagonal = n_reco_bins == n_gen_bins

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
        if draw_diagonal:
            i_0 = 0
            j_0 = i
            max_diff_i = len(gen_pt_edges) - 2 - i_0
            max_diff_j = len(reco_pt_edges) - 2 - j_0
            if max_diff_i < max_diff_j:
                i_f = len(gen_pt_edges) - 2
                j_f = j_0 - i_0 + i_f
            else:
                j_f = len(reco_pt_edges) - 2
                i_f = i_0 - j_0 + j_f
            diagonal_line(ax, (i_0, j_0), (i_f, j_f))

    for i in range(0, len(gen_pt_edges) - 1):
        # horizontal lines
        if i > 0:
            min_ = 0
            max_ = len(edges_detector)
            span = len(gen_msd_edges) - 1
            span_arr = [span * i, span * i]
            ax.plot(span_arr, [min_, max_], "k--", alpha=0.6)

        if draw_diagonal:
            # diagonal lines lower half
            i_0 = i
            j_0 = 0
            max_diff_i = len(gen_pt_edges) - 2 - i_0
            max_diff_j = len(reco_pt_edges) - 2 - j_0
            if max_diff_i < max_diff_j:
                i_f = len(gen_pt_edges) - 2
                j_f = j_0 - i_0 + i_f
            else:
                j_f = len(reco_pt_edges) - 2
                i_f = i_0 - j_0 + j_f
            diagonal_line(ax, (i_0, j_0), (i_f, j_f))

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
    ax.text(0.05 * x_shape, 0.9 * y_shape, r"$W$(qq)+jets")
    ax.text(0.05 * x_shape, 0.85 * y_shape, extratext)
    fig.tight_layout()
    plt.savefig(
        output_name,
        bbox_inches="tight",
        pad_inches=0.01,
    )


if __name__ == "__main__":
    import correctionlib
    import os
    events_sel = ak.from_parquet("WJetsToQQ_tinyTree.parquet")

    events_sel["pt_raw"] = events_sel.Jets.pt[:, 0]
    events_sel["pt"] = events_sel.pt_raw * events_sel.jecfactor[:, 0]
    events_sel["ptgen"] = events_sel.pt_gen_ak8
    events_sel["mjet_raw"] = events_sel.mjet
    events_sel["mjet"] = events_sel.mjet_raw * events_sel.jecfactor[:, 0]
    events_sel["mjetgen"] = events_sel.msd_gen_ak8
    events_sel["rho"] = 2 * np.log(events_sel.mjet / events_sel.pt)
    events_sel["rhogen"] = 2 * np.log(events_sel.mjetgen / events_sel.ptgen)

    polynomial_msd_correction_set = correctionlib.CorrectionSet.from_file(
        "jms_corrections_quadratic_47e3e54d1c.json"
    )

    pt_reco_ax = hist.axis.Variable(
        np.array([500, 575, 650, 725, 800, 1000, 1200, np.inf]), name="ptreco", label=r"$p_{T,\mathrm{reco}}$ [GeV]"
    )

    # with correction on msd reco
    for correction in [True, False]:

        def msd_corr(pt):
            if correction:
                return 1.0 / polynomial_msd_correction_set["response_g_jec"].evaluate(pt)
            else:
                return 1.0

        correction_str = "_msdcorrected" if correction else ""
        h_2d_fine = hist.Hist(
            hist.axis.Variable(
                np.array([0, 650, 800, 1200, np.inf]), name="ptgen", label=r"$p_{T,\mathrm{gen}}$ [GeV]"
            ),
            hist.axis.Regular(54, 30.0, 300.0, name="mJgen", label=r"$m_{SD,\mathrm{gen}}$ [GeV]"),
            pt_reco_ax,
            hist.axis.Regular(54, 30.0, 300.0, name="mJreco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"),
            storage=hist.storage.Weight(),
        )

        h_2d_fine.fill(
            **{
                "ptgen": events_sel.ptgen,
                "ptreco": events_sel.pt,
                "mJgen": events_sel.mjetgen,
                "mJreco": events_sel.mjet * msd_corr(events_sel.pt),
                "weight": events_sel.weight,
            }
        )

        plot_migration_matrix(
            h_2d_fine,
            zlog=True,
            extratext="msd corrected" if correction else "",
            output_name=f"unfolding_binning_plots/migration_matrix_fine_binning{correction_str}.pdf",
        )

    cmd = (
        "convert -delay 100 "
        + "unfolding_binning_plots/migration_matrix_fine_binning.pdf "
        + "unfolding_binning_plots/migration_matrix_fine_binning_msdcorrected.pdf "
        + "unfolding_binning_plots/migration_matrix_fine_binning.gif"
    )
    os.system(cmd)

    def msd_corr(pt):
        return 1.0 / polynomial_msd_correction_set["response_g_jec"].evaluate(pt)

    # final choice of binning
    mjetgen_binning = np.array([0.0, 44.5, 68.0, 80.5, 92.0, 132.5, 300])
    # mjetgen_binning = np.array([0, 48.0, 67.0, 81.0, 92.5, 134.5, np.inf])

    h_2d = hist.Hist(
        hist.axis.Variable(np.array([0, 650, 800, 1200, np.inf]), name="ptgen", label=r"$p_{T,\mathrm{gen}}$ [GeV]"),
        hist.axis.Variable(mjetgen_binning, name="mJgen", label=r"$m_{SD,\mathrm{gen}}$ [GeV]"),
        pt_reco_ax,
        hist.axis.Regular(54, 30.0, 300.0, name="mJreco", label=r"$m_{SD,\mathrm{reco}}$ [GeV]"),
        storage=hist.storage.Weight(),
    )

    h_2d.fill(
        **{
            "ptgen": events_sel.ptgen,
            "ptreco": events_sel.pt,
            "mJgen": events_sel.mjetgen,
            "mJreco": events_sel.mjet * msd_corr(events_sel.pt),
            "weight": events_sel.weight,
        }
    )

    plot_migration_matrix(h_2d, zlog=True, output_name="unfolding_binning_plots/migration_matrix_final_binning.pdf")

    unfolding_1d_plotter = Unfolding1DPlotter("mjet", msd_corr)
    for pt_low, pt_high in pt_reco_ax:
        pt_low_str, pt_high_str = map(lambda b: str(b).replace(".0", ""), [pt_low, pt_high])
        f, ax = unfolding_1d_plotter.plot_migration_metric(
            mjetgen_binning,
            events_sel[(events_sel.pt >= pt_low) & (events_sel.pt < pt_high)],
            threshold=0.5,
        )
        ax.text(
            1 / 3 * ax.get_xlim()[1],
            2 / 3 * ax.get_ylim()[1],
            r"$%s \leq p_T^\mathrm{reco} < %s$" % (pt_low_str, pt_high_str),
            fontsize=18,
        )
        ax.set_ylim(0.0, 1.0)
        f.savefig(
            "unfolding_binning_plots/" + "migration_metric_final_binning_pt" + f"{pt_low_str}To{pt_high_str}.pdf",
            bbox_inches="tight",
        )
