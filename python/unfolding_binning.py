import awkward as ak
import numpy as np
import hist
from unfolding_plotting import Unfolding1DPlotter, migration_metric
from JMS_from_MC import iterative_fit
from collections.abc import Callable
from typing import Union
from typing_extensions import Literal
import matplotlib.pyplot as plt
import mplhep as hep
import seaborn as sns
from matplotlib.colors import LogNorm
import unfolding_plotting
import correctionlib

hep.style.use("CMS")

logging_styles = {
    "green_bold": lambda x: f"\033[92m\033[1m{x}\033[0m",
    "orange_bold": lambda x: f"\033[93m\033[1m{x}\033[0m",
    "red_bold_underline": lambda x: f"\033[91m\033[1m\033[4m{x}\033[0m",
}


def reco_binning(gen_edges: np.ndarray) -> np.ndarray:
    """For now just return binning finer by a factor of 2, i.e. split every bin in half"""
    gen_edges_ = gen_edges[:-1] if gen_edges[-1] == np.inf else gen_edges
    reco_edges = np.concatenate((gen_edges_, (gen_edges_[:-1] + gen_edges_[1:]) / 2))
    if gen_edges[-1] == np.inf:
        reco_edges = np.concatenate((reco_edges, np.array([np.inf])))
    reco_edges.sort()
    return reco_edges


class BinOptimizer(object):
    def __init__(
        self,
        events: ak.Array,
        initial_binning: np.ndarray,
        threshold: float,
        variable: str = "mjet",
        plot_dir: str = ".",
    ):
        self.logger = logging.getLogger()
        self.events = events
        self.initial_binning = initial_binning
        self.threshold = threshold
        self.variable = variable
        self.plot_dir = plot_dir

    def save_plot(self, fax, outname):
        self.logger.info(("saving " + self.plot_dir + "/" + outname))
        fax[0].savefig((self.plot_dir + "/" + outname), bbox_inches="tight")

    def create_new_binning_(
        self,
        binning: np.ndarray,
        events: ak.Array = None,
        flow: bool = False,
        bin_merger: Literal["sequential", "bidirectional"] = "sequential",
        split_value: float = 80.0,
        mjet_reco_correction: Union[float, Callable] = 1.0,
    ):
        """
        Creating new iteration of binning in optimization process.
        Parameters:
        binning: np.ndarray
        np.ndarray initial binning to start with
        flow: bool
        whether or not to use flow bins in the hist.Hist in backend
        bin_merger: Literal["sequential", "bidirectional"] = "sequential"
        split_value: float
        if `bin_merger == 'bidirectional'` start at `split_value` and merge bins sequentially
        left and right of `split_value`
        mjet_reco_correction: Union[float, Callable]
        mjet correction factor. Either flat factor or callable. The latter will be evaluated on `self.events['pt']`
        """
        if events is None:
            events = self.events

        unfolding_plotter = Unfolding1DPlotter(self.variable, mjet_reco_correction)
        migmat = unfolding_plotter.build_migration_matrix(
            hist.axis.Variable(binning, name=f"{self.variable}_reco", overflow=flow),
            hist.axis.Variable(binning, name=f"{self.variable}_gen", overflow=flow),
            events,
        )

        stability = migration_metric(migmat, f"{self.variable}_reco", flow=flow)[0]
        purity = migration_metric(migmat, f"{self.variable}_gen", flow=flow)[0]

        if bin_merger == "sequential":
            for bin_ind in range(len(binning) - 1):
                if purity[bin_ind] < self.threshold or stability[bin_ind] < self.threshold:
                    return np.delete(binning, bin_ind + 1) if bin_ind < len(binning) - 2 else binning
        elif bin_merger == "bidirectional":
            split_bin = np.abs(binning - split_value).argmin()

            for bin_ind in reversed(range(split_bin + 1)):
                if (purity[bin_ind - 1] < self.threshold or stability[bin_ind - 1] < self.threshold) and bin_ind > 1:
                    return np.delete(binning, bin_ind - 1)  # if bin_ind == 0 else binning
                for bin_ind in range(split_bin, len(binning) - 1):
                    if (purity[bin_ind] < self.threshold or stability[bin_ind] < self.threshold) and bin_ind < len(
                        binning
                    ) - 2:
                        return np.delete(binning, bin_ind + 1)  # if bin_ind < len(binning)-2 else binning
        return binning

    def optimize_binning_(
        self,
        optimizer: Callable,
    ):
        binning = self.initial_binning  # np.arange(500,1500,5)
        # print(binning)
        counter = 0
        while counter < 1000:
            counter += 1
            n_old_bins = len(binning)
            binning = optimizer(
                {
                    "binning": binning,
                }
            )
            n_new_bins = len(binning)
            if counter % 50 == 0:
                self.logger.info(f"{n_old_bins} -> {n_new_bins}")
                # print(binning)
            if n_old_bins == n_new_bins:
                break
        self.logger.info(f"Final binning: {binning}")
        self.logger.info(f"bin widths: {(binning[1:] - binning[:-1])}")
        return binning

    def optimize_binning(
        self,
        events: ak.Array = None,
        mjet_reco_correction: Union[float, Callable] = 1.0,
        bin_merger: Literal["sequential", "bidirectional"] = "sequential",
        split_value: float = 80.0,
        flow: bool = False,
        plots: bool = False,
    ):
        if events is None:
            events = self.events

        if split_value < 0:
            unfolding_plotter = Unfolding1DPlotter(self.variable, mjet_reco_correction)
            migmat = unfolding_plotter.build_migration_matrix(
                hist.axis.Variable(self.initial_binning, name=f"{self.variable}_reco", overflow=flow),
                hist.axis.Variable(self.initial_binning, name=f"{self.variable}_gen", overflow=flow),
                events,
            )

            # postfit parameters of last iteration
            postfit_parameters = iterative_fit(migmat[::sum, :])[0][-1]["popt"]
            # hadrcoded: second to last parameter is always mean
            split_value = postfit_parameters[-2]

        # create optimizer with or without correction_factor
        def optimizer(x):
            return self.create_new_binning_(
                **x,
                events=events,
                flow=flow,
                bin_merger=bin_merger,
                split_value=split_value,
                mjet_reco_correction=mjet_reco_correction,
            )

        # optimize initial mjetgen bining and create mjetreco binning
        final_binning = self.optimize_binning_(optimizer)
        optimization_plot = None
        if plots:
            optimization_plot = self.make_optimization_plots(
                final_binning,
                events,
                mjet_reco_correction=mjet_reco_correction,
                flow=flow,
            )
        return final_binning, optimization_plot

    def optimization_pt_scan(
        self,
        pt_reco_binning: np.ndarray,
        outname_prefix: str = None,
        events: ak.Array = None,
        mjet_reco_correction: Union[Union[np.ndarray, float], Callable] = 1.0,
        bin_merger: Literal["sequential", "bidirectional"] = "sequential",
        split_value: float = 80.0,
        flow: bool = False,
    ):
        final_binnings = []
        for ipt in range(len(pt_reco_binning) - 1):
            pt_low_str, pt_high_str = map(lambda b: str(b).replace(".0", ""), pt_reco_binning[ipt:ipt+2])
            self.logger.info(logging_styles["green_bold"](f"{pt_low_str} < pT < {pt_high_str}"))
            pt_low, pt_high = pt_reco_binning[ipt:ipt+2]
            events_pt_bin = events_sel[(events_sel.pt >= pt_low) & (events_sel.pt < pt_high)]
            mjet_reco_correction_pt_bin = mjet_reco_correction
            if isinstance(mjet_reco_correction_pt_bin, np.ndarray):
                mjet_reco_correction_pt_bin = mjet_reco_correction[ipt]

            pt_bin_optimized_bins, pt_bin_optimization_plot = mjet_bin_optimizer.optimize_binning(
                events_pt_bin,
                bin_merger=bin_merger,
                split_value=split_value,
                plots=True,
                mjet_reco_correction=mjet_reco_correction_pt_bin,
            )
            pt_bin_optimization_plot[1][1].set_title(
                r"$%s \leq p_T^\mathrm{reco} < %s~$GeV" % (pt_low_str, pt_high_str), fontsize=24
            )
            final_binnings.append(pt_bin_optimized_bins)
            if outname_prefix:
                self.save_plot(
                    pt_bin_optimization_plot, f"{outname_prefix}_optization_ptbin_{pt_low_str}_{pt_high_str}.pdf"
                )
        return final_binnings

    def make_optimization_plots(
        self,
        binning: np.ndarray,
        events: ak.Array = None,
        mjet_reco_correction: Union[float, Callable] = 1.0,
        flow=False,
        grid: bool = True,
    ):
        if events is None:
            events = self.events

        # PLOTS
        # build migration matrix with optimized binning
        unfolding_plotter = Unfolding1DPlotter(self.variable, mjet_reco_correction)
        migmat = unfolding_plotter.build_migration_matrix(
            hist.axis.Variable(binning, name=f"{self.variable}_reco", overflow=flow),
            hist.axis.Variable(binning, name=f"{self.variable}_gen", overflow=flow),
            events,
        )
        # build migration matrix with initial binning
        coarse_initial_binning = self.initial_binning[::4]
        migmat_initial = unfolding_plotter.build_migration_matrix(
            hist.axis.Variable(coarse_initial_binning, name=f"{self.variable}_reco", overflow=flow),
            hist.axis.Variable(coarse_initial_binning, name=f"{self.variable}_gen", overflow=flow),
            events,
        )

        # setup figure, grid and first subplots
        if grid:
            f = plt.figure(figsize=(18, 18))
            grid = f.add_gridspec(3, 2, hspace=0.3, height_ratios=[1, 1, 0.05])
            bottom_grid = f.add_gridspec(3, 2, hspace=0.4, height_ratios=[1, 1, 0.05])
            ax01 = f.add_subplot(grid[0])
            ax02 = f.add_subplot(grid[1])
        else:
            (f1, ax1), (f2, ax2), (f01, ax01), (f02, ax02) = tuple(unfolding_plotting.fax() for i in range(4))
            f = [f1, f2, f01, f02]

        # plot stability and purity
        hep.histplot(
            migration_metric(migmat, f"{self.variable}_reco", flow=flow),
            label="stability",
            ax=ax01,
            **{"color": unfolding_plotting.stability_color},
        )
        hep.histplot(
            migration_metric(migmat, f"{self.variable}_gen", flow=flow),
            label="purity",
            ax=ax01,
            **{"color": unfolding_plotting.purity_color},
        )
        ax01.legend(loc="lower center")
        ax01.plot(ax01.get_xlim(), [self.threshold, self.threshold], "k--", alpha=0.6)

        # plot mjetgen with initial and optimized binning and mjetreco with initial binning
        hep.histplot(
            migmat_initial[::sum, :],
            ax=ax02,
            label=r"$%s^{\mathrm{gen}}$" % unfolding_plotting.label_tex_dict[self.variable],
            # density=True,
            binwnorm=True,
        )
        hep.histplot(
            migmat_initial[:, ::sum],
            ax=ax02,
            label=r"$%s^{\mathrm{reco}}$" % unfolding_plotting.label_tex_dict[self.variable],
            # density=True,
            binwnorm=True,
        )

        hep.histplot(
            migmat[::sum, :],
            ax=ax02,
            label=r"$%s{^{\mathrm{gen}}$ (optimized)" % unfolding_plotting.label_tex_dict[self.variable],
            # density=True,
            binwnorm=True,
        )
        if isinstance(mjet_reco_correction, float) and mjet_reco_correction != 1.0:
            ax02.text(
                1 / 3 * ax02.get_xlim()[1],
                ax02.get_ylim()[1] / 2,
                r"$m_{SD}$ corr. = %.2f" % mjet_reco_correction,
                fontsize=18,
            )
        if isinstance(mjet_reco_correction, Callable):
            ax02.text(
                1 / 3 * ax02.get_xlim()[1],
                ax02.get_ylim()[1] / 2,
                r"polynominal $m_{SD}$ corr $f(p_T)$",
                fontsize=18,
            )

        ax02.legend()
        ax02.set_xlabel(r"$%s~$[GeV]" % unfolding_plotting.label_tex_dict[self.variable])
        ax02.set_ylabel(r"Events / GeV")
        unfolding_plotting.cms_label(ax01)

        # MIGMAT renormed plot
        ax1 = f.add_subplot(bottom_grid[2])
        ax2 = f.add_subplot(bottom_grid[3])
        cax1 = f.add_subplot(bottom_grid[4])
        cax2 = f.add_subplot(bottom_grid[5])
        migmat_mjet_arr, bins_x, bins_y = migmat_initial.to_numpy()

        migmat_mjet_genax_normed = unfolding_plotting.safe_div(
            migmat_mjet_arr.T, np.sum(migmat_mjet_arr.T, axis=1)[:, None]
        ).T
        migmat_mjet_recoax_normed = unfolding_plotting.safe_div(
            migmat_mjet_arr, np.sum(migmat_mjet_arr, axis=1)[:, None]
        )

        sns.heatmap(migmat_mjet_genax_normed.T, ax=ax1, norm=LogNorm(), cbar=False)
        f.colorbar(ax1.get_children()[0], cax=cax1, orientation="horizontal")
        ax1.set_title("renormed along gen-axis", fontsize=20)
        sns.heatmap(migmat_mjet_recoax_normed.T, ax=ax2, norm=LogNorm(), cbar=False)
        ax2.set_title("renormed along reco-axis", fontsize=20)
        f.colorbar(ax2.get_children()[0], cax=cax2, orientation="horizontal")

        every_nth_tick = 20
        tick_positions = np.arange(0, len(bins_x), every_nth_tick)
        tick_labels = bins_x[::every_nth_tick]
        for ax in [ax1, ax2]:
            ax.invert_yaxis()
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            ax.set_yticks(tick_positions)
            ax.set_yticklabels(tick_labels)
            ax.set_xlabel(r"$%s^{\mathrm{reco}}~$[GeV]" % unfolding_plotting.label_tex_dict[self.variable])
            ax.set_ylabel(r"$%s^{\mathrm{gen}}~$[GeV]" % unfolding_plotting.label_tex_dict[self.variable])
        return f, (ax01, ax02, ax1, ax2)


if __name__ == "__main__":
    import awkward as ak
    import numpy as np
    import argparse
    import logging

    parser = argparse.ArgumentParser()

    optimization_step_choices = ["pt", "mjet", "mjetSimpleCorr", "mjetCorr", "all"]

    parser.add_argument("--input", "-i", type=str, default="WJetsToQQ_tinyTree.parquet")
    parser.add_argument("--optimize", nargs="+", choices=optimization_step_choices + ["menu"], default=["menu"])
    parser.add_argument("--threshold", "-t", type=float, default=0.5)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logFormatter = logging.Formatter("%(asctime)s %(funcName)s %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("unfolding_binning.log")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    if args.verbose:
        consoleHandler.setLevel(logging.INFO)
    else:
        consoleHandler.setLevel(logging.CRITICAL)

    logger.addHandler(consoleHandler)

    if "menu" in args.optimize:
        from simple_term_menu import TerminalMenu

        optimization_step_menu = TerminalMenu(
            optimization_step_choices,
            multi_select=True,
            show_multi_select_hint=True,
            cursor_index=optimization_step_choices.index("all"),
        )
        optimization_step_indx = optimization_step_menu.show()
        args.optimize = list(optimization_step_menu.chosen_menu_entries)

    if "all" in args.optimize:
        args.optimize = optimization_step_choices

    events_sel = ak.from_parquet(args.input)

    events_sel["pt_raw"] = events_sel.Jets.pt[:, 0]
    events_sel["pt"] = events_sel.pt_raw * events_sel.jecfactor[:, 0]
    events_sel["ptgen"] = events_sel.pt_gen_ak8
    events_sel["mjet_raw"] = events_sel.mjet
    events_sel["mjet"] = events_sel.mjet_raw * events_sel.jecfactor[:, 0]
    events_sel["mjetgen"] = events_sel.msd_gen_ak8
    events_sel["rho"] = 2 * np.log(events_sel.mjet / events_sel.pt)
    events_sel["rhogen"] = 2 * np.log(events_sel.mjetgen / events_sel.ptgen)

    threshold_str = str(args.threshold).replace(".", "p")

    plot_dir = "unfolding_binning_plots"

    def save_plot(fax, outname):
        fax[0].savefig((plot_dir + "/" + outname), bbox_inches="tight")

    pt_gen_binning = np.array([0, 550, 650, 800, 1200])  # , np.inf])
    pt_reco_binning = np.array([500, 550, 650, 725, 800, 1000, 1200, np.inf])

    # pt binning
    if "pt" in args.optimize:
        logger.info(logging_styles["red_bold_underline"]("optimize pt binning (simple)"))
        initial_pt_binning = np.arange(500, 1500, 5)

        unfolding_plotter_pt = Unfolding1DPlotter("pt")
        save_plot(
            unfolding_plotter_pt.plot_migration_metric(pt_gen_binning, events_sel, w=9, h=9),
            "pt_metrics.pdf",
        )
        save_plot(
            unfolding_plotter_pt.plot_migration_matrix(pt_reco_binning[:-1], pt_gen_binning[2:], events_sel, w=9, h=9),
            "pt_migration_matrix.pdf",
        )
        save_plot(
            unfolding_plotter_pt.plot_distributions(pt_reco_binning[:-1], pt_gen_binning[2:], events_sel, w=9, h=9),
            "pt_distributions.pdf",
        )

        pt_bin_optimizer = BinOptimizer(events_sel, initial_pt_binning, args.threshold, variable="pt")

        optimized_pt_binning = pt_bin_optimizer.optimize_binning(plots=True)[0]
        logger.info(f"optimized pt binning: {optimized_pt_binning}")
        logger.info(f"chosen pt binning: {pt_gen_binning}")

    if any("mjet" in step for step in args.optimize):
        # mjet binning
        initial_mjet_binning = np.arange(30, 300, 0.5)
        mjet_bin_optimizer = BinOptimizer(
            events_sel, initial_mjet_binning, args.threshold, variable="mjet", plot_dir=plot_dir
        )

    if "mjet" in args.optimize:
        logger.info(logging_styles["red_bold_underline"]("optimize mjet binning (no correction)"))
        unfolding_plotter_mjet = Unfolding1DPlotter("mjet")
        save_plot(
            unfolding_plotter_mjet.plot_migration_metric(initial_mjet_binning, events_sel, w=9, h=9),
            "mjet_metrics.pdf",
        )
        save_plot(
            unfolding_plotter_mjet.plot_migration_matrix(
                initial_mjet_binning, initial_mjet_binning, events_sel, w=9, h=9
            ),
            "mjet_migration_matrix.pdf",
        )
        save_plot(
            unfolding_plotter_mjet.plot_distributions(initial_mjet_binning, initial_mjet_binning, events_sel, w=9, h=9),
            "mjet_distributions.pdf",
        )

        save_plot(
            unfolding_plotter_mjet.plot_migration_metric_nbins_comparison(30, 300, events_sel, [100, 50, 25]),
            "mjet_metrics_nbins_comparison.pdf"
        )

        optimized_mjet_binning, _ = mjet_bin_optimizer.optimize_binning(bin_merger="bidirectional", split_value=85.0)
        logger.info(f"final binning: {optimized_mjet_binning}")
        save_plot(
            unfolding_plotter_mjet.plot_migration_metric(optimized_mjet_binning, events_sel, w=9, h=9),
            "mjet_metrics_optimized.pdf",
        )
        save_plot(
            unfolding_plotter_mjet.plot_migration_matrix(
                reco_binning(optimized_mjet_binning),
                optimized_mjet_binning,
                events_sel,
                w=9,
                h=9,
            ),
            "mjet_migration_matrix_optimized.pdf",
        )

        logger.info(logging_styles["orange_bold"]("pt-scan bidirectional"))
        no_correction_bidirectional_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning, "mjet_binning", events_sel, bin_merger="bidirectional", split_value=85.0
        )
        logger.info(str(no_correction_bidirectional_binning))

    if "mjetSimpleCorr" in args.optimize:
        logger.info(logging_styles["red_bold_underline"]("optimize mjet binning (simple correction)"))
        simple_pt_dependent_msd_corr = np.array(
            [0.89010989, 0.91011236, 0.93103448, 0.95294118, 0.95294118, 0.97590361, 0.95294118]
        )

        # sequential bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan sequential"))
        simple_correction_sequential_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_simple_msdcorr_sequential",
            events_sel,
            mjet_reco_correction=simple_pt_dependent_msd_corr,
        )
        logger.info(str(simple_correction_sequential_binning))

        # bidirectional bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan bidirectional"))
        simple_correction_bidirectional_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_simple_msdcorr_bidirectional",
            events_sel,
            bin_merger="bidirectional",
            split_value=-1.,
            mjet_reco_correction=simple_pt_dependent_msd_corr,
        )
        logger.info(str(simple_correction_bidirectional_binning))

    if "mjetCorr" in args.optimize:
        logger.info(logging_styles["red_bold_underline"]("optimize mjet binning (\"advanced\" correction (response))"))

        polynomial_msd_correction_set = correctionlib.CorrectionSet.from_file(
            "jms_corrections_quadratic_47e3e54d1c.json"
        )

        # optimization with dedicated JMS from MC (response)
        def polynomial_msd_corr_response(pt):
            return 1.0 / polynomial_msd_correction_set["response_g_jec"].evaluate(pt)

        # sequential bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan sequential"))
        correction_response_sequential_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_polynomial_response_msdcorr_sequential",
            events_sel,
            mjet_reco_correction=polynomial_msd_corr_response,
        )
        logger.info(str(correction_response_sequential_binning))

        # bidirectional bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan bidirectional"))
        correction_response_bidirectional_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_polynomial_response_msdcorr_bidirectional",
            events_sel,
            bin_merger="bidirectional",
            split_value=-1.,
            mjet_reco_correction=polynomial_msd_corr_response,
        )
        logger.info(str(correction_response_bidirectional_binning))

        logger.info(logging_styles["red_bold_underline"]("optimize mjet binning (\"advanced\" correction (means))"))

        # optimization with dedicated JMS from MC (response)
        def polynomial_msd_corr_means(pt):
            return 1.0 / polynomial_msd_correction_set["means_g_jec"].evaluate(pt)

        # sequential bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan sequential"))
        correction_means_sequential_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_polynomial_means_msdcorr_sequential",
            events_sel,
            mjet_reco_correction=polynomial_msd_corr_means,
        )
        logger.info(str(correction_means_sequential_binning))

        # bidirectional bin merging
        logger.info(logging_styles["orange_bold"]("pt-scan bidirectional"))
        correction_means_bidirectional_binning = mjet_bin_optimizer.optimization_pt_scan(
            pt_reco_binning,
            "mjet_binning_polynomial_means_msdcorr_bidirectional",
            events_sel,
            bin_merger="bidirectional",
            split_value=-1.,
            mjet_reco_correction=polynomial_msd_corr_means,
        )
        logger.info(str(correction_means_bidirectional_binning))
