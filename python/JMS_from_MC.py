#!/usr/bin/env pythonJMS.sh
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from coffea.nanoevents.methods import vector
import hist
import os
from math import ceil
from scipy.optimize import curve_fit
from scipy.stats import moyal
import logging
from typing import Union, Callable, List
from typing_extensions import Literal
import mplhep as hep
import matplotlib as mpl
from scipy.stats import chisquare
from scipy.stats import kstest
import correctionlib.convert
import correctionlib.schemav2 as cs
from hashlib import sha512
import argparse
import re

mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=[
        "#e66101",
        "#fdb863",
        "#5e3c99",
        "#b2abd2",
        "#a6cee3",
        "#1f78b4",
        "#b2df8a",
        "#33a02c",
    ]
)
hep.style.use("CMS")
logger = logging.getLogger(__name__)


def cms_label(ax, fs=20, year=2017):
    if isinstance(year, str):
        year_substr = re.search("[1678]{2}", year)
        if year_substr:
            year = 2000 + int(year_substr.group())
    # hep.cms.label(label=", Work in Progress", year=year, ax=ax, fontsize=fs)
    hep.cms.label("Preliminary", year=year, ax=ax, fontsize=fs, data=False)


def fax(w=9, h=9):
    return plt.subplots(figsize=(w, h))


def calc_mean(x, y):
    return sum(x * y) / sum(y)


def calc_sigma(x, y):
    mean = calc_mean(x, y)
    return np.sqrt(sum(y * (x - mean) ** 2) / sum(y))


def gauss(x, p0, p1, x0, sigma):
    return p0 + p1 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def landau(x, H, A, x0, sigma):
    return H + A * moyal.pdf(x, x0, sigma)


def landau_gauss(x, p0, p1, x0, sigma):
    # return moyal.pdf(x,x0,sigma)*gauss(x,H,A,x0,sigma)
    return p0 + p1 * np.exp(-((x - x0) ** 2) / (2 * sigma**2)) * moyal.pdf(
        x, x0, sigma
    )


def GOF(obs, fit, ndof):
    chi2, pval = chisquare(obs, fit, ndof)
    kstest_result = kstest(obs, fit)
    print(f"chi2/ndof = {chi2}/{ndof} = {chi2/ndof}")
    print(f"p-val from chi2: {pval}")
    print(kstest_result)
    return kstest_result.pvalue, pval


class JMSFitter(object):
    def __init__(self, x, xerr, ax, poly_dim=2):
        self.x = x
        self.xerr = xerr
        self.ax = ax
        self.jms_artists = []
        self.fit_artists = []
        self.fit_results = []
        self.jms_names = []
        self.poly_dim = poly_dim

    def fit_jms(self, obs):

        jms_corr, fit_diag = correctionlib.convert.ndpolyfit(
            points=[self.x],
            values=obs,
            weights=np.ones_like(self.x),
            varnames=["pt"],
            degree=(self.poly_dim,),
        )
        return (jms_corr, fit_diag)

    def add_jms(self, name, obs, label):
        self.jms_names.append(name)

        ndof = len(obs) - (self.poly_dim + 1)
        self.fit_results.append(self.fit_jms(obs))

        def corr_eval(pt):
            return self.fit_results[-1][0].to_evaluator().evaluate(pt)

        ks_result = GOF(obs, corr_eval(self.x), ndof)
        pt = np.linspace(500, 1200, 300)

        self.jms_artists.append(
            self.ax.errorbar(
                self.x,
                obs,
                xerr=self.xerr,
                linestyle="",
                marker="*",
                label=(
                    label
                    + " | KS p-val: "
                    + str(round(ks_result[0], 2))
                    + r"| $\chi^2$ p-val: "
                    + str(round(ks_result[1], 2))
                ),
            )
        )

        color = self.jms_artists[-1][0].get_color()

        self.fit_artists.append(self.ax.plot(pt, corr_eval(pt), color=color))

    def export(self):
        corrections = []
        for ifit in range(len(self.fit_results)):
            self.fit_results[ifit][0].name = self.jms_names[ifit]
            corrections.append(self.fit_results[ifit][0])
        return corrections


def find_bin(arr, x):
    return np.argmin(abs(arr - x))


def iterative_fit(
    h: hist.Hist,
    iterations: int = 2,
    fit_func: Callable = gauss,
    plot: Union[plt.Axes, bool] = False,
):
    width_multiplier = 1
    edges = h.axes[0].edges
    x = h.axes[0].centers
    y = h.values()
    mean = calc_mean(x, y)
    sigma = calc_sigma(x, y)

    ax = plot  # plot is either True/False or mpl Axes object
    if ax:
        if isinstance(plot, plt.Axes):
            ax = plot
        else:
            f, ax = fax()
        hep.histplot((y, edges), label="init", ax=ax)

    fit_results = []
    for i in range(iterations):
        sigma = abs(sigma)
        xmin = x[0]
        xmax = x[-1]
        if i > 0:
            xmin = mean - sigma * width_multiplier
            xmax = mean + sigma * width_multiplier
            min_bin = find_bin(x, xmin)
            max_bin = find_bin(x, xmax)
            x = x[min_bin:max_bin]
            y = y[min_bin:max_bin]
            edges = edges[min_bin:max_bin+1]
        # print(sigma)
        if len(y) < 4:
            logger.log(
                0,
                "Number of bins got reduced to less than 4."
                + " skipping fit-iteration.",
            )
            continue
        try:
            popt, pcov = curve_fit(
                fit_func, x, y, p0=[min(y), max(y), calc_mean(x, y), calc_sigma(x, y)]
            )
        except RuntimeError as e:
            logger.exception(f"Fit iteration #{i} failed.")
            logger.exception(e)

            fit_results.append(
                {
                    "popt": [np.nan, np.nan, np.nan, np.nan],
                    "pcov": np.nan,
                    "xlim": (xmin, xmax),
                }
            )
            continue

        mean = popt[2]
        sigma = popt[3]
        fit_results.append(
            {
                "popt": popt,
                "pcov": pcov,
                "xlim": (xmin, xmax),
            }
        )
        if plot:
            hep.histplot((y, edges), linestyle="--", label=f"crop {i}", ax=ax)
            x_fine = np.linspace(xmax, xmin, 1000)
            ax.plot(
                x_fine, fit_func(x_fine, *popt), label=f"fit {i} mean={popt[2]:.2f}"
            )
    if plot:
        ax.set_xlim(*fit_results[0]["xlim"])
    return fit_results, ax


class JMSExtractor(object):
    def __init__(self, tree_file, year):
        self.tree = tree_file
        self.year = year
        self.fs = 18
        self.groomings = ["g", "u"]
        self.plot_dir = "."
        self.tex_alias = {
            "jec": "(JEC)",
            "nojec": "(no JEC)",
            "g": "groomed",
            "u": "ungroomed",
        }

    @property
    def tree(self):
        return self._tree

    @tree.setter
    def tree(self, tree_file, triggersf=True):
        # getting raw tree without High-Level behaviour
        raw_tree = ak.from_parquet(tree_file)
        if triggersf:
            raw_tree.weight = raw_tree.weight * raw_tree.triggersf

        raw_tree.Jets = ak.zip(
            {
                "pt": raw_tree.Jets.pt,
                "eta": raw_tree.Jets.eta,
                "phi": raw_tree.Jets.phi,
                "mass": raw_tree.Jets.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )

        self._tree = raw_tree

    @tree.deleter
    def tree(self):
        del self._tree

    @property
    def pt_binning(self):
        return self._pt_binning

    @pt_binning.setter
    def pt_binning(self, edges):
        self._pt_binning = edges

    @pt_binning.deleter
    def pt_binning(self):
        del self._pt_binning

    def pt_centers(self):
        return (self.pt_binning[0:-1] + self.pt_binning[1:]) / 2

    def pt_widths(self):
        return (-self.pt_binning[0:-1] + self.pt_binning[1:]) / 2

    def pt_reco_tex_str(self, i):
        return r"$%.0f < p_{T,\mathrm{reco}} \leq %.0f$ GeV" % (
            self.pt_binning[:-1][i],
            self.pt_binning[1:][i],
        )

    @property
    def plot_dir(self):
        return self._plot_dir

    @plot_dir.setter
    def plot_dir(self, path):
        if not os.path.isdir(path):
            os.makedirs(path)
        self._plot_dir = path

    def construct_hists(self, groomed=True, JEC=True):
        m_min = 30
        m_max = 300
        n_bins = 2 * (m_max - m_min)

        h = hist.Hist(
            hist.axis.Variable(self.pt_binning, name="pt_reco"),
            hist.axis.Regular(n_bins, m_min, m_max, name="msd_reco"),
            hist.axis.Regular(n_bins, m_min, m_max, name="msd_gen"),
            storage=hist.storage.Weight(),
        )

        mreco = self.tree.mjet if groomed else self.tree.Jets.mass[:, 0]
        if JEC:
            mreco = mreco * self.tree.jecfactor[:, 0]

        h.fill(
            pt_reco=self.tree.Jets.pt[:, 0] * self.tree.jecfactor[:, 0],
            msd_reco=mreco,
            msd_gen=self.tree.msd_gen_ak8,
            weight=self.tree.weight,
        )

        h_response = hist.Hist(
            hist.axis.Variable(self.pt_binning, name="pt_reco"),
            hist.axis.Regular(300, 0, 2, name="msd_response"),
            storage=hist.storage.Weight(),
        )

        m_over_30 = (
            (mreco > m_min)
            & (self.tree.msd_gen_ak8 > m_min)
            & (mreco < m_max)
            & (self.tree.msd_gen_ak8 < m_max)
        )

        h_response.fill(
            pt_reco=(self.tree.Jets.pt[:, 0] * self.tree.jecfactor[:, 0])[m_over_30],
            msd_response=(mreco / self.tree.msd_gen_ak8)[m_over_30],
            weight=self.tree.weight[m_over_30],
        )
        return h, h_response

    def create_hist_dict(self):
        self.hists = {
            f"{grooming}_{jec}": self.construct_hists(
                groomed=(grooming == "g"), JEC=(jec == "jec")
            )
            for grooming in ["g", "u"]
            for jec in ["jec", "nojec"]
        }

    def save_control_plots(self):
        groomed_jec_combinations = {
            f"{grooming}_{jec}": f"{self.tex_alias[grooming]} {self.tex_alias[jec]}"
            for grooming in self.groomings
            for jec in ["jec", "nojec"]
        }

        # plot msd responses
        f, ax = fax()
        pt_id = {"pt_reco": sum}
        for hid, hlabel in groomed_jec_combinations.items():
            hep.histplot(self.hists[hid][1][pt_id], ax=ax, label=hlabel, density=True)

        ax.legend(loc="upper left", fontsize=self.fs)
        cms_label(ax, year=self.year)
        ax.set_xlabel(r"$m_{SD,\mathrm{reco}}/m_{SD,\mathrm{gen}}$")
        ax.set_ylabel(r"$\Delta N / N$")
        f.savefig(f"{self.plot_dir}/msd_responses_{self.year}.pdf", bbox_inches="tight")

        # plot msd distributions
        f, ax = fax()
        for hid, hlabel in groomed_jec_combinations.items():
            hep.histplot(
                self.hists[hid][0][{**pt_id, "msd_gen": sum}],
                ax=ax,
                label="reco" + hlabel,
            )
        hep.histplot(
            self.hists["g_jec"][0][{**pt_id, "msd_reco": sum}], ax=ax, label="gen"
        )
        ax.legend(loc="upper left", fontsize=self.fs)
        cms_label(ax, year=self.year)
        ax.set_xlabel(r"$m_{SD}$ [GeV]")
        ax.set_ylabel(r"Events/$\Delta m_{SD}$")
        f.savefig(f"{self.plot_dir}/msd_distributions_{self.year}.pdf", bbox_inches="tight")

    def extract_jms(self):

        jecs = ["jec", "nojec"]
        fit_variables = ["mean_reco", "mean_gen", "response"]
        xlabels = {
            "mean_reco": r"$m_{SD,\mathrm{reco}}$",
            "mean_gen": r"$m_{SD,\mathrm{gen}}$",
            "response": r"$m_{SD,\mathrm{reco}}/m_{SD,\mathrm{gen}}$",
        }
        fitfunc_dict = {"means_reco": gauss, "means_gen": gauss, "response": gauss}

        self.jms_from_mc = {}
        iterations = 2
        for grooming in self.groomings:
            for jec in jecs:
                h_ = self.hists[f"{grooming}_{jec}"]
                n_pt_bins = h_[0].axes["pt_reco"].size
                fit_results = {var: [] for var in fit_variables}
                for fit_variable in fit_variables:
                    f = plt.figure(figsize=(20, 12))
                    columns = 3
                    rows = ceil(n_pt_bins / columns)
                    grid = f.add_gridspec(rows, columns, hspace=0.3)
                    grid_overlap = f.add_gridspec(1, 1)
                    cms_label_ax = f.add_subplot(grid_overlap[0])
                    cms_label_ax.axis("off")
                    cms_label(cms_label_ax, year=self.year)

                    for ipt in range(n_pt_bins):
                        ax_ipt = f.add_subplot(grid[ipt])

                        if fit_variable == "mean_reco":
                            fit_result = iterative_fit(
                                h_[0][{"pt_reco": ipt, "msd_gen": sum}],
                                iterations,
                                fit_func=fitfunc_dict["means_reco"],
                                plot=ax_ipt,
                            )
                        elif fit_variable == "mean_gen":
                            fit_result = iterative_fit(
                                h_[0][{"pt_reco": ipt, "msd_reco": sum}],
                                iterations,
                                fit_func=fitfunc_dict["means_gen"],
                                plot=ax_ipt,
                            )
                        elif fit_variable == "response":
                            fit_result = iterative_fit(
                                h_[1][{"pt_reco": ipt}],
                                iterations,
                                fit_func=fitfunc_dict["response"],
                                plot=ax_ipt,
                            )
                        ax_ipt.text(
                            ax_ipt.get_xlim()[1] * 0.6,
                            ax_ipt.get_ylim()[1] * 0.8,
                            self.pt_reco_tex_str(ipt),
                            fontsize=self.fs - 4,
                        )
                        ax_ipt.set_xlabel(xlabels[fit_variable])
                        ax_ipt.legend(loc="upper left", fontsize=self.fs - 3)

                        fit_results[fit_variable].append(
                            fit_result[0][-1]["popt"][-2]
                            if len(fit_result[0]) > 0
                            else -999.0
                        )

                    f.savefig(
                        f"{self.plot_dir}/{grooming}_{jec}_{fit_variable}_control_plots_{self.year}.pdf",
                        bbox_inches="tight",
                    )
                self.jms_from_mc[f"response_{grooming}_{jec}"] = np.array(
                    fit_results["response"]
                )
                self.jms_from_mc[f"means_{grooming}_{jec}"] = np.array(
                    fit_results["mean_reco"]
                ) / np.array(fit_results["mean_gen"])
        return

    def plot_extracted_jms(self):
        for grooming in self.groomings:
            f, ax = fax(12, 9)

            ax.errorbar(
                self.pt_centers(),
                self.jms_from_mc[f"means_{grooming}_jec"],
                xerr=self.pt_widths(),
                label="dividing means jec",
            )
            ax.errorbar(
                self.pt_centers(),
                self.jms_from_mc[f"means_{grooming}_nojec"],
                xerr=self.pt_widths(),
                label="dividing means nojec",
            )

            ax.errorbar(
                self.pt_centers(),
                self.jms_from_mc[f"response_{grooming}_jec"],
                xerr=self.pt_widths(),
                label="fitting response jec",
            )
            ax.errorbar(
                self.pt_centers(),
                self.jms_from_mc[f"response_{grooming}_nojec"],
                xerr=self.pt_widths(),
                label="fitting response nojec",
            )

            ax.set_ylim(0.8, 1.4)
            ax.legend()
            grooming_text = "groomed" if (grooming == "g") else "ungroomed"
            ax.set_ylabel(f"JMS of {grooming_text} jet mass")
            ax.set_xlabel(r"$p_{T,\mathrm{reco}}$ [GeV]")
            ax.plot(ax.get_xlim(), [1.0, 1.0], "k--", alpha=0.6)
            cms_label(ax, year=self.year)
            f.savefig(
                f"{self.plot_dir}/{grooming}_extracted_jms_{self.year}.pdf", bbox_inches="tight"
            )

    def create_corrections(
            self,
            methods: List[Literal["means", "response"]] = ["means", "response"],
            jecs: List[Literal["jec", "nojec"]] = ["jec", "nojec"],
            poly_dim: int = 2
    ):
        corrections = []
        for grooming in self.groomings:
            grooming_text = "groomed" if grooming == "g" else "ungroomed"
            f, ax = fax(12, 9)
            jms_fits = JMSFitter(self.pt_centers(), self.pt_widths(), ax, poly_dim=poly_dim)
            for method in methods:
                for jec in jecs:
                    jms_name = f"{method}_{grooming}_{jec}"
                    jms_fits.add_jms(
                        f"{jms_name}_{self.year}",
                        self.jms_from_mc[jms_name],
                        (method + " " + grooming_text + " " + self.tex_alias[jec]),
                    )
            ax.legend()
            cms_label(ax, year=self.year)
            ax.set_ylim(0.8, 1.4)
            ax.set_xlabel(r"$p_{T,\mathrm{reco}}$ [GeV]")
            ax.set_ylabel(f"JMS of {grooming_text} jet mass")
            f.savefig(f"{self.plot_dir}/{grooming}_fitted_jms_{self.year}.pdf", bbox_inches="tight")
            corrections += jms_fits.export()
        return corrections


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--tagger", default="n2ddt", choices=["n2ddt", "pNetddt"])
    parser.add_argument("--ndpoly", type=int, default=2, help="degree of polynomial used for fitting the JMS")
    parser.add_argument("--output-suffix", type=str, default="quadratic",
                        help="descriptive suffix used both for jms-correctionset json and control-plot directory.")

    args = parser.parse_args()

    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
    corrections = []

    for year in years:
        wjets_jms_extractor = JMSExtractor(f"WJetsToQQ_tinyTree_{year}_{args.tagger}.parquet", year=year)
        wjets_jms_extractor.plot_dir = f"jms_from_mc_plots_{args.tagger}_{args.output_suffix}"
        wjets_jms_extractor.pt_binning = np.array([500, 550, 650, 725, 800, 1000, 1200])
        print("set pt binning to", wjets_jms_extractor.pt_binning)
        print("constructing and filling hists")
        wjets_jms_extractor.create_hist_dict()
        print("making control plots")
        wjets_jms_extractor.save_control_plots()
        wjets_jms_extractor.extract_jms()
        wjets_jms_extractor.plot_extracted_jms()
        corrections += wjets_jms_extractor.create_corrections(poly_dim=args.ndpoly)

    correction_set = cs.CorrectionSet(
        schema_version=2,
        description=(
            "pT-reco dependent reco-msd corrections derived on groomed (g),"
            "ungroomed (u) mass as well as with (jec) and without (nojec)"
            "applying the JEC on the jet mass."
        ),
        corrections=corrections,
    )

    corrections_set_json = correction_set.json(exclude_unset=True)
    sha512sum = sha512(corrections_set_json.encode("utf-8")).hexdigest()
    fname = f"jms_corrections_{args.tagger}_{args.output_suffix}_{sha512sum[-10:]}.json"
    with open(fname, "w") as fout:
        fout.write(corrections_set_json)
