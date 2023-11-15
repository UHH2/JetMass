#!/usr/bin/env pythonJMS.sh
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from scipy.special import binom
from utils import year_alias
import json

hep.style.use("CMS")


def bernstein(x, lamb, nu):
    return binom(nu, lamb) * (x**lamb) * (1 - x) ** (nu - lamb)


def bernstein_poly(pt, rho, n_pt, n_rho, parameters):
    f = 0
    for i_pt in range(n_pt + 1):
        for i_rho in range(n_rho + 1):
            f += parameters[i_pt][i_rho] * bernstein(rho, i_rho, n_rho) * bernstein(pt, i_pt, n_pt)
    return f


def plot_year_bernstein(fit_dir: str, year: str, data: bool = False):

    config = json.load(open(f"{fit_dir}/WJets{year}Unfoldingconfig.json"))
    fit_results = json.load(open(f"{fit_dir}/FullRunIIfitResult.json"))
    orders = config.get("BernsteinOrders", None)

    pt_edges = np.array([575.0, 650.0, 725.0, 800.0, 1000.0, 1200.0, config.get("pt_cutoff", 1400.0)])
    min_msd, max_msd = (50, 300)
    msd_edges = np.linspace(min_msd, max_msd, 51)
    pt_min = pt_edges[0]
    pt_max = pt_edges[-1]

    rho_min = 2 * np.log(min_msd / pt_max)
    rho_max = -2.1

    pt_pts, msd_pts = np.meshgrid(
        pt_edges[:-1] + 0.3 * np.diff(pt_edges), msd_edges[:-1] + 0.5 * np.diff(msd_edges), indexing="ij"
    )
    rho_pts = 2 * np.log(msd_pts / pt_pts)

    pt_scaled = (pt_pts - pt_min) / (pt_max - pt_min)
    rho_scaled = (rho_pts - rho_min) / (rho_max - rho_min)
    valid_rho_bins = (rho_scaled <= 1) & (rho_scaled >= 0)
    rho_scaled[~valid_rho_bins] = 0

    def param_name(ipt, irho, year):
        return f"tf_params_{year}_pt_par{ipt}_rho_par{irho}"

    tf_param_names = [param_name(ipt, irho, year) for ipt in range(orders[0] + 1) for irho in range(orders[1] + 1)]
    tf_params = {name: fit_results[name] for name in tf_param_names}
    bernstein_params = np.array([tf_params[name][0] for name in tf_params]).reshape(orders[0] + 1, orders[1] + 1)

    def eval_bernstein(pt, rho):
        return bernstein_poly(pt, rho, orders[0], orders[1], bernstein_params)

    eval_bernstein_map = np.vectorize(eval_bernstein)
    bernstein_map = eval_bernstein_map(pt_scaled.flatten(), rho_scaled.flatten()).reshape(pt_scaled.shape)
    bernstein_map[~valid_rho_bins] = np.nan

    f, ax = plt.subplots(figsize=(10, 9))
    for ipt in range(len(pt_edges) - 1):
        msd_edges_invalid_min = msd_edges[:-1][~valid_rho_bins[ipt]]
        msd_edges_invalid_max = msd_edges[1:][~valid_rho_bins[ipt]]
        if len(msd_edges_invalid_min) == 0 or len(msd_edges_invalid_max) == 0:
            continue
        msd_min = np.min(msd_edges_invalid_min)
        msd_max = np.max(msd_edges_invalid_max)
        plt.fill_between(
            [msd_min, msd_max],
            pt_edges[ipt],
            pt_edges[ipt + 1],
            facecolor="none",
            hatch="xxx",
            edgecolor=(0, 0, 0, 0.4),
            linewidth=0.0,
        )

    precision = 100
    z_min = np.floor(precision*np.min(bernstein_map[valid_rho_bins]))/precision + 3/(precision*10)
    z_max = np.ceil(precision*np.max(bernstein_map[valid_rho_bins]))/precision + 3/(precision*10)
    z_offset = np.max([abs(z_max-1.), abs(z_min-1.)])

    cmap = ax.pcolormesh(msd_edges, pt_edges, bernstein_map, cmap="RdGy_r", vmin=1 - z_offset, vmax=1 + z_offset)

    fs = 17
    hep.cms.label("Work in Progress", year=year_alias[year], data=data, fontsize=fs+(0 if data else -3), ax=ax, loc=0)
    cb = f.colorbar(cmap)
    cb.set_label("$R_{p/f}$", loc="center")
    ax.set_xlabel(r"$m_\mathrm{SD}$ [GeV]")
    ax.set_ylabel(r"$p_T$ [GeV]")
    pt_ticklabels = [600, 700, 800, 900, 1000, 1100, 1200, 1300, r"$\infty$"]
    pt_ticklabels_positions = [600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400]
    ax.set_yticks(pt_ticklabels_positions, pt_ticklabels)
    f.savefig(f"{fit_dir}/plots/pretty_bernstein_{year}.pdf")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("fitdir")
    parser.add_argument("--year", default="RunII")
    parser.add_argument("--data", action="store_true")
    args = parser.parse_args()
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"] if args.year == "RunII" else [args.year]
    for year in years:
        plot_year_bernstein(args.fitdir, year, data=args.data)
