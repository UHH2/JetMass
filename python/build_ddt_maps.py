#!/usr/bin/env pythonJMS.sh
import numpy as np
import hist
from coffea.util import load
from utils import numpy_to_th2
import scipy.ndimage.filters as filters


def build_ddt_map(hqcd: hist.Hist, qcd_eff: float, year: str = "UL17", gauss_radius: int = -1):

    # extract values
    vals_var = hqcd[{"year": year}].values()

    # derive cumulative sum along n2 axis
    cum_sum = np.cumsum(vals_var, axis=2)
    # find max-val and normalise to 1
    # in some (rho,pt)-bins there are no entries/negative default values so cum_sum would sum to negative or zero
    max_val = cum_sum[:, :, -1]
    norma = cum_sum / np.maximum(1e-100, max_val[:, :, np.newaxis])

    # find bins with (qcd-eff)-percentile
    res = np.apply_along_axis(lambda norma: norma.searchsorted(qcd_eff), axis=2, arr=norma)

    # create mask from bins that have been set to default small value
    mask_zero = np.apply_along_axis(lambda norma: norma.searchsorted(1e-100), axis=2, arr=norma)
    # apply mask and set default small values to zero
    res = np.where(mask_zero == norma.shape[2], 0, res)

    # make sure we don't exceed n2-bins
    res = np.where(res > norma.shape[2] - 1, norma.shape[2] - 1, res)

    quantile_map = hqcd.axes["discriminator"].edges[res]

    if gauss_radius <= 0:
        return quantile_map
    smooth_map = filters.gaussian_filter(quantile_map, gauss_radius)
    return smooth_map


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", default="coffea_hists/qcd_pt_v_rho_v_n2_allyears.coffea")
    parser.add_argument("--output", "-o", default="../Histograms/ddtmaps")
    parser.add_argument("--working-points", nargs="+", default=[0.025, 0.05])
    args = parser.parse_args()

    hists = load(args.input)

    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]

    maps = {"numpy": {}, "root": {}}

    for qcd_eff in args.working_points:
        for gauss_radius in [4]:
            smooth = gauss_radius > 0
            for year in years:
                for corrected in ["", "_corrected_pt", "_corrected_pt_mass"]:
                    hqcd = hists[f"leadingJet{corrected}"]
                    print(f"creating map for {qcd_eff} {gauss_radius} {year} {corrected}")
                    discddt = build_ddt_map(hqcd, qcd_eff, year=year, gauss_radius=gauss_radius)

                    pt_edges = hqcd.axes["pt"].edges
                    rho_edges = hqcd.axes["rho"].edges

                    smooth_subtr = "_smooth_" + str(gauss_radius).replace(".", "p") if smooth else ""
                    gauss_radius_str = str(qcd_eff).replace(".", "p")
                    name = f"discddt_map_{year}{smooth_subtr}_{gauss_radius_str}{corrected}"

                    # save as numpy hist
                    maps["numpy"][name] = (discddt.T, rho_edges, pt_edges)
                    # save as TH2 in root
                    maps["root"][name] = numpy_to_th2(discddt, rho_edges, pt_edges, "discddt", "rho", "pt")

    np.save(f"{args.output}.npy", maps["numpy"], allow_pickle=True)
    import uproot

    outfile = uproot.recreate(f"{args.output}.root")
    for name_, map_ in maps["root"].items():
        outfile[name_] = map_
    outfile.close()
