#!/usr/bin/env python3
import argparse
import numpy as np
from jetmass_scale_fit_utils import correlation_coefficient
import json
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="path to fitDiagnostics.root file containing fit result with cov-matrix")
    parser.add_argument("--verbose", "-v", help="print corr coeffcients to stdout", action="store_true")

    args = parser.parse_args()

    fname = os.path.abspath(args.input)
    model_dir = os.path.dirname(fname)
    config = json.load(open("{}/config.json".format(model_dir)))
    unfolding_bins = config["unfolding_bins"]
    pois = [
        "r_ptgen%i_msdgen%i" % (iptgen, imsdgen)
        for iptgen in range(len(unfolding_bins["ptgen"]) - 1)
        for imsdgen in range(len(unfolding_bins["msdgen"]) - 1)
    ]
    corr_coeff = correlation_coefficient(fname, pois)
    if args.verbose:
        print("average global correlation coeffcient:", np.mean(corr_coeff))
    np.save("{}/corr_coeff.npy".format(model_dir), corr_coeff)
