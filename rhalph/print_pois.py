#!/usr/bin/env python3
from __future__ import print_function
import ROOT
import os
import numpy as np
import prettytable

def extract_pois(fname, poi_pattern, fit_result_name="fit_s"):
    file_ = ROOT.TFile.Open(fname, "READ")
    try:
        fit_result = file_.Get(fit_result_name)
    except BaseException as e:
        print(e)
        print("could not find fit result ({}) in provided file.".format(fit_result))
        exit(-1)

    poi_result_tmp = {}
    pois = []
    for p in fit_result.floatParsFinal():
        if poi_pattern in p.GetName():
            poi_result_tmp[p.GetName()] = [p.getVal(), p.getErrorHi(), p.getErrorLo()]
            pois.append(p.GetName())
    pois.sort()
    return [np.array([poi_result_tmp[name][i] for name in pois]) for i in [0, 1, 2]], pois


def correlation_coefficient(fname, pois):
    file_ = ROOT.TFile.Open(fname, "READ")
    fit_s = file_.Get("fit_s")

    try:
        n_pars = fit_s.floatParsFinal().getSize()

        poi_ind = [fit_s.floatParsFinal().index(p) for p in pois]

        cov = np.array([[fit_s.correlationMatrix()[row][col] for col in range(n_pars)] for row in range(n_pars)])
        cov_inv = np.linalg.inv(cov)
        return np.sqrt(1 - 1 / (cov_inv[poi_ind, poi_ind] * cov[poi_ind, poi_ind]))
    except BaseException:
        return np.array([-1])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fitDiagnostics", help="fitDiagnostics file containing fitResults.", nargs="+", default=[]
    )
    parser.add_argument("-p", "--poipattern", default="r_")

    args = parser.parse_args()

    for fd in sorted(args.fitDiagnostics):
        fit_dir = os.path.dirname(os.path.abspath(fd))
        # print(fit_dir.split("/")[-1])
        print(fd)
        try:
            fr, names = extract_pois(fd, args.poipattern)
            corr = correlation_coefficient(fd, names)
        except:
            print("fit failed.")
            continue
        table = prettytable.PrettyTable()
        table.field_names = names
        table.add_row(fr[0])
        print(table)
        print(corr.sum() / len(corr))
