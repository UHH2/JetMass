#!/usr/bin/env python
from __future__ import print_function
import UHH2.JetMass.plotter as plotter
import ROOT
import os
import numpy as np
from copy import deepcopy
plotter.draw_extra_text = True
plotter.private_work = True


def plot_variable(selection, varname, hist_file_path, output_dir, scaleQCD=True, variable_specifics={}):
    f_hists = ROOT.TFile(hist_file_path, "READ")
    ROOT.TH1.AddDirectory(0)
    logY = variable_specifics.get("logy", False)
    plotter.logY = logY
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    hist_dir = "{}_%s_{}".format(selection, varname)
    try:
        h_data, mc_hists = plotter.get_hists(
            f_hists, plotter.mc_samples[selection], hist_dir, selection, scaleQCD=scaleQCD,
            rebin_factor=variable_specifics.get("rebin", 1)
        )
    except BaseException as e:
        print("Could not get hists ({} {}) from file. Probably they do not exists.".format(selection, varname))
        print(e)
        return
    legend_entries = []
    h_data.GetXaxis().SetTitle(variable_specifics.get("xlabel", "x"))
    if "xrange" in variable_specifics:
        h_data.GetXaxis().SetRangeUser(*variable_specifics["xrange"])
    bkg_stack = ROOT.THStack()
    mc_samples = deepcopy(plotter.mc_samples[selection])
    if selection == "W":
        mc_samples.remove("QCD")
        mc_samples.append("QCD")

    for sample in mc_samples:
        mc_hists[sample].SetFillColorAlpha(plotter.colors.get(sample), 0.8)
        mc_hists[sample].SetLineColor(plotter.colors.get(sample))
        legend_label = plotter.legend_names[sample]
        if selection == "W" and sample == "QCD":
            legend_label += " \\times {:.2f}".format(plotter.qcd_scale)
        legend_entries.append(
            (mc_hists[sample], legend_label, "f")
        )
        bkg_stack.Add(mc_hists[sample], "Hist")
    legend_entries.append((h_data, "Data", "p"))

    ymax = h_data.GetMaximum()
    mc_max_vals = [mc_hists[sample].GetMaximum() for sample in mc_samples]
    ymin_log = min(mc_max_vals)/5.
    print(ymin_log)
    if logY:
        print(ymin_log, 10**(np.ceil(np.log10(ymax))/(2./3.)))
        h_data.GetYaxis().SetRangeUser(ymin_log, 10**(np.log10(ymax)/(2/3.)))
    else:
        h_data.GetYaxis().SetRangeUser(0., ymax/(2./3.))

    plotter.xTitle = variable_specifics["xlabel"]
    plotter.plot_data_mc(
        h_data=h_data,
        h_mc=bkg_stack,
        plot_title=hist_dir.replace("_%s", ""),
        out_dir=output_dir,
        legend_entries=legend_entries,
        # additional_text=additional_text,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--year", default="2017", choices=[
        "2017",
        "UL16preVFP", "UL16postVFP",
        "UL17",
        "UL18"
    ])
    parser.add_argument("--origQCDScale", action="store_true")
    parser.add_argument("--input", "-i", default="templates_2017_1d.root")
    parser.add_argument("--output", "-o", default="")
    args = parser.parse_args()
    plotter.year = args.year
    plotter.rebin = False

    variables = {
        "pt": {"xlabel": r"p_{T} [GeV]", "logy": True, "rebin": 5},
        "eta": {"xlabel": r"\eta", "logy": False},
        "phi": {"xlabel": r"\phi", "logy": False},
        "rho": {"xlabel": r"\rho", "logy": True},
        "npv": {"xlabel": "N_{PV}", "logy": True},
        "n2ddt": {"xlabel": r"N_{2}^{DDT}", "xrange": (-0.2, 0.5), "logy": True},
        "n2": {"xlabel": r"N_{2}", "xrange": (0, 1.0), "logy": True},
        "pNet_TvsQCD": {"xlabel": r"TvsQCD", "xrange": (0, 1.0), "logy": True},
        "pNet_WvsQCD": {"xlabel": r"WvsQCD", "xrange": (0, 1.0), "logy": True},
        "pNet_MD_WvsQCD": {"xlabel": r"MD WvsQCD", "xrange": (0, 1.0), "logy": True},
        "tau21": {"xlabel": r"\frac{\tau_{2}}{\tau_{1}}", "xrange": (0, 1.0), "logy": True},
        "tau32": {"xlabel": r"\frac{\tau_{3}}{\tau_{2}}", "xrange": (0, 1.0), "logy": True},
    }

    for selection in ["W", "top"]:
        for varname, varsettings in variables.items():
            plot_variable(
                selection,
                varname,
                args.input,
                args.output,
                scaleQCD=not args.origQCDScale,
                variable_specifics=varsettings,
            )
