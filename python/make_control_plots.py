#!/usr/bin/env python
from __future__ import print_function
import UHH2.JetMass.plotter as plotter
import ROOT
import os
import numpy as np
from copy import deepcopy
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

plotter.draw_extra_text = True
plotter.private_work = True


def plot_variable(
    selection, varname, hist_file_path, output_dir, scaleQCD=True, variable_specifics={}, matched_samples=False
):
    f_hists = ROOT.TFile(hist_file_path, "READ")
    ROOT.TH1.AddDirectory(0)
    logY = variable_specifics.get("logy", False)
    plotter.logY = logY
    mc_samples = deepcopy(
        plotter.mc_samples[selection] if matched_samples else plotter.mc_samples_nomatching[selection]
    )
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    hist_dir = "{}_%s_{}".format(selection, varname)
    try:
        h_data, mc_hists = plotter.get_hists(
            f_hists, mc_samples, hist_dir, selection, scaleQCD=scaleQCD, rebin_factor=variable_specifics.get("rebin", 1)
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
    if selection == "W":
        mc_samples.remove("QCD")
        mc_samples.append("QCD")

    for sample in mc_samples:
        mc_hists[sample].SetFillColorAlpha(plotter.colors.get(sample), 0.8)
        mc_hists[sample].SetLineColor(plotter.colors.get(sample))
        legend_label = plotter.legend_names[sample]
        if selection == "W" and sample == "QCD":
            legend_label += " \\times {:.2f}".format(plotter.qcd_scale)
        legend_entries.append((mc_hists[sample], legend_label, "f"))
        bkg_stack.Add(mc_hists[sample], "Hist")
    legend_entries.append((h_data, "Data", "p"))

    ymax = max(bkg_stack.GetMaximum(), h_data.GetMaximum())
    mc_max_vals = [mc_hists[sample].GetMaximum() for sample in mc_samples]
    ymin_log = min(mc_max_vals) / 5.0
    if logY:
        print(ymin_log, 10 ** (np.ceil(np.log10(ymax)) / (2.0 / 3.0)))
        h_data.GetYaxis().SetRangeUser(ymin_log, 10 ** (np.log10(ymax) / (2 / 3.0)))
    else:
        h_data.GetYaxis().SetRangeUser(0.0, ymax / (2.0 / 3.0))

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
    parser = argparse.ArgumentParser(
        description=(
            "Create control plots from 1d hists created either from UHH2 code "
            "(make sure to hadd files before using hadd_control_plots.py)"
        )
    )

    parser.add_argument("--year", default="2017", choices=["2017", "UL16preVFP", "UL16postVFP", "UL17", "UL18"])
    parser.add_argument("--origQCDScale", action="store_true")
    parser.add_argument("--input", "-i", default="templates_2017_1d.root")
    parser.add_argument("--output", "-o", default="")
    parser.add_argument("--uhh2", action="store_true", help='use dict of hists from UHH2 framework')
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

    if args.uhh2:
        variables = {
            "npv": {"xlabel": r"N_{PV}", "logy": True},
            "weights": {"xlabel": r"weights", "logy": True},
            "met": {"xlabel": r"E_{T}^{miss}", "logy": True},
            "ht": {"xlabel": r"H_{T}", "logy": True},
            "htlep": {"xlabel": r"H_{T,\text{lep}}", "logy": True},
            "st": {"xlabel": r"S_{T}", "logy": True},
            "nelectron": {"xlabel": r"N_{e}", "logy": True},
            "nmuon": {"xlabel": r"N_{\mu}", "logy": True},
            "pt_muon": {"xlabel": r"p_{T, \mu}", "logy": True},
            "eta_muon": {"xlabel": r"\eta_{\mu}", "logy": False},
            "phi_muon": {"xlabel": r"\varphi_{\mu}", "logy": False},
            "isolation_muon": {"xlabel": r"Iso_{\mu}", "logy": False},
            "ptrel_muon": {"xlabel": r"p_{T,rel}^{\mu}", "logy": False},
            "charge_muon": {"xlabel": r"Q_{\mu}", "logy": False},
            "deltaRmin_muon": {"xlabel": r"\Delta R_{min}(muon, jet)", "logy": False},
            # "deltaRmin_ptrel_muon":{"xlabel":"\Delta R_min", "logy": False},
            "njets": {"xlabel": r"N_{AK4}", "logy": True},
            "pt_jet": {"xlabel": r"p_{T, AK4}", "logy": True},
            "eta_jet": {"xlabel": r"\eta_{AK4}", "logy": False},
            "phi_jet": {"xlabel": r"\varphi_{AK4}", "logy": False},
            "mass_jet": {"xlabel": r"m_{AK4}", "logy": False},
            "ntopjets": {"xlabel": r"N_{AK8}", "logy": True},
            "nsubjets": {"xlabel": r"N_{subjets}", "logy": True},
            "pt_topjet": {"xlabel": r"p_{T, AK8}", "logy": True},
            "eta_topjet": {"xlabel": r"\eta_{AK8}", "logy": False},
            "phi_topjet": {"xlabel": r"\varphi_{AK8}", "logy": False},
            "mass_topjet": {"xlabel": r"m_{AK8}", "logy": False},
            "deltaR_lepton_allJets": {"xlabel": r"\Delta R (leptons, jets)", "logy": True},
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
