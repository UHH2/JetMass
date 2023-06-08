import ROOT
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom
# from pylab import pcolor, show, colorbar, xticks, yticks
# from pylab import *

sys.path.append("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/")
import cms_style  # noqa: E402


left_margin = 0.14
right_margin = 0.16
top_margin = 0.08
bottom_margin = 0.12


def bernstein(x, lamb, nu):
    return binom(nu, lamb) * (x**lamb) * (1 - x) ** (nu - lamb)


def bernstein_poly(pt, rho, n_pt, n_rho, parameters):
    f = 0
    for i_pt in range(n_pt + 1):
        for i_rho in range(n_rho + 1):
            f += parameters[i_pt][i_rho] * bernstein(rho, i_rho, n_rho) * bernstein(pt, i_pt, n_pt)
    return f


def extract_fit_pars(file_path, order, suffix="params", TFSuffix="tf"):
    """
    file_path: path of root-file containing fit result (i.e. fitDiagnostics.root)
    order: order of BernsteinPolynomials (n_pt,n_rho). Used to derive the number of parameter to be extracted.
    suffix: determines the TF whose params will be extracted.
      'params'-> 1 TF (default),
      'MCtempl' -> 2TF-MCtemplFit,
      'dataResidual' -> 2TF-DataResidualFit
    """
    fit_diagnostics = ROOT.TFile(file_path)
    fitargs = fit_diagnostics.Get("fit_s").floatParsFinal()
    bernstein_pars = []
    bernstein_par_names = []
    for i in range(fitargs.getSize()):
        if (
            "tf" in fitargs.at(i).GetName()
            and suffix in fitargs.at(i).GetName()
            and TFSuffix in fitargs.at(i).GetName()
        ):
            bernstein_pars.append(fitargs.at(i).getVal())
            bernstein_par_names.append(fitargs.at(i).GetName())
    return np.array(bernstein_pars).reshape(order[0] + 1, order[1] + 1), np.array(bernstein_par_names).reshape(
        order[0] + 1, order[1] + 1
    )


def plot_qcd_fail_parameters(config={"ModelName": "WMassModel"}):
    model_dir = config.get("ModelDir", config["ModelName"])
    out_dir = model_dir + "/plots/"
    min_msd, max_msd = (50, 200)
    binwidth = 10
    if "binning" in config:
        w_binning = config["binning"]
        min_msd = w_binning[0]
        max_msd = w_binning[1]
        binwidth = w_binning[2]

    nbins = int(np.floor((max_msd - min_msd) / binwidth))
    msd_bin_edges = np.linspace(min_msd, nbins * binwidth + min_msd, nbins + 1)

    pt_cutoff = config.get("pt_cutoff", 1400)
    pt_bin_edges = np.array(
        [
            (pt_cutoff if pt == "Inf" else pt)
            for pt in config.get("pt_edges", [500.0, 550.0, 600.0, 675.0, 800.0, 1200.0])
        ]
    )

    result = ROOT.TFile(model_dir+ "/fitDiagnostics.root", "READ").Get("fit_s")
    args = result.floatParsFinal()
    qcd_param_names = [
        name
        for name in args.contentsString().split(",")
        if ("qcdparam" in name and config.get("TFSuffix", "qcdparam") in name)
    ]
    # pt_bins = set([int(name.split('_')[2].replace('ptbin','')) for name in qcd_param_names ])
    # msd_bins = set([int(name.split('_')[3].replace('msdbin','')) for name in qcd_param_names ])
    pt_bins = set([int(name.split("ptbin")[1].split("_")[0]) for name in qcd_param_names])
    msd_bins = set([int(name.split("msdbin")[1].split("_")[0]) for name in qcd_param_names])
    th2_qcd_params = ROOT.TH2F("qcdparams", "qcdparam", len(msd_bins), msd_bin_edges, len(pt_bins), pt_bin_edges)
    th2_qcd_params.GetYaxis().SetTitle("p_{T}")
    th2_qcd_params.GetXaxis().SetTitle("m_{SD}")

    for pt in pt_bins:
        for msd in msd_bins:
            param_name = "qcdparam_ptbin%i_msdbin%i_%s" % (pt, msd, config["year"])
            if param_name not in qcd_param_names:
                th2_qcd_params.SetBinContent(msd, pt, 0)
            else:
                param = args.find(param_name)
                val = param.getValV()
                th2_qcd_params.SetBinContent(msd + 1, pt + 1, val)

    cms_style.cms_style()
    c = ROOT.TCanvas("qcdparams", "qcdparams", 1500, 1000)
    cms_style.ratio_plot = False
    plotpad, _, _ = cms_style.setup_pads(c)
    ROOT.gStyle.SetPalette(ROOT.kCool)
    ROOT.gStyle.SetNumberContours(250)
    ROOT.gPad.SetRightMargin(0.18)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    cms_style.setup_hist(th2_qcd_params)
    th2_qcd_params.SetMarkerSize(1)

    max_val = max(th2_qcd_params.GetMaximum(), abs(th2_qcd_params.GetMinimum()))
    th2_qcd_params.GetZaxis().SetRangeUser(-1.1 * max_val, 1.1 * max_val)
    th2_qcd_params.GetZaxis().SetTitle("qcdparam")
    th2_qcd_params.Draw("colztext")
    cms_style.draw_lumi(c, 41.8, do_extra_text=False, out_of_frame=True, do_cms_text=False, private_work=False)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for file_ext in [".pdf", ".png"]:
        c.SaveAs(out_dir + "qcdparams" + file_ext)


def plot_qcd_bernstein3D(config, parameter_suffix="params"):
    model_dir = config.get("ModelDir", config["ModelName"])
    order = (
        tuple(config.get("InitialQCDFitOrders", [2, 2]))
        if parameter_suffix == "MCtempl"
        else tuple(config.get("BernsteinOrders", [2, 2]))
    )
    bernstein_pars, _ = extract_fit_pars(model_dir + "/fitDiagnostics.root", order, parameter_suffix)

    param_list_str = "{"
    for ptbin_params in bernstein_pars:
        param_list_str += "{"
        for param in ptbin_params:
            print(param)
            param_list_str += str(param) + ","
        param_list_str = param_list_str[:-1]
        param_list_str += "},"
    param_list_str = param_list_str[:-1]
    param_list_str += "}"
    final_param_list_str = ""
    for substr in param_list_str.split(","):
        if "e" in substr:
            substr = substr.replace("e", "*10^(") + ")"
        final_param_list_str += substr + ","
    final_param_list_str = final_param_list_str[:-1]
    print("Bernstein: %i,%i,%s" % (order[1], order[0], final_param_list_str))
    out_name = model_dir + "/plots/bernstein_3d_%s.pdf" % parameter_suffix
    os.system("module load mathematica/11.1")
    os.system(
        './fitplotter/bernstein_3d_plot.wls "%s" %i %i "%s"' % (out_name, order[1], order[0], final_param_list_str)
    )


def plot_qcd_bernstein(config={"ModelName": "WMassModel"}, do_3d_plot=True):
    model_dir = config.get("ModelDir", config["ModelName"])
    parameter_suffixes = ["params"]
    if config["InitialQCDFit"] == "True":
        parameter_suffixes = ["MCtempl", "dataResidual"] + parameter_suffixes
    out_dir = model_dir + "/plots/"

    pt_edges = set()
    for channel in config["channels"].keys():
        if config["channels"][channel].get("QcdEstimation", "False") == "True":
            this_bin = config["channels"][channel]["pt_bin"]
            print(this_bin)
            pt_edges.add(float(this_bin.split("to")[0]))
            pt_edges.add(float(this_bin.split("to")[1]))
    pt_edges = list(pt_edges)
    pt_edges.sort()
    pt_cutoff = config.get("pt_cutoff", 1400)
    pt_edges = np.array(
        [
            (pt_cutoff if pt == "Inf" else pt)
            for pt in config.get("pt_edges", [500.0, 550.0, 600.0, 675.0, 800.0, 1200.0])
        ]
    )

    # pt_min = 500.
    # pt_max = 1200.
    pt_min = pt_edges[0]
    pt_max = pt_edges[-1]

    # pt_bins = np.array([500,550,600,675,800,1200],dtype="f")
    pt_bins = np.array(pt_edges, dtype="f")

    min_msd, max_msd = (50, 200)
    binwidth = 10
    if "binning" in config:
        w_binning = config["binning"]
        min_msd = w_binning[0]
        max_msd = w_binning[1]
        binwidth = w_binning[2]

    # rho_min = -6.
    rho_min = 2 * np.log(min_msd / pt_max)
    rho_max = -2.1

    nbins = int(np.floor((max_msd - min_msd) / binwidth))
    msd_bins = np.linspace(min_msd, nbins * binwidth + min_msd, nbins + 1)
    msd_bins = msd_bins.astype("f")
    pt_pts, msd_pts = np.meshgrid(
        pt_bins[:-1] + 0.3 * np.diff(pt_bins), msd_bins[:-1] + 0.5 * np.diff(msd_bins), indexing="ij"
    )
    rho_pts = 2 * np.log(msd_pts / pt_pts)
    pt_scaled = (pt_pts - pt_min) / (pt_max - pt_min)
    rho_scaled = (rho_pts - rho_min) / (rho_max - rho_min)
    valid_rho_bins = (rho_scaled <= 1) & (rho_scaled >= 0)
    rho_scaled[~valid_rho_bins] = 0

    bernstein_maps = {}

    for parameter_suffix in parameter_suffixes:
        order = (
            tuple(config.get("InitialQCDFitOrders", [2, 2]))
            if parameter_suffix == "MCtempl"
            else tuple(config.get("BernsteinOrders", [2, 2]))
        )

        if parameter_suffix == "params" and config.get("InitialQCDFit", "False") == "True":
            bernstein_maps["params"] = bernstein_maps["MCtempl"] * bernstein_maps["dataResidual"]
        else:
            bernstein_pars, _ = extract_fit_pars(
                model_dir + "/fitDiagnostics.root", order, parameter_suffix, config.get("TFSuffix", "tf")
            )

            if parameter_suffix == "MCtempl":
                import sys

                sys.path.append(
                    "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/"
                    "src/UHH2/JetMass/rhalph/rhalphalib"
                )
                import rhalphalib as rl

                rl.util.install_roofit_helpers()
                qcd_fit = (
                    ROOT.TFile(
                        (
                            model_dir
                            + "/qcdfit_"
                            + config["year"]
                            + "_"
                            + config["ModelName"]
                            + config.get("TFSuffix", "")
                            + ".root"
                        )
                    )
                    .Get("w")
                    .genobj("fitresult_qcdmodel_simPdf_qcdmodel_observation")
                )
                mean_values = qcd_fit.valueArray()
                fitargs = qcd_fit.floatParsFinal()
                qcdFit_pars = []
                qcdFit_parnames = []
                for i in range(fitargs.getSize()):
                    if "tf" in fitargs.at(i).GetTitle():
                        qcdFit_pars.append(fitargs.at(i).getVal())
                        qcdFit_parnames.append(fitargs.at(i).GetTitle())
                mean_values = np.array(qcdFit_pars)

                deco_param_values = bernstein_pars.reshape(-1)
                decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult("", qcd_fit, qcdFit_parnames)._transform
                parameter_for_plot = decoVector.dot(np.array(deco_param_values)) + mean_values
                bernstein_pars = parameter_for_plot.reshape(bernstein_pars.shape)

            def eval_bernstein(pt, rho):
                return bernstein_poly(pt, rho, order[0], order[1], bernstein_pars)

            eval_bernstein_map = np.vectorize(eval_bernstein)
            bernstein_map = eval_bernstein_map(pt_scaled.flatten(), rho_scaled.flatten()).reshape(pt_scaled.shape)
            bernstein_map[~valid_rho_bins] = 0
            bernstein_maps[parameter_suffix] = bernstein_map

        bernstein_map = bernstein_maps[parameter_suffix]
        # TH2F
        tf_name = "TF"
        if config.get("InitialQCDFit", "False") == "False":
            tf_name += " (%i,%i)" % order
        else:
            if parameter_suffix == "params":
                mc_order = tuple(config.get("InitialQCDFitOrders", [2, 2]))
                data_order = tuple(config.get("BernsteinOrders", [2, 2]))
                tf_name += "(MCtempl (%i,%i) #times dataResidual (%i,%i))" % (
                    mc_order[0],
                    mc_order[1],
                    data_order[0],
                    data_order[1],
                )
            else:
                tf_name += " %s (%i,%i)" % (parameter_suffix, order[0], order[1])

        th2_name = tf_name
        th2_map = ROOT.TH2F(
            th2_name, th2_name + ";m_{SD} [GeV];p_{T} [GeV]", len(msd_bins) - 1, msd_bins, len(pt_bins) - 1, pt_bins
        )
        th2_map.GetZaxis().SetTitle("TF")
        minRpf = 10
        maxRpf = -1
        for i in range(th2_map.GetNbinsX()):
            for j in range(th2_map.GetNbinsY()):
                map_value = bernstein_map[j][i]
                maxRpf = max(maxRpf, map_value)
                if map_value > 0.0:
                    minRpf = min(minRpf, map_value)
                th2_map.SetBinContent(i + 1, j + 1, map_value)

        th2_map.GetZaxis().SetRangeUser(0.9 * minRpf, 1.1 * maxRpf)
        cms_style.cms_style()
        c = ROOT.TCanvas(th2_name, th2_name, 900, 600)
        cms_style.ratio_plot = False
        plotpad, _, _ = cms_style.setup_pads(c)

        latex = ROOT.TLatex()
        latex.SetNDC()
        top_margin = plotpad.GetTopMargin()
        text_size = 0.6 * top_margin
        y_pos = 1 - top_margin + 0.4 * top_margin

        # Lumi Text
        latex.SetTextFont(42)
        latex.SetTextAlign(22)
        latex.SetTextSize(text_size)
        latex.DrawLatex(0.5, y_pos, th2_name)

        plotpad.cd()

        ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        ROOT.gStyle.SetNumberContours(250)
        ROOT.gPad.SetRightMargin(0.18)
        ROOT.gStyle.SetPaintTextFormat("1.3f")
        cms_style.setup_hist(th2_map)
        th2_map.SetMarkerSize(1)
        th2_map.Draw("textcolz")
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for file_ext in [".pdf", ".png"]:
            c.SaveAs(
                out_dir
                + "TF_th2%s"
                % (("" if parameter_suffix == "params" else "_" + parameter_suffix) + config.get("TFSuffix", ""))
                + file_ext
            )

        # MATPLOTLIB 2D
        plt.switch_backend("agg")
        plt.rc("text", usetex=True)

        figure, ax = plt.subplots()
        ax.set_title("Bernstein map (%i,%i)" % (order[0], order[1]))
        ax.set_xlabel(r"$m_{SD}$ [GeV]")
        ax.set_ylabel(r"$p_{T}$ [GeV]")
        masked_bernstein_map = np.ma.array(bernstein_map, mask=~valid_rho_bins)
        levels = np.linspace(np.min(masked_bernstein_map) - 0.01, np.max(masked_bernstein_map) + 0.01, 500)
        plt_cont = ax.contourf(msd_pts, pt_pts, masked_bernstein_map, levels=levels, cmap="rainbow")

        z_bar = figure.colorbar(plt_cont, format="%.2f")
        z_bar.set_label(r"TF", rotation=270, labelpad=15)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for f_type in [".pdf", ".png"]:
            plt.savefig(out_dir + "bernstein%s" % (parameter_suffix + config.get("TFSuffix", "") + f_type))

        # Mathematica 3D
        if do_3d_plot and (config["InitialQCDFit"] == "True" and parameter_suffix != "params"):
            plot_qcd_bernstein3D(config, parameter_suffix)
