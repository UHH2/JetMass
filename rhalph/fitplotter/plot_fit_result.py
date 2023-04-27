import ROOT
from ROOT import gROOT, gStyle
import os

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)

gStyle.SetTextFont(43)

gStyle.SetTitleOffset(0.86, "X")
gStyle.SetTitleOffset(1.8, "Y")
# gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadLeftMargin(0.19)
# gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.08)
# gStyle.SetPadRightMargin(0.08)
gStyle.SetPadRightMargin(0.1)
gStyle.SetMarkerSize(0.5)
gStyle.SetHistLineWidth(2)
gStyle.SetTitleSize(0.05, "XYZ")
gStyle.SetLabelSize(0.04, "XYZ")
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)

leftMargin = 0.14


colors = {
    "qcd": 867,
    "QCD": 867,
    "WJets": 413,
    "WJetsUnmatched": 419,
    "WJetsMatched": 413,
    "ZJets": 800,
    "DYJets": 800,
    "ZJetsUnmatched": 797,
    "ZJetsMatched": 800,
    "WUnmatched": 419,
    "WMatched": 413,
    "ZUnmatched": 797,
    "ZMatched": 800,
    "SingleTop": 800,
    "ST_tWch": 41,
    "ST_tch": 43,
    "ST_tWch_top": 41,
    "ST_tWch_antitop": 42,
    "ST_tW_top": 41,
    "ST_tW_antitop": 42,
    "ST_tch_top": 43,
    "ST_tch_antitop": 44,
    "ST_sch": 40,
    "TTbar": 810,
    "TTbar_Hadronic": 810,
    "TTbar_SemiLeptonic": 804,
    "TTToHadronic": 810,
    "TTToSemiLeptonic": 804,
    "TTbar_had": 810,
    "TTbar_semilep": 804,
    "TTbar_dilep": 803,
    "WJets": 413,
    "other": 867,
}

signal_scale = -1

obs_draw_option = "PE1"
fit_draw_option = "H"

xLabelSize = 18.0
yLabelSize = 18.0
xTitleSize = 20.0
yTitleSize = 22.0
xTitleOffset = 2.8
yTitleOffset = 1.5

logX = False
logY = False
stacked_plot = False
ratio_plot = True
yplot = 0.7
yratio = 0.3
ymax = 1.0
xmax = 1.0
xmin = 0.0

y_range = [None, None]
# y_range = [0,2*10**4]
x_range = [None, None]
YRangeUser = all(list(map(lambda a: a is not None, y_range)))
XRangeUser = all(list(map(lambda a: a is not None, x_range)))


def draw_lumi(plotpad, lumi=41.2, extra_text=True, out_of_frame=True, cms_text=False, private_work=False):
    lumi_text = "%.1f fb^{-1} (13 TeV)" % float(lumi)

    latex = ROOT.TLatex()
    latex.SetNDC()

    top_margin = plotpad.GetTopMargin()
    right_margin = plotpad.GetRightMargin()
    left_margin = plotpad.GetLeftMargin()

    text_padding = 0.1

    lumi_text_size = 0.6

    cms_text_size = 0.6 if private_work else 0.75

    extra_text_size = cms_text_size * 0.76
    extra_text_rel_X = 0.12

    y_pos = 1 - top_margin + text_padding * top_margin
    # Lumi Text
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumi_text_size * top_margin)
    latex.DrawLatex(1 - right_margin, y_pos, lumi_text)

    if not private_work and not cms_text:
        return
    # CMS Text
    latex.SetTextFont(42 if private_work else 61)
    latex.SetTextAlign(11)
    latex.SetTextSize(cms_text_size * top_margin)
    cms_text = "private work" if private_work else "CMS"
    if out_of_frame:
        latex.DrawLatex(left_margin, y_pos, cms_text)
    else:
        y_pos -= top_margin + top_margin * text_padding
        latex.DrawLatex(left_margin * (1 + text_padding), y_pos, cms_text)

    if extra_text and not private_work:
        extra_text = "Preliminary Simulation private"
        latex.SetTextFont(52)
        extra_x_pos = (
            left_margin + extra_text_rel_X * (1 - left_margin - right_margin)
            if out_of_frame
            else left_margin * (1 + text_padding)
        )
        counter = 1
        if not out_of_frame:
            y_pos -= extra_text_size * top_margin
        for extra_t in extra_text.split(" "):
            latex.SetTextSize(extra_text_size * top_margin * (1 - 0.1 * (counter - 1)))
            latex.SetText(extra_x_pos, y_pos, extra_t)
            latex.DrawLatex(extra_x_pos, y_pos, extra_t)
            if out_of_frame:
                extra_x_pos += latex.GetXsize() * 1.1
            else:
                y_pos -= extra_text_size * top_margin
            counter += 1


def setup_hist(hist, xTitle=""):
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(yTitleSize)
    hist.GetYaxis().SetTitleOffset(yTitleOffset)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(yLabelSize)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0, "y")
    if ratio_plot:
        hist.GetXaxis().SetTitleSize(0.0)
        hist.GetXaxis().SetLabelSize(0.0)
    else:
        hist.GetXaxis().SetTitle(xTitle)
        hist.GetXaxis().SetTitleFont(43)
        hist.GetXaxis().SetTitleSize(xTitleSize)
        hist.GetXaxis().SetTitleOffset(xTitleOffset)
        hist.GetXaxis().SetLabelFont(43)
        hist.GetXaxis().SetLabelSize(xLabelSize)

    if YRangeUser:
        hist.GetYaxis().SetRangeUser(y_range[0], y_range[1])
    if XRangeUser:
        hist.GetXaxis().SetRangeUser(x_range[0], x_range[1])


def setup_ratio_hist(ratioHist):
    ratioHist.GetYaxis().CenterTitle()
    ratioHist.GetYaxis().SetTitleFont(43)
    ratioHist.GetYaxis().SetTitleSize(yTitleSize)
    ratioHist.GetYaxis().SetTitleOffset(yTitleOffset)
    ratioHist.GetYaxis().SetLabelFont(43)
    ratioHist.GetYaxis().SetLabelSize(yLabelSize)
    ratioHist.GetYaxis().SetNdivisions(506)

    ratioHist.GetXaxis().SetTitle("m_{SD} [GeV]")
    ratioHist.GetXaxis().SetTitleFont(43)
    ratioHist.GetXaxis().SetTitleSize(xTitleSize)
    ratioHist.GetXaxis().SetTitleOffset(xTitleOffset)
    ratioHist.GetXaxis().SetLabelFont(43)
    ratioHist.GetXaxis().SetLabelSize(xLabelSize)
    ratioHist.GetXaxis().SetTickLength(0.08)
    ratioHist.GetXaxis().SetNdivisions(506)


def setup_pads(c):
    if ratio_plot:
        plotpad = ROOT.TPad("plotpad", "Plot", xmin, ymax - yplot, xmax, ymax)
        ratiopad = ROOT.TPad("ratiopad", "Ratio", xmin, ymax - yplot - yratio, xmax, ymax - yplot)
    else:
        plotpad = ROOT.TPad("plotpad", "Plot", xmin, ymax - yplot - yratio, xmax, ymax)
        ratiopad = None

    plotpad.SetTopMargin(0.08)
    plotpad.SetLeftMargin(leftMargin)
    plotpad.SetRightMargin(0.05)
    plotpad.SetTicks()
    plotpad.Draw()

    if ratio_plot:
        plotpad.SetBottomMargin(0.016)
        ratiopad.SetTopMargin(0.016)
        ratiopad.SetBottomMargin(0.35)
        ratiopad.SetLeftMargin(leftMargin)
        ratiopad.SetRightMargin(0.05)
        ratiopad.SetTicks()
        ratiopad.Draw()
    else:
        plotpad.SetBottomMargin(0.1)

    if logY:
        plotpad.SetLogy()
        c.SetLogy()
    if logX:
        plotpad.SetLogx()
        if ratio_plot:
            ratiopad.SetLogx()
        c.SetLogx()
    return plotpad, ratiopad


def plot_fit_result(config={"ModelName": "WMassModel"}, logY=False, fit_shapes_root="fit_shapes.root"):
    f_shapes = ROOT.TFile(config["ModelName"] + "/" + fit_shapes_root, "READ")
    ROOT.TH1.AddDirectory(0)

    out_dir = config["ModelName"] + "/plots/" + fit_shapes_root.replace(".root", "/") + ("/logY/" if logY else "")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for channel_str, channel in config["channels"].items():
        print("plotting", channel_str)
        signals = channel["signal"]
        backgrounds = [bg for bg in channel["samples"] if bg not in signals]
        backgrounds = list(
            map(
                lambda bg: "qcd"
                if ("QCD" in bg and "QcdEstimation" in channel and channel["QcdEstimation"] == "True")
                else bg,
                backgrounds,
            )
        )
        regions = channel["regions"] if "regions" in channel else [""]
        for region in regions:
            for suffix in ["prefit", "postfit"]:
                plot_title = "%s %s %s" % (channel_str, region, suffix)
                c = ROOT.TCanvas(plot_title, plot_title, 600, 600)
                # legend = ROOT.TLegend(0.60,0.65,0.9,0.9)

                legend_bbox = (0.4, 0.7, 0.9, 0.9)

                legend = ROOT.TLegend(*legend_bbox)
                legend.SetFillStyle(0)
                plotpad, ratiopad = setup_pads(c)
                plotpad.cd()

                hist_dir = channel_str + region + "_" + suffix
                h_obs = f_shapes.Get(hist_dir + "/data_obs")
                total_b = h_obs.Clone()
                total_b.Reset()
                total_sb = h_obs.Clone()
                total_sb.Reset()

                shapes_signal = [f_shapes.Get(hist_dir + "/" + signal) for signal in signals]
                shapes_background = [f_shapes.Get(hist_dir + "/" + background) for background in backgrounds]
                for ibackground in range(len(backgrounds)):
                    shapes_background[ibackground].SetLineColor(colors[backgrounds[ibackground]])
                    total_b.Add(shapes_background[ibackground], 1)
                    total_sb.Add(shapes_background[ibackground], 1)
                total_b.SetLineColor(46)
                total_b.SetLineWidth(2)

                for isignal in range(len(signals)):
                    shapes_signal[isignal].SetLineColor(colors[signals[isignal]])
                    total_sb.Add(shapes_signal[isignal])
                    if signal_scale > 0:
                        shapes_signal[isignal].Scale(signal_scale)
                total_sb.SetLineColor(32)

                minima = [hist.GetMinimum() for hist in shapes_signal] + [
                    hist.GetMinimum() for hist in shapes_background
                ]
                while 0.0 in minima:
                    minima.remove(0.0)
                if len(minima) == 0:
                    minima.append(1.0)
                Y_minimum = min(minima)
                Y_maximum = max(h_obs.GetMaximum(), total_b.GetMaximum(), total_sb.GetMaximum())
                h_obs.SetLineColor(1)
                h_obs.SetMarkerStyle(8)
                h_obs.SetTitle(plot_title)
                h_obs.Draw(obs_draw_option)
                h_obs.GetYaxis().SetRangeUser(Y_minimum * 0.9 if logY else 0, (10.0 if logY else 1.1) * Y_maximum)
                setup_hist(h_obs)
                h_obs.Draw(obs_draw_option + "SAME")
                if logY:
                    plotpad.SetLogy()
                legend.AddEntry(h_obs, channel["obs"], "p")

                if stacked_plot:
                    stack = ROOT.THStack()
                    total_err = h_obs.Clone()
                    total_err.Reset()
                    for ibackground in range(len(backgrounds)):
                        total_err.Add(shapes_background[ibackground])
                        shapes_background[ibackground].SetFillColor(colors[backgrounds[ibackground]])
                        stack.Add(shapes_background[ibackground])
                        legend.AddEntry(shapes_background[ibackground], backgrounds[ibackground], "f")

                    for isignal in range(len(signals)):
                        total_err.Add(shapes_signal[isignal])
                        shapes_signal[isignal].SetFillColor(colors[signals[isignal]])
                        stack.Add(shapes_signal[isignal])
                        legend.AddEntry(
                            shapes_signal[isignal],
                            signals[isignal] + ("*%.1f" % signal_scale if signal_scale > 0 else ""),
                            "f",
                        )
                    stack.Draw("HISTSAME")
                    total_err.SetMarkerStyle(0)
                    total_err.SetFillStyle(3204)
                    total_err.SetFillColor(922)
                    total_err.SetLineColor(1)
                    total_err.Draw("E2SAME")
                else:
                    for ibackground in range(len(backgrounds)):
                        shapes_background[ibackground].SetLineWidth(2)
                        shapes_background[ibackground].SetLineStyle(2)
                        shapes_background[ibackground].Draw(fit_draw_option + "SAME")
                        legend.AddEntry(shapes_background[ibackground], backgrounds[ibackground], "l")
                    total_b.Draw(fit_draw_option + "SAME")
                    legend.AddEntry(total_b, "Fit b", "l")

                    for isignal in range(len(signals)):
                        shapes_signal[isignal].SetLineWidth(2)
                        shapes_signal[isignal].Draw(fit_draw_option + "SAME")
                        legend.AddEntry(
                            shapes_signal[isignal],
                            signals[isignal] + ("*%.1f" % signal_scale if signal_scale > 0 else ""),
                            "l",
                        )
                    total_sb.SetLineWidth(2)
                    total_sb.Draw(fit_draw_option + "SAME")
                    legend.AddEntry(total_sb, "Fit s+b", "l")

                h_obs.Draw(obs_draw_option + "SAME")
                legend.SetNColumns(2)
                legend.Draw("SAME")

                normal_ratio = True
                ratiopad.cd()
                if normal_ratio:
                    ratio_hist = h_obs.Clone()
                    ratio_hist.GetYaxis().SetTitle("#frac{obs.}{fit}")
                    ratio_hist.GetYaxis().SetRangeUser(0.5, 1.5)
                    ratio_hist.Divide(total_sb)
                    ratio_hist.SetLineColor(1)
                    ratio_hist.SetLineWidth(2)
                    ratio_hist.SetMarkerColor(1)
                    ratio_hist.SetMarkerStyle(8)
                    setup_ratio_hist(ratio_hist)
                    ratio_hist.GetYaxis().SetRangeUser(0.3, 1.7)

                    ratio_hist_1 = h_obs.Clone()
                    ratio_hist_1.GetYaxis().SetTitle("#frac{obs.}{total b}")
                    ratio_hist_1.Divide(total_b)
                    ratio_hist_1.SetLineColor(46)
                    ratio_hist_1.SetLineWidth(2)
                    ratio_hist_1.SetMarkerStyle(1)

                    ratio_hist.Draw("PE1X0")
                    # ratio_hist_1.Draw('EPSAME')

                    ratio_stat_err = total_b.Clone()
                    ratio_stat_err.Divide(total_b)
                    ratio_stat_err.SetFillStyle(3204)
                    ratio_stat_err.SetFillColor(922)
                    ratio_stat_err.SetLineColor(922)
                    ratio_stat_err.SetMarkerSize(0)

                    ratio_stat_err.SetLineStyle(ROOT.kDashed)

                    ratio_stat_err.Draw("E2SAME")
                    ratio_hist.Draw("PE1X0SAME")

                else:
                    # constructing data-BG/sigma
                    ratio_hist = h_obs.Clone()
                    ratio_hist.GetYaxis().SetTitle("#frac{X - total B}{#sigma_{Data}}")
                    ratio_hist.Add(total_b, -1.0)

                    obs_sigma = ratio_hist.Clone()
                    for i in range(obs_sigma.GetNbinsX() + 1):
                        obs_sigma.SetBinContent(i, ratio_hist.GetBinError(i))
                        obs_sigma.SetBinError(i, 1)
                    ratio_hist.Divide(obs_sigma)

                    # constructing fit-BG/sigma
                    ratio_fit = total_sb.Clone()
                    ratio_fit.Add(total_b, -1)

                    fit_sigma = ratio_fit.Clone()
                    for i in range(fit_sigma.GetNbinsX()):
                        fit_sigma.SetBinContent(i, ratio_fit.GetBinError(i))
                    ratio_fit.Divide(obs_sigma)

                    signal_ratios = []
                    for signal_shape in shapes_signal:
                        signal_sigma = signal_shape.Clone()
                        for i in range(signal_shape.GetNbinsX()):
                            signal_sigma.SetBinContent(i, signal_shape.GetBinError(i))
                            signal_sigma.SetBinError(i, 1)

                        new_signal_ratio = signal_shape.Clone()
                        new_signal_ratio.Divide(obs_sigma)
                        signal_ratios.append(new_signal_ratio)

                    ratio_hist.GetXaxis().SetTitleSize(xTitleSize)
                    setup_ratio_hist(ratio_hist)
                    ratio_hist.Draw("EP")
                    ratio_hist_truncated = ratio_hist.Clone()
                    for i in range(ratio_hist.GetNbinsX() + 1):
                        if ratio_hist.GetBinError(i) > 10:
                            ratio_hist_truncated.SetBinError(i, 1)
                    # ratioYMin = ratio_hist_truncated.GetYaxis().GetYmin()
                    # ratioYMax = ratio_hist_truncated.GetYaxis().GetYmax()
                    ratioYMax = 1.5 * max(
                        abs(ratio_hist_truncated.GetMinimum()), abs(ratio_hist_truncated.GetMaximum())
                    )
                    ratio_hist.GetYaxis().SetRangeUser(-ratioYMax, ratioYMax)
                    ratio_hist.Draw("EPSAME")
                    # ratioHist.GetYaxis().SetRangeUser(-1,1.7)
                    ratio_fit.Draw("HistSAME")
                    for signal_ratio in signal_ratios:
                        signal_ratio.Draw("HistSAME")
                ratioXMin = ratio_hist.GetXaxis().GetXmin()
                ratioXMax = ratio_hist.GetXaxis().GetXmax()

                zeropercent = ROOT.TLine(ratioXMin, 1, ratioXMax, 1)
                plus10percent = ROOT.TLine(ratioXMin, 1.1, ratioXMax, 1.1)
                minus10percent = ROOT.TLine(ratioXMin, 0.9, ratioXMax, 0.9)
                plus10percent.SetLineStyle(ROOT.kDashed)
                minus10percent.SetLineStyle(ROOT.kDashed)
                if normal_ratio:
                    zeropercent.Draw()
                    plus10percent.Draw()
                    minus10percent.Draw()

                plotpad.cd()
                plotpad.RedrawAxis()
                latex = ROOT.TLatex()
                latex.SetNDC(1)
                latex.SetTextFont(42)
                bin_text_size = 2
                if "w" in channel_str.lower():
                    pt_bin = tuple(channel["pt_bin"].split("to"))
                    bin_text = "%s-region  - %s GeV #leq p_{T} < %s GeV"
                    text_length_factor = 1.0 / len(bin_text % ("pass", "XXXX", "XXXX"))
                    latex.SetTextSize(bin_text_size * text_length_factor)
                    latex.DrawLatex(
                        legend_bbox[0],
                        legend_bbox[1] - bin_text_size * text_length_factor,
                        bin_text % (region, pt_bin[0], pt_bin[1]),
                    )
                if "lumi" in config:
                    draw_lumi(plotpad, config["lumi"])
                else:
                    draw_lumi(plotpad)

                c.RedrawAxis()
                c.SaveAs(out_dir + "/" + hist_dir + ("_logY" if logY else "") + ".pdf")
