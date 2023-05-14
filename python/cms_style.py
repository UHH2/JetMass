import ROOT

ROOT.gROOT.SetBatch(True)

color_list = [
    857,
    629,
    797,
    414,
    434,
    885,
    905,
    805,
    394,
    835,
    862,
    595,
    880,
    615,
    627,
    803,
    402,
    927,
]

colors = {
    "QCD": 867,
    "WUnmatched": 419,
    "WMatched": 413,
    "ZUnmatched": 797,
    "ZMatched": 800,
    "SingleTop": 800,
    "TTbar": 810,
    "WJets": 413,
    "other": 867,
}


obs_draw_option = "PE1"
draw_option = "H"

ratio_plot = True
additional_pad = False
label_size_modifier = 1.0
xLabelSize = 18.0 * label_size_modifier
yLabelSize = 18.0 * label_size_modifier
xTitleSize = 22.0 * label_size_modifier
yTitleSize = 22.0 * label_size_modifier
xTitleOffset_ratio = 2.8
xTitleOffset = 1.0
yTitleOffset = 1.7

bottom_right_margin_modifier = 1.0


logX = False
logY = False

yplot = 0.7
yratio = 0.3
ymax = 1.0
xmax = 1.0
xmin = 0.0


cms_text = "CMS"
extra_text = "Preliminary Simulation"
extra_text_rel_X = 0.12
font_size_modifier = 1.0
text_padding = 0.1

additional_text_size = 0.3 * font_size_modifier
additional_text_ypos = 1.0
extra_right_margin = 0.0


def draw_lumi(
    plotpad,
    year="2017",
    lumi=41.8,
    do_extra_text=True,
    out_of_frame=True,
    do_cms_text=False,
    private_work=False,
    parent_pad=None,
):
    global additional_text_size
    global additional_text_ypos
    global cms_text
    lumi_text = "(%s) %.2f fb^{-1} (13 TeV)" % (str(year), float(lumi))
    if private_work:
        cms_text = "Private work"
        # extra_text = "________________(CMS_data/simulation)"
        extra_text = "________________(CMS_simulation)"
    latex = ROOT.TLatex()
    latex.SetNDC()

    top_margin = plotpad.GetTopMargin()
    right_margin = plotpad.GetRightMargin()
    left_margin = plotpad.GetLeftMargin()

    # text_padding = 0.4

    lumi_text_size = 0.6 * font_size_modifier
    additional_text_size = lumi_text_size * top_margin

    cms_text_size = (0.6 if private_work else 0.75) * font_size_modifier

    extra_text_size = cms_text_size * 0.76

    y_pos = 1 - top_margin + text_padding * top_margin
    additional_text_ypos = y_pos - top_margin - top_margin * text_padding
    # Lumi Text
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumi_text_size * top_margin)
    if lumi > 0.0:
        latex.DrawLatex(1 - right_margin - extra_right_margin, y_pos, lumi_text)

    if not private_work and not do_cms_text:
        return
    # CMS Text
    if not private_work:
        extra_text_size = cms_text_size
        latex.SetTextFont(61)
    latex.SetTextAlign(11)
    latex.SetTextSize(cms_text_size * top_margin)

    # if(not private_work):
    if out_of_frame:
        latex.DrawLatex(left_margin, y_pos, cms_text)
    else:
        y_pos -= top_margin + top_margin * text_padding
        latex.DrawLatex(left_margin * (1 + text_padding), y_pos, cms_text)

    if do_extra_text:
        # latex.SetTextFont(42 if private_work else 52)
        if not private_work:
            latex.SetTextFont(52)
        extra_x_pos = (
            left_margin * (1 + text_padding)
            if private_work
            else left_margin + extra_text_rel_X * (1 - left_margin - right_margin)
            if out_of_frame
            else left_margin * (1 + text_padding)
        )
        counter = 1
        if not out_of_frame:
            y_pos -= extra_text_size * top_margin
        for extra_t in extra_text.split(" "):
            extra_t = extra_t.replace("_", " ")
            latex.SetTextSize(extra_text_size * top_margin * (1 - 0.1 * (counter - 1)))
            latex.SetText(extra_x_pos, y_pos, extra_t)
            latex.DrawLatex(extra_x_pos, y_pos, extra_t)
            if out_of_frame:
                extra_x_pos += latex.GetXsize() * 1.1
            else:
                y_pos -= extra_text_size * top_margin
            counter += 1


def cms_style():
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetMarkerSize(2.5)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0, "y")

    cmsStyle = ROOT.TStyle("CMS", "CMS approved plots style")

    cmsStyle.SetLegendBorderSize(0)
    cmsStyle.SetFrameBorderMode(0)
    cmsStyle.SetCanvasBorderMode(0)
    cmsStyle.SetPadBorderMode(0)
    cmsStyle.SetPadColor(0)
    cmsStyle.SetCanvasColor(0)
    cmsStyle.SetTitleColor(1)
    cmsStyle.SetStatColor(0)
    cmsStyle.SetFrameFillColor(0)
    cmsStyle.SetEndErrorSize(2)

    cmsStyle.SetPaperSize(20, 26)
    cmsStyle.SetPadTopMargin(0.055)
    cmsStyle.SetPadRightMargin(0.055)

    cmsStyle.SetPadBottomMargin(0.15 * label_size_modifier)
    cmsStyle.SetPadLeftMargin(0.15 * label_size_modifier)

    # cmsStyle.SetTextFont(132)
    cmsStyle.SetTextFont(42)
    cmsStyle.SetTextSize(0.08)
    cmsStyle.SetLabelFont(42, "x")
    cmsStyle.SetLabelFont(42, "y")
    cmsStyle.SetTitleOffset(1.15, "x")
    cmsStyle.SetTitleOffset(1.15, "y")
    cmsStyle.SetTitleOffset(1.15, "z")
    cmsStyle.SetLabelFont(42, "z")
    cmsStyle.SetLabelSize(0.04, "x")
    cmsStyle.SetTitleSize(0.05, "x")
    cmsStyle.SetNdivisions(506, "x")
    cmsStyle.SetLabelSize(0.04, "y")
    cmsStyle.SetTitleSize(0.05, "y")
    cmsStyle.SetNdivisions(506, "y")
    cmsStyle.SetLabelSize(0.05, "z")
    cmsStyle.SetTitleSize(0.06, "z")
    cmsStyle.SetNdivisions(506, "z")

    cmsStyle.SetMarkerStyle(8)
    cmsStyle.SetHistLineWidth(2)
    cmsStyle.SetLineStyleString(2, "[12 12]")

    cmsStyle.SetOptTitle(0)
    cmsStyle.SetOptStat(0)
    cmsStyle.SetOptFit(0)

    cmsStyle.SetPalette(1)
    cmsStyle.SetOptTitle(0)

    ROOT.gROOT.SetStyle("Plain")
    ROOT.gROOT.SetStyle("CMS")


def set_style():
    ROOT.gStyle.SetMarkerSize(2.5)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0, "y")

    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(0)
    ROOT.gStyle.SetCanvasColor(0)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetStatColor(0)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetEndErrorSize(2)

    ROOT.gStyle.SetPaperSize(20, 26)
    ROOT.gStyle.SetPadTopMargin(0.055)
    ROOT.gStyle.SetPadRightMargin(0.055)

    ROOT.gStyle.SetPadBottomMargin(0.15 * label_size_modifier)
    ROOT.gStyle.SetPadLeftMargin(0.15 * label_size_modifier)

    # ROOT.gStyle.SetTextFont(132)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetTextSize(0.08)
    ROOT.gStyle.SetLabelFont(42, "x")
    ROOT.gStyle.SetLabelFont(42, "y")
    ROOT.gStyle.SetTitleOffset(1.15, "x")
    ROOT.gStyle.SetTitleOffset(1.15, "y")
    ROOT.gStyle.SetTitleOffset(1.15, "z")
    ROOT.gStyle.SetLabelFont(42, "z")
    ROOT.gStyle.SetLabelSize(0.04, "x")
    ROOT.gStyle.SetTitleSize(0.05, "x")
    ROOT.gStyle.SetNdivisions(506, "x")
    ROOT.gStyle.SetLabelSize(0.04, "y")
    ROOT.gStyle.SetTitleSize(0.05, "y")
    ROOT.gStyle.SetNdivisions(506, "y")
    ROOT.gStyle.SetLabelSize(0.05, "z")
    ROOT.gStyle.SetTitleSize(0.06, "z")
    ROOT.gStyle.SetNdivisions(506, "z")

    ROOT.gStyle.SetMarkerStyle(8)
    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetLineStyleString(2, "[12 12]")

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetOptTitle(0)


def setup_hist(hist):
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(yTitleSize)
    hist.GetYaxis().SetTitleOffset(yTitleOffset)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(yLabelSize)
    if ratio_plot:
        hist.GetXaxis().SetTitleSize(0.0)
        hist.GetXaxis().SetLabelSize(0.0)
    else:
        # hist.GetXaxis().SetTitle(xTitle)
        hist.GetXaxis().SetTitleFont(43)
        hist.GetXaxis().SetTitleSize(xTitleSize)
        hist.GetXaxis().SetTitleOffset(xTitleOffset)
        hist.GetXaxis().SetLabelFont(43)
        hist.GetXaxis().SetLabelSize(xLabelSize)
        hist.GetXaxis().SetNdivisions(506)
        if hasattr(hist, "GetZaxis"):
            hist.GetZaxis().SetTitleFont(43)
            hist.GetZaxis().SetTitleSize(xTitleSize)
            hist.GetZaxis().SetTitleOffset(1.2 * xTitleOffset)
            hist.GetZaxis().SetLabelFont(43)
            hist.GetZaxis().SetLabelSize(xLabelSize)
            hist.GetZaxis().SetNdivisions(506)

    # if(YRangeUser):
    #     hist.GetYaxis().SetRangeUser(y_range[0],y_range[1])
    # if(XRangeUser):
    #     hist.GetXaxis().SetRangeUser(x_range[0],x_range[1])


def setup_ratio_hist(ratioHist):
    ratioHist.GetYaxis().CenterTitle()
    ratioHist.GetYaxis().SetTitleFont(43)
    ratioHist.GetYaxis().SetTitleSize(yTitleSize)
    ratioHist.GetYaxis().SetTitleOffset(yTitleOffset)
    ratioHist.GetYaxis().SetLabelFont(43)
    ratioHist.GetYaxis().SetLabelSize(yLabelSize)
    ratioHist.GetYaxis().SetNdivisions(506)

    # ratioHist.GetXaxis().SetTitle("m_{SD} [GeV]")
    ratioHist.GetXaxis().SetTitleFont(43)
    ratioHist.GetXaxis().SetTitleSize(xTitleSize)
    ratioHist.GetXaxis().SetTitleOffset(xTitleOffset_ratio)
    ratioHist.GetXaxis().SetLabelFont(43)
    ratioHist.GetXaxis().SetLabelSize(xLabelSize)
    ratioHist.GetXaxis().SetTickLength(0.08)
    ratioHist.GetXaxis().SetNdivisions(506)


def setup_add_hist(addHist):
    addHist.GetYaxis().CenterTitle()
    addHist.GetYaxis().SetTitleFont(43)
    addHist.GetYaxis().SetTitleSize(yTitleSize)
    addHist.GetYaxis().SetTitleOffset(yTitleOffset)
    addHist.GetYaxis().SetLabelFont(43)
    addHist.GetYaxis().SetLabelSize(yLabelSize)
    addHist.GetYaxis().SetNdivisions(503)

    addHist.GetXaxis().SetTitleFont(43)
    addHist.GetXaxis().SetTitleSize(0)
    addHist.GetXaxis().SetTitleOffset(xTitleOffset)
    addHist.GetXaxis().SetLabelFont(43)
    addHist.GetXaxis().SetLabelSize(0)
    addHist.GetXaxis().SetTickLength(0.08)
    addHist.GetXaxis().SetNdivisions(506)


def setup_pads(c, logY=False):
    global yplot, yratio, yadd
    if additional_pad:
        yplot = 0.6
        yratio = 0.25
        yadd = 0.15
        addpad = ROOT.TPad(
            "addpad", "AddPlot", xmin, ymax - yplot - yadd, xmax, ymax - yplot
        )
    else:
        yplot = 0.7
        yratio = 0.3
        yadd = 0.0
        addpad = None
    if ratio_plot:
        plotpad = ROOT.TPad("plotpad", "Plot", xmin, ymax - yplot, xmax, ymax)
        ratiopad = ROOT.TPad(
            "ratiopad",
            "Ratio",
            xmin,
            ymax - yplot - yratio - yadd,
            xmax,
            ymax - yplot - yadd,
        )
    else:
        plotpad = ROOT.TPad("plotpad", "Plot", xmin, ymax - yplot - yratio, xmax, ymax)
        ratiopad = None

    plotpad.SetTopMargin(0.08)
    plotpad.SetLeftMargin(0.14 * label_size_modifier)
    plotpad.SetRightMargin(0.05 * (bottom_right_margin_modifier + 1.2))
    plotpad.SetTicks()
    plotpad.Draw()

    if ratio_plot:
        plotpad.SetBottomMargin(0.016 * label_size_modifier)
        ratiopad.SetTopMargin(0.016)
        ratiopad.SetBottomMargin(0.35 * bottom_right_margin_modifier)
        ratiopad.SetLeftMargin(0.14 * label_size_modifier)
        ratiopad.SetRightMargin(0.05 * (bottom_right_margin_modifier + 1.2))
        ratiopad.SetTicks()
        ratiopad.Draw()
    else:
        plotpad.SetBottomMargin(
            0.1 * label_size_modifier * bottom_right_margin_modifier
        )

    if additional_pad:
        addpad.SetTopMargin(0.016)
        addpad.SetBottomMargin(
            0.016 * label_size_modifier * bottom_right_margin_modifier
        )
        addpad.SetLeftMargin(0.14 * label_size_modifier)
        addpad.SetRightMargin(0.05 * (bottom_right_margin_modifier + 0.2))
        addpad.SetTicks()
        addpad.Draw()

    if logY:
        plotpad.SetLogy()
        c.SetLogy()
    if logX:
        plotpad.SetLogx()
        if ratio_plot:
            ratiopad.SetLogx()
        c.SetLogx()
    return plotpad, ratiopad, addpad
