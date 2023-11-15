from __future__ import print_function
import ROOT
import numpy as np
import cms_style
import collections
import logging
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

logger = logging.getLogger()

cms_style.set_style()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)


leftMargin = 0.14
ultrawide = False
legend_on_extern_canvas = False
legend_columns = 1
legend_bbox = (0.60, 0.6, 0.9, 0.9)

new_binning_dict = {
    "top": np.linspace(0, 350, 71),
    "W": np.linspace(40, 160, 13),
    "WfromTop": np.linspace(0, 200, 40),
}
# new_binning_dict = {'top':np.linspace(-10,350,73),'W':np.linspace(-10,160,18),'WfromTop':np.linspace(-10,200,42)}
# myBinning
# new_binning_dict = {'top':np.linspace(0,350,71),'W':np.linspace(0,500,101),'WfromTop':np.linspace(0,300,61)}

binning_dict = {
    "ATLAS": {
        "top": np.linspace(0, 350, 71),
        "W": np.linspace(40, 160, 13),
        "WfromTop": np.linspace(0, 200, 40),
    },
    # "CMS": {'top': np.linspace(0, 350, 71), 'W': np.linspace(0, 500, 101), 'WfromTop': np.linspace(0, 300, 61)},
    # "CMS": {'top': np.linspace(50, 200, 16), 'W': np.linspace(0, 500, 101), 'WfromTop': np.linspace(0, 300, 61)},
    # "CMS": {'top': np.linspace(50, 300, 51), 'W': np.linspace(50, 200, 16), 'WfromTop': np.linspace(0, 300, 61)},
    "CMS": {
        # "top": np.linspace(0, 300, 61),
        # "W": np.linspace(50, 300, 51),
        "top": np.linspace(30, 300, 55),
        "W": np.linspace(30, 300, 55),
        "Zbb": np.linspace(50, 300, 51),
    }
    # "CMS": {'top': np.linspace(0, 500, 11), 'W': np.linspace(0, 500, 11), 'Zbb': np.linspace(0, 500, 11)}
    # "CMS": {'top': np.linspace(50, 300, 26), 'W': np.linspace(50, 200, 16), 'WfromTop': np.linspace(0, 300, 61)}
}

rebin = True

colors = {
    "QCD": 867,
    "qcd": 867,
    "QCD_lowPt": 867,
    "QCD_highPt": 867,
    "WJets": 413,
    "WUnmatched": 419,
    "WMatched": 413,
    "WJetsUnmatched": 419,
    "WJetsMatched": 413,
    "WJetsMatched0p4": 413,
    "WJetsMatched0p2": 413,
    "WJetsMatched0p1": 413,
    "WJetsMatched_fakes": 430,
    "ZJets": 797,
    "DYJets": 797,
    "ZJetsUnmatched": 794,
    "ZJetsMatched": 797,
    "ZUnmatched": 794,
    "ZMatched": 797,
    "SingleTop": 800,
    "ST": ROOT.kBlue - 2,
    "ST_tWch": 41,
    "ST_tch": 43,
    "ST_tW_top": 41,
    "ST_tW_antitop": 42,
    "ST_tWch_top": 41,
    "ST_tWch_antitop": 42,
    "ST_tch_top": 43,
    "ST_tch_antitop": 44,
    "ST_sch": 40,
    "ST_tW": 41,
    "ST_t": 43,
    "ST_t_top": 43,
    "ST_t_antitop": 44,
    "ST_s": 40,
    "ttbar": 810,
    "TTbar": 810,
    "TTbar_Hadronic": 810,
    "TTbar_SemiLeptonic": 804,
    "TTToHadronic": 810,
    "TTToSemiLeptonic": 804,
    "TTTo2L2Nu": 803,
    "TTbar_had": 810,
    "TTbar_semilep": 804,
    "TTbar_dilep": 803,
    "TTbar_semilep_mergedTop": 810,
    "TTbar_semilep_mergedW": 633,
    "TTbar_semilep_mergedQB": 634,
    "TTbar_semilep_semiMergedTop": 634,
    "TTbar_semilep_notMerged": 635,
    "TTToSemiLeptonic_mergedTop": ROOT.kRed,
    "TTToSemiLeptonic_mergedW": 798,
    "TTToSemiLeptonic_mergedQB": 792,
    "TTToSemiLeptonic_notMerged": 635,
    "TTbarNonSemiLep": 804,
    "TTbarMergedTop": ROOT.kRed,
    "TTbarMergedW": 798,
    "TTbarMergedQB": 792,
    "TTbarNotMerged": ROOT.kRed - 2,
    "other": 921,
    "other_vjets": 921,
    "other_ttbar": 921,
}
MC_stat_err_color = 922
MC_stat_err_fillstyle = 3235  # 3204

legend_names = {
    "qcd": "QCD",
    "QCD": "QCD",
    "WJets": "W+Jets",
    "WUnmatched": "W+Jets (not merged)",
    "WMatched": "W+Jets (merged W)",
    "WJetsUnmatched": "W+Jets (not merged)",
    "WJetsMatched": "W+Jets (merged W)",
    "WJetsMatched_fakes": "W+Jets (merged W) (fakes)",
    # 'WUnmatched': "W+Jets (Unmatched)",
    # 'WMatched': "W+Jets (Matched)",
    # 'WJetsUnmatched': "W+Jets (Unmatched)",
    # 'WJetsMatched': "W+Jets (Matched)",
    "ZJets": "Z+Jets",
    "DYJets": "DY+Jets",
    "ZJetsUnmatched": "Z+Jets (not merged)",
    "ZJetsMatched": "Z+Jets (merged Z)",
    "ZUnmatched": "Z+Jets (not merged)",
    "ZMatched": "Z+Jets (merged Z)",
    # 'ZJetsUnmatched': "Z+Jets (Unmatched)",
    # 'ZJetsMatched': "Z+Jets (Matched)",
    # 'ZUnmatched': "Z+Jets (Unmatched)",
    # 'ZMatched': "Z+Jets (Matched)",
    "ST": "Single top",
    "SingleTop": "Single top",
    "ST_tWch": "Single anti-/top (tW-ch.)",
    "ST_tch": "Single anti-/top (t-ch.)",
    "ST_tW_top": "Single top (tW-ch.)",
    "ST_tW_antitop": "Single antitop (tW-ch.)",
    "ST_tWch_top": "Single top (tW-ch.)",
    "ST_tWch_antitop": "Single antitop (tW-ch.)",
    "ST_tch_top": "Single top (t-ch.)",
    "ST_tch_antitop": "Single antitop (t-ch.)",
    "ST_sch": "Single anti-/top (s-ch.)",
    "ST_tW": "Single anti-/top (tW-ch.)",
    "ST_t": "Single anti-/top (t-ch.)",
    "ST_s": "Single anti-/top (s-ch.)",
    "ttbar": "t#bar{t}",
    "TTbar": "t#bar{t}",
    "TTbar_Hadronic": "t#bar{t} (hadronic)",
    "TTbar_SemiLeptonic": "t#bar{t} (semi-leptonic)",
    "TTToHadronic": "t#bar{t} (hadronic)",
    "TTToSemiLeptonic": "t#bar{t} (semi-leptonic)",
    "TTbar_had": "t#bar{t} (hadronic)",
    "TTbar_semilep": "t#bar{t} (semi-leptonic)",
    "TTbar_dilep": "t#bar{t} (di-leptonic)",
    "TTbar_semilep_mergedTop": "t#bar{t} (semi-leptonic) (fully merged)",
    "TTbar_semilep_mergedW": "t#bar{t} (semi-leptonic) (merged W)",
    "TTbar_semilep_mergedQB": "t#bar{t} (semi-leptonic) (mergedQB)",
    "TTbar_semilep_semiMergedTop": "t#bar{t} (semi-leptonic) (semi merged)",
    "TTbar_semilep_notMerged": "t#bar{t} (semi-leptonic) (not merged)",
    "TTToSemiLeptonic_mergedTop": "t#bar{t} (fully merged)",
    "TTToSemiLeptonic_mergedW": "t#bar{t} (merged W)",
    "TTToSemiLeptonic_mergedQB": "t#bar{t} (mergedQB)",
    "TTToSemiLeptonic_semiMergedTop": "t#bar{t} (semi merged)",
    "TTToSemiLeptonic_notMerged": "t#bar{t} (not merged)",
    "TTTo2L2Nu": "t#bar{t} (di-leptonic)",
    "TTbarNonSemiLep": "t#bar{t} (hadronic & di-leptonic)",
    "TTbarMergedTop": "t#bar{t} (fully merged)",
    "TTbarMergedW": "t#bar{t} (merged W)",
    "TTbarMergedQB": "t#bar{t} (mergedQB)",
    "TTbarNotMerged": "t#bar{t} (not merged)",
    "other": "other",
    "other_ttbar": "other",
    "other_vjets": "other",
    "others": "other",
}


samples = {
    "W": [
        "QCD",
        "ST_tW_top",
        "ST_tW_antitop",
        "TTToSemiLeptonic",
        "TTToHadronic",
        "ZJetsUnmatched",
        "ZJetsMatched",
        "WJetsUnmatched",
        "WJetsMatched",
    ],
    "Zbb": [
        "QCD",
        "ST_tW_top",
        "ST_tW_antitop",
        "TTToSemiLeptonic",
        "TTToHadronic",
        "ZJetsUnmatched",
        "ZJetsMatched",
        "WJetsUnmatched",
        "WJetsMatched",
    ],
    "WfromTop": [
        "QCD",
        "ST_sch",
        "ST_tch",
        "ST_tWch",
        "DYJets",
        "WJets",
        "TTbar_had",
        "TTbar_dilep",
        "TTbar_semilep",
    ],
    # "top": ["QCD",  "ST_sch",  "ST_tch", "ST_tWch",  "DYJets",  "WJets",  "TTbar_had", "TTbar_dilep",
    # "TTbar_semilep_mergedTop",  "TTbar_semilep_mergedW",  "TTbar_semilep_mergedQB",  "TTbar_semilep_notMerged"],
    "top": [
        "QCD",
        "ST_s",
        "ST_t",
        "ST_tW",
        "DYJets",
        "WJets",
        "TTToHadronic",
        "TTTo2L2Nu",
        "TTToSemiLeptonic_mergedTop",
        "TTToSemiLeptonic_mergedW",
        "TTToSemiLeptonic_mergedQB",
        "TTToSemiLeptonic_semiMergedTop",
        "TTToSemiLeptonic_notMerged",
    ],
}
# "TTbar_semilep_semiMergedTop"

mc_samples = {
    "W": [
        "other_vjets",
        "QCD",
        "ZJets",
        # "ZJetsUnmatched",
        # "ZJetsMatched",
        "WJetsUnmatched",
        "WJetsMatched",
    ],
    "Zbb": [
        "QCD",
        "ST_tW_top",
        "ST_tW_antitop",
        "TTToSemiLeptonic",
        "TTToHadronic",
        "ZJetsUnmatched",
        "ZJetsMatched",
        "WJetsUnmatched",
        "WJetsMatched",
    ],
    "WfromTop": ["other_ttbar", "WJets", "ST", "TTbar"],
    # "top": ["other_ttbar",  "WJets",  "ST",  "TTbar"],
    # "top": ["other_ttbar",  "WJets",  "ST",  "TTbarNonSemiLep", "TTbar_semilep_notMerged",  "TTbar_semilep_mergedQB",
    # "TTbar_semilep_mergedW",  "TTbar_semilep_mergedTop" ],
    # "top": [
    #     "WJets",
    #     "DYJets",
    #     "TTToHadronic",
    #     "TTToSemiLeptonic",
    #     "TTTo2L2Nu",
    #     "ST_t",
    #     "ST_tW",
    #     "ST_s",
    #     "QCD",
    # ],
    "top": [
        "other_ttbar",
        "WJets",
        "ST",
        "TTbarNotMerged",
        "TTToSemiLeptonic_mergedQB",
        "TTToSemiLeptonic_mergedW",
        "TTToSemiLeptonic_mergedTop",
    ],
}
mc_samples_nomatching = {
    "W": [
        "other_vjets",
        "QCD",
        "ZJets",
        "WJets",
    ],
    "top": [
        "other_ttbar",
        "WJets",
        "TTToHadronic",
        "TTToSemiLeptonic",
        "TTTo2L2Nu",
        "ST",
    ],
}

merged_hists = {
    "TTbar": [
        "TTbar_had",
        "TTbar_dilep",
        "TTbar_semilep_mergedTop",
        "TTbar_semilep_mergedW",
        "TTbar_semilep_mergedQB",
        "TTbar_semilep_semiMergedTop",
        "TTbar_semilep_notMerged",
    ],
    # 'TTbarMergedTop': ['TTbar_mergedTop'],
    # 'TTbarMergedW': ['TTbar_mergedW'],
    # 'TTbarMergedQB': ['TTbar_mergedQB'],
    # "TTbarNotMerged": ["TTbar_notMerged"],
    # 'TTbarMergedTop': ['TTbar_semilep_mergedTop', 'TT_Mtt_mergedTop'],
    # 'TTbarMergedW': ['TTbar_semilep_mergedW', 'TT_Mtt_mergedW'],
    # 'TTbarMergedQB': ['TTbar_semilep_mergedQB', 'TT_Mtt_mergedQB'],
    # "TTbarNonSemiLep": ["TTbar_had", "TTbar_dilep"],
    # "TTbarNotMerged": ["TTbar_had", "TTbar_dilep", "TTbar_semilep_notMerged", 'TT_Mtt_notMerged'],
    # 'TTbarMergedTop': ['TTbar_semilep_mergedTop'],
    # 'TTbarMergedW': ['TTbar_semilep_mergedW'],
    # 'TTbarMergedQB': ['TTbar_semilep_mergedQB'],
    # "TTbarNonSemiLep": ["TTbar_had", "TTbar_dilep"],
    # "TTbarNotMerged": ["TTbar_had", "TTbar_dilep", "TTbar_semilep_notMerged"],
    "TTbarMergedTop": ["TTToSemiLeptonic_mergedTop"],
    "TTbarMergedW": ["TTToSemiLeptonic_mergedW"],
    "TTbarMergedQB": ["TTToSemiLeptonic_mergedQB"],
    "TTbarNonSemiLep": ["TTToHadronic", "TTTo2L2Nu"],
    "TTbarNotMerged": ["TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonic_notMerged"],
    # "ST": ["ST_sch",  "ST_tch", "ST_tWch"],
    "ST": ["ST_s", "ST_t", "ST_tW"],
    "other_ttbar": ["QCD", "DYJets"],
    # "other_vjets": ["ST_tW_top", "ST_tW_antitop", "TTToSemiLeptonic", "TTToHadronic"]
    "other_vjets": ["TTToSemiLeptonic", "TTToHadronic"],
}

selection_tex = {
    "top": "t#bar{t}#rightarrow #mu + jets",
    "W": "W(q#bar{q}) + jets",
    "Zbb": "Zbb + jets",
}

regions = {
    "top": ["_pass", "_passW", "_fail"],
    "W": ["_pass", "_fail"],
    "Zbb": ["_pass", "_fail"],
}

region_tex = {
    "top": {
        "": "",
        "pass": "pass-top region",
        "passW": "pass-W region",
        "fail": "fail region",
    },
    "W": {"": "", "pass": "pass region", "fail": "fail region"},
    "Zbb": {"": "", "pass": "pass region", "fail": "fail region"},
}

pt_bins_dict = {
    # "top": ["inclusive", "200to250", "250to300", "200to300", "300to400", "400to500", "500to600", "300to500",
    # "300to350", "350to400", "400to450", "450to500", "500to550", "550to600", "600toInf", "400toInf"],
    # "W":  ["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf"],
    "top": {
        "jms": [
            "inclusive",
            "200to300",
            "300to400",
            "400to500",
            "500to650",
            "650toInf",
        ],
        "unfolding": [
            "inclusive",
            "200to300",
            "300to400",
            "400to500",
            "500to650",
            "650toInf",
        ],
    },
    "W": {
        "jms": [
            "inclusive",
            "500to650",
            "650to800",
            "800to1200",
            "1200toInf",
        ],
        "unfolding": [
            "inclusive",
            "500to575",
            "575to650",
            "650to725",
            "725to800",
            "800to1000",
            "1000to1200",
            "1200toInf",
        ],
    },
    "Zbb": {
        "jms": [
            "inclusive",
            "500to650",
            "650to800",
            "800to1200",
            "1200toInf",
        ],
        "unfolding": [
            "inclusive",
            "500to575",
            "575to650",
            "650to725",
            "725to800",
            "800to1000",
            "1000to1200",
            "1200toInf",
        ],
    }
    # "W": ["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200to1600", "1600to2000",
    # "2000to3000", "3000to4000", "4000to5000", "5000to6000", "6000toInf"],
    # "top": ["inclusive", "200to250", "250to300", "300to400", "400to500", "500toInf", "300to500", "300to350",
    # "350to400", "400to450", "450to500", "500to550", "550to600", "600toInf"],
    # # "top": ["550to600"],
    # "W": ["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf"],
    # "W": ["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf",
    #   "500to600", "600to750", "750to900", "900to1200", "500to1200"],
    # "WfromTop": ["inclusive", "200to300", "300to400", "400to500", "500toInf",
    #       "200to225", "225to250", "250to275", "275to300", "300to325", "325to350", "350to400", "400to500", "200to350"],
}

pt_bins_tex_dict = {
    selection: {workflow: {} for workflow in pt_bins_dict[selection].keys()} for selection in pt_bins_dict.keys()
}
for selection, workflows in pt_bins_dict.items():
    for workflow, bins in workflows.items():
        for pt_bin in bins:
            pt_bins_tex_dict[selection][workflow][pt_bin] = (
                ""
                if (pt_bin == "inclusive")
                else " %s GeV #leq p_{T} < %s GeV" % (pt_bin.split("to")[0], pt_bin.split("to")[1])
            )
        min_pt = min(map(lambda x: int(x.split("to")[0]) if "to" in x else 9999, bins))
        pt_bins_tex_dict[selection][workflow]["inclusive"] = " p_{T} #geq %s GeV" % (min_pt)
nbjet_bins_dict = {"top": [""], "W": [""], "Zbb": ["", "Nbjeteq0", "Nbjetgt0"]}

nbjet_bins_tex_dict = {
    "": "",
    "Nbjeteq0": "N_{b-tag outside AK8} = 0",
    "Nbjetgt0": "N_{b-tag outside AK8} > 0",
}
lumis = {
    "2017": 41528.9954019578 / 1000.0,
    "UL16preVFP": 19301.591954787407 / 1000.0,
    "UL16postVFP": 16626.734093195286 / 1000.0,
    "UL17": 41479.68052876168 / 1000.0,
    "UL18": 59832.47533908866 / 1000.0,
}

year_alias = {
    "UL16preVFP": "legacy 2016 (early)",
    "UL16postVFP": "legacy 2016 (late)",
    "UL17": "legacy 2017",
    "UL18": "legacy 2018",
}

lumis["RunII"] = sum([lumi for year, lumi in lumis.items() if "UL" in year])
lumis["UL16"] = sum([lumi for year, lumi in lumis.items() if "UL16" in year])


obs_draw_option = "PE1X0"
draw_option = "H"
ratio_draw_options = "PE1X0"

obs_line_color = 1
obs_marker_style = 8
obs_marker_size = 0.5
year = "2017"
# extra_text = "Preliminary"
extra_text = "Work_in_progress"
lumi_text_padding = 0.4
additional_text_padding = 0.4
additional_text_size_modifier = 1.4
draw_extra_text = True
private_work = False

ratio_plot = True
ratio_hist_yTitle = "#frac{Data}{Simulation}"
xTitle = "m_{SD} [GeV]"
yTitle = "Events"
xLabelSize = 18.0
yLabelSize = 18.0
xTitleSize = 20.0
yTitleSize = 22.0
xTitleOffset = 2.8
yTitleOffset = 1.5

logX = False
logY = False

yplot = 0.7
yratio = 0.3
ymax = 1.0
xmax = 1.0
xmin = 0.0

qcd_scale = 1.0

y_range = [None, None]
y_range_ratio = [None, None]
x_range = [None, None]
YRangeUser = all(list(map(lambda a: a is not None, y_range)))
XRangeUser = all(list(map(lambda a: a is not None, x_range)))


def get_hists(
    f_hists,
    samples,
    hist_dir="_%s_mjet_inclusive_pass",
    selection="top",
    pseudo_data=False,
    include_merged_hists=True,
    scaleQCD=True,
    pseudo_data_file=None,
    rebin_factor=1,
):
    # print(f_hists, samples, hist_dir,selection)
    hist_dir = selection + hist_dir if (hist_dir[0] == "_") else hist_dir
    print(hist_dir)
    new_binning = None
    if rebin:
        new_binning = binning_dict["CMS"].get(selection, None)

    h_data = None
    if not pseudo_data:
        for data_name in ["Data", "data_obs", "data"]:
            h_data = f_hists.Get(str(hist_dir % data_name))
            # print(hist_dir % data_name)
            try:
                h_data.GetName()
                break
            except BaseException as e:
                logger.info("tried getting data hist with " + data_name + "(which failed miserably)")
                logger.exception(e)

        if pseudo_data_file:
            suffix_aliase = {
                "prefit": "shapes_prefit",
                "postfit": "shapes_fit_s",
            }
            channel_region, suffix = hist_dir.split("/")[0].split("_")
            pseudo_data_hist_dir = suffix_aliase[suffix]+"/"+channel_region+"/data"
            h_data = pseudo_data_file.Get(str(pseudo_data_hist_dir))

        if "tgraph" in str(type(h_data)).lower():
            x = np.array([1.0])
            y = np.array([1.0])
            h_data.GetPoint(1, x, y)
            binwidth = x[0]
            h_data.GetPoint(0, x, y)
            binwidth -= x[0]
            minX = x[0] - binwidth / 2
            h_data.GetPoint(h_data.GetN() - 1, x, y)
            maxX = x[0] - binwidth / 2 + binwidth
            h_data_temp = ROOT.TH1F("data", "data", h_data.GetN(), minX, maxX)
            for i in range(h_data.GetN()):
                h_data.GetPoint(i, x, y)

                h_data_temp.SetBinContent(i + 1, y)
                h_data_temp.SetBinError(i + 1, h_data.GetErrorY(i))
            h_data = h_data_temp
            if pseudo_data_file:
                h_data.Scale(binwidth)
        if rebin_factor > 1:
            h_data = h_data.Rebin(rebin_factor)
        if new_binning is not None:
            h_data = h_data.Rebin(len(new_binning) - 1, "", new_binning)
            # h_data.GetYaxis().SetTitle(
            #     "Events / %i GeV" % (int(new_binning[1] - new_binning[0]))
            # )
            h_data.GetYaxis().SetTitle(
                "Events"
            )
        elif yTitle is not None:
            h_data.GetYaxis().SetTitle(yTitle)

        h_qcd_from_data = h_data.Clone()

    # mc_hists = {}
    mc_hists = collections.OrderedDict()

    for sample in samples:
        this_hist = None

        def get_sample_hist(sample_name, return_empty=False):
            this_subhist_ = None
            try:
                this_subhist_ = f_hists.Get(str(hist_dir % sample_name)).Clone()
            except BaseException as e:
                this_subhist_ = None
                logger.exception(
                    ("hist " + str(hist_dir % sample_name) + " not found. Trying adding the year as suffix.")
                )
                logger.exception(e)
                sample_name_year = "{}_{}".format(sample_name, year)
                try:
                    this_subhist_ = f_hists.Get(str(hist_dir % sample_name_year)).Clone()
                except BaseException as e:
                    this_subhist_ = None
                    logger.exception(
                        " hist" + str(hist_dir % sample_name_year) + "not found. Taking empty hist or None instead."
                    )
                    logger.exception(e)
            return this_subhist_

        if sample in merged_hists and include_merged_hists:
            for subsample in merged_hists[sample]:
                this_subhist = get_sample_hist(subsample)
                if this_subhist is None:
                    this_subhist = h_data.Clone()
                    this_subhist.Reset()
                if this_hist is None:
                    this_hist = this_subhist.Clone()
                    this_hist.SetTitle(this_hist.GetTitle().replace(subsample, sample))
                    this_hist.SetName(this_hist.GetName().replace(subsample, sample))
                else:
                    this_hist.Add(this_subhist.Clone())
        else:
            this_hist = get_sample_hist(sample)
        if this_hist is None:
            logger.critical("Histogram of sample {} could not be found. Skipping.".format(sample))
            continue
        if rebin_factor > 1:
            this_hist = this_hist.Rebin(rebin_factor)
        if new_binning is not None:
            this_hist = this_hist.Rebin(len(new_binning) - 1, "", new_binning)
            # this_hist.GetYaxis().SetTitle(
            #     "Events / %i GeV" % (int(new_binning[1] - new_binning[0]))
            # )
            this_hist.GetYaxis().SetTitle(
                "Events"
            )
        elif yTitle is not None:
            this_hist.GetYaxis().SetTitle(yTitle)

        if "QCD" not in sample and not pseudo_data:
            h_qcd_from_data.Add(this_hist, -1.0)

        mc_hists.update({sample: this_hist})

    if pseudo_data:
        for sample_name, h in mc_hists.items():
            if h_data is None:
                h_data = h.Clone()
                h_data.SetName(h_data.GetName().replace(sample_name, "PseudoData"))
            else:
                h_data.Add(h, 1.0)
    else:
        if scaleQCD and (selection in ["W", "Zbb"]):
            norm = mc_hists["QCD"].Integral()
            global qcd_scale
            qcd_scale = (h_qcd_from_data.Integral() / norm) if norm > 0 else 1.0
            mc_hists["QCD"].Scale(qcd_scale)

    return (h_data, mc_hists)


def plot_data_mc(
    h_data=None,
    h_mc=None,
    plot_title="",
    out_dir=None,
    legend_entries=[],
    additional_text="",
    signal_mc=[],
    additional_hists=[],
    additional_data=[],
    cutflow=False,
    ratio_y_range=(0.75, 1.25),
    contains_data=True,
):

    cms_style.bottom_right_margin_modifier = 1.5 if cutflow else 1.0
    # cms_style.cms_style()
    if h_data is None and h_mc is None and len(additional_hists) == 0:
        raise ValueError("No Histograms were provided!")
    c_width = 1200 if ultrawide else (900 if cutflow else 600)
    # print(plot_title)
    c = ROOT.TCanvas(plot_title, plot_title, c_width, 600)
    legend = (
        ROOT.TLegend(0, 0, 1, 1)
        if legend_on_extern_canvas
        else ROOT.TLegend(*legend_bbox)
    )
    legend.SetFillStyle(0)
    if legend_columns > 1:
        legend.SetNColumns(legend_columns)
    cms_style.ratio_plot = ratio_plot
    cms_style.additional_pad = len(signal_mc) > 0
    plotpad, ratiopad, additional_pad = cms_style.setup_pads(c, logY=logY)

    minY = y_range[0]
    maxY = y_range[1]

    # if(private_work):
    #     cms_style.cms_text = ""
    #     cms_style.extra_text = "private work"
    # else:
    cms_style.extra_text = extra_text
    cms_style.extra_text_rel_X = 0.1
    cms_style.font_size_modifier = 0.6
    # if ratio_plot:
    cms_style.text_padding = lumi_text_padding
    # cms_style.text_padding = lumi_text_padding

    cms_style.draw_lumi(
        plotpad,
        year=year_alias.get(year, year),
        lumi=lumis[year],
        do_extra_text=draw_extra_text,
        out_of_frame=True,
        do_cms_text=draw_extra_text,
        private_work=private_work,
        data=contains_data
    )
    plotpad.cd()

    bkg_err = None
    bkg_stack = None

    if type(h_mc) is ROOT.THStack:
        bkg_stack = h_mc
        for h in h_mc:
            if bkg_err is None:
                bkg_err = h.Clone()
            else:
                bkg_err.Add(h)
    elif type(h_mc) is ROOT.TH1F:
        bkg_stack = h_mc.Clone()
        bkg_err = h_mc.Clone()
    elif type(h_mc) is list or type(h_mc) in [dict, collections.OrderedDict]:
        bkg_stack = ROOT.THStack()
        bkg_err = None
        h_mc_list = h_mc if type(h_mc) is list else h_mc.values()
        for h in h_mc_list:
            bkg_stack.Add(h.Clone(), "Hist")
            if bkg_err is None:
                bkg_err = h.Clone()
            else:
                bkg_err.Add(h)
    elif h_mc is None:
        pass
    else:
        raise NotImplementedError(
            "Type conversion for MC Hist(s) is not implementd (provided type: "
            + str(type(h_mc))
            + ")"
        )

    for entry in reversed(legend_entries):
        legend.AddEntry(*entry)
    if bkg_err is not None:
        bkg_err.SetFillStyle(MC_stat_err_fillstyle)
        bkg_err.SetFillColor(MC_stat_err_color)
        bkg_err.SetLineWidth(0)
        bkg_err.SetMarkerSize(0)
        legend.AddEntry(bkg_err, "Simul. stat. #otimes syst. unc.", "f")

    frame_hist = (None, None)
    if h_data is not None:
        h_data.SetLineColor(obs_line_color)
        h_data.SetMarkerStyle(obs_marker_style)
        h_data.SetMarkerSize(obs_marker_size)
        h_data.SetTitle(plot_title)
        h_data.Draw(obs_draw_option)
        frame_hist = (h_data, obs_draw_option)
        h_data.GetXaxis().SetTitle(xTitle)
        cms_style.setup_hist(h_data)
        if bkg_stack is not None:
            bkg_stack.Draw(draw_option + "SAME")
    else:
        if bkg_stack is not None:
            if maxY is not None:
                bkg_stack.SetMaximum(maxY)
            if minY is not None:
                bkg_stack.SetMinimum(minY)

            bkg_stack.Draw(draw_option)
            frame_hist = (bkg_stack.GetHistogram(), draw_option)
            # frame_hist[0].Draw(draw_option)
            # for i in range(1,frame_hist[0].GetNbinsX()+1):
            #     print('changing label:')
            #     print(frame_hist[0].GetXaxis().GetBinLabel(i))
            #     frame_hist[0].GetXaxis().ChangeLabel(i,30.)

            # frame_hist[0].Draw(draw_option)
            # bkg_stack.Draw(draw_option+"SAME")
            # frame_hist = (bkg_stack.GetHistogram(),draw_option)
            # frame_hist = (bkg_stack,draw_option)
            # bkg_stack.GetYaxis().SetRangeUser(0.9*pow(10,round(np.log10(max_val))-4),1.4*max_val)
            # plotpad.RedrawAxis()

            cms_style.setup_hist(bkg_stack)
        else:
            first_add_hist = additional_hists.pop(0)
            first_add_hist.Draw(draw_option)
            frame_hist = (first_add_hist, draw_option)

    max_val = frame_hist[0].GetMaximum()

    if h_data is not None:
        max_val = max(max_val, h_data.GetMaximum())
        # max_val = h_data.GetMaximum()
    if bkg_err is not None:
        max_val = max(max_val, bkg_err.GetMaximum())

    for h_add_data in additional_data:
        max_val = max(max_val, h_add_data.GetMaximum())
        h_add_data.Draw(obs_draw_option + "SAME")

    if bkg_err is not None:
        bkg_err.Draw("E2SAME")
    for h in additional_hists:
        max_val = max(max_val, h.GetMaximum())
        h.Draw(draw_option + "SAME")

    if frame_hist is not None:
        if None in y_range:
            minY = 0.9 * pow(10, round(np.log10(max_val)) - 4) if logY else 0.0
            maxY = 1.4 * max_val
        # frame_hist[0].GetYaxis().SetRangeUser(minY, maxY)
        # for i in range(1,frame_hist[0].GetNbinsX()+1):
        #     print('changing label:')
        #     print(frame_hist[0].GetXaxis().GetBinLabel(i))
        #     frame_hist[0].GetXaxis().ChangeLabel(i,30.)
        #     frame_hist[0].GetXaxis().ChangeLabel(i,-1,-1,-1,-1,-1,"#pi")

        frame_hist[0].Draw(frame_hist[1] + "SAME")

    if legend_on_extern_canvas:
        legend_canvas = ROOT.TCanvas("legend", "legend", 100, 150)
        legend_canvas.cd()
        legend.Draw()
        legend_canvas.SaveAs(out_dir + "/legend.pdf")
        plotpad.cd()
    else:
        # legend.Draw("SAME")
        legend.Draw()

    if additional_pad is not None:

        bkg_names = [h.GetName() for h in bkg_stack]
        bkg_dict = dict(zip(bkg_names, [h for h in bkg_stack]))
        signal_dict = {}

        for signal_name in signal_mc:
            for bkg_name in bkg_names:
                if signal_name in bkg_name:
                    signal_dict.update({signal_name: bkg_dict[bkg_name]})
                    bkg_dict.pop(bkg_name)

        h_data_minus_bg = h_data.Clone()
        for h in bkg_dict.values():
            h_data_minus_bg.Add(h, -1.0)

        h_signal = ROOT.THStack()
        for signal_name in signal_mc:
            h_signal.Add(signal_dict[signal_name], "HIST")

        additional_pad.cd()

        cms_style.setup_add_hist(h_data_minus_bg)
        h_data_minus_bg.Draw(obs_draw_option)
        h_data_minus_bg.GetYaxis().SetTitle("Data-Bkg")
        h_data_minus_bg.GetYaxis().SetRangeUser(-1, 1.1 * h_signal.GetMaximum())
        h_signal.Draw(draw_option + "SAME")
        h_data_minus_bg.Draw(obs_draw_option + "SAME")

    if ratio_plot:
        def setup_yrange(h):
            y_min_ratio, y_max_ratio = y_range_ratio

            ratio_hist_max = 1.
            for ibin in range(1, h.GetNbinsX()+1):
                binc = h.GetBinContent(ibin)
                binerr = h.GetBinError(ibin)
                if binc == 0:
                    continue
                ratio_hist_max = max(abs(1-ratio_hist_max), abs(1-binc)) + 1
                ratio_hist_max = max(abs(1-ratio_hist_max), abs(binerr)) + 1
            if y_min_ratio is None:
                y_min_ratio = ratio_hist_max
            if y_max_ratio is None:
                y_max_ratio = ratio_hist_max
            # make symmetric
            yrange = max(abs(1-y_min_ratio), abs(1-y_max_ratio))
            return 1-yrange, 1+yrange

        ratiopad.cd()
        ratio_hist = []
        if bkg_err is not None:
            ratio_hist.append(h_data.Clone())
            ratio_hist[0].Divide(bkg_err)

            ratio_hist[0].SetLineColor(1)
            ratio_hist[0].SetMarkerColor(1)
            ratio_hist[0].SetLineWidth(2)
            ratio_hist[0].SetMarkerStyle(8)
            ratio_hist[0].SetMarkerSize(0.5)
            cms_style.setup_ratio_hist(ratio_hist[0])
            ratio_y_min, ratio_y_max = ratio_y_range
            ratio_hist[0].GetYaxis().SetRangeUser(ratio_y_min, ratio_y_max)
            ratio_hist[0].Draw("PE1X0")
            arrows = []
            arrow_length = (1 - ratio_y_min) / 5
            arrow_color = ROOT.kRed - 7  # 13
            for ibin in range(1, ratio_hist[0].GetNbinsX() + 1):
                bin_content = ratio_hist[0].GetBinContent(ibin)
                bin_center = ratio_hist[0].GetBinCenter(ibin)
                if bin_content > ratio_y_max:
                    arrows.append(
                        ROOT.TArrow(bin_center, ratio_y_max - arrow_length, bin_center, ratio_y_max, 0.01, "|>")
                    )
                elif bin_content < ratio_y_min and bin_content > 0:
                    arrows.append(
                        ROOT.TArrow(bin_center, ratio_y_min + arrow_length, bin_center, ratio_y_min, 0.01, "|>")
                    )
                for arrow in arrows:
                    arrow.SetLineColor(arrow_color)
                    arrow.SetFillColor(arrow_color)
                    arrow.Draw()

            ratio_stat_err = bkg_err.Clone()
            ratio_stat_err.Divide(bkg_err)
            ratio_stat_err.SetFillStyle(MC_stat_err_fillstyle)
            ratio_stat_err.SetFillColor(MC_stat_err_color)
            ratio_stat_err.SetLineColor(MC_stat_err_color)
            ratio_stat_err.SetMarkerSize(0)
            ratio_stat_err.SetLineStyle(ROOT.kDashed)

            ratio_stat_err.Draw("E2SAME")
            ratio_hist[0].Draw("PE1X0SAME")
        else:
            for h in additional_hists:
                ratio_hist.append(h.Clone())
                ratio_hist[-1].Divide(h_data)

            cms_style.setup_ratio_hist(ratio_hist[0])
            ratio_hist[0].GetYaxis().SetRangeUser(*setup_yrange(ratio_hist[0]))
            ratio_hist[0].Draw(ratio_draw_options)
            for h in ratio_hist[1:]:
                h.Draw(ratio_draw_options + "SAME")

        for h in ratio_hist:
            h.GetXaxis().SetTitle(h_data.GetXaxis().GetTitle())
            h.GetYaxis().SetTitle(ratio_hist_yTitle)

        ratioXMin = ratio_hist[0].GetXaxis().GetXmin()
        ratioXMax = ratio_hist[0].GetXaxis().GetXmax()

        zeropercent = ROOT.TLine(ratioXMin, 1, ratioXMax, 1)
        plus10percent = ROOT.TLine(ratioXMin, 1.1, ratioXMax, 1.1)
        minus10percent = ROOT.TLine(ratioXMin, 0.9, ratioXMax, 0.9)
        plus10percent.SetLineStyle(ROOT.kDashed)
        minus10percent.SetLineStyle(ROOT.kDashed)

        zeropercent.Draw()
        # plus10percent.Draw()
        # minus10percent.Draw()

        plotpad.cd()
    additional_latex = ROOT.TLatex()
    additional_latex.SetNDC(1)
    additional_latex.SetTextFont(42)
    additional_latex.SetTextSize(
        cms_style.additional_text_size * additional_text_size_modifier
    )
    # additional_latex.SetTextAlign(31)
    additional_latex.SetTextAlign(11)
    for i in range(len(additional_text.split("\\"))):
        additional_latex.DrawLatex(
            plotpad.GetLeftMargin() * (1 + additional_text_padding),
            cms_style.additional_text_ypos
            - 1.1 * i * cms_style.additional_text_size * additional_text_size_modifier,
            additional_text.split("\\")[i],
        )

    plotpad.RedrawAxis()
    c.RedrawAxis()
    if out_dir is not None:
        for file_ext in [".pdf", ".png"]:
            c.SaveAs(out_dir + "/" + plot_title.replace(" ", "_") + file_ext)
        # c.SaveAs(out_dir+'/'+plot_title.replace(" ","_")+'.png')
    return (
        c,
        (plotpad, ratiopad, additional_pad),
        plot_title.replace(" ", "_") + ".pdf",
    )


if __name__ == "__main__":
    pass
