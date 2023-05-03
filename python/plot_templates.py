#!/usr/bin/env python
from __future__ import print_function
import UHH2.JetMass.plotter as plotter
import ROOT
import os

response_samples = {
    "W": [
        "ST_tW_top",
        "ST_tW_antitop",
        "TTToSemiLeptonic",
        "TTToHadronic",
        "ZJetsUnmatched",
        "ZJetsMatched",
        "WJetsUnmatched",
        "WJetsMatched",
        "QCD",
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
    "top": [
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
}

# pt_bins_dict = {
#     "top": ["inclusive", "300to400", "400to500", "500toInf",
#             "350to400", "500to600", "600to1000", "350to1000"],

#     "W": ["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf",
#           "500to600", "600to750", "750to900", "900to1200", "500to1200"],

#     "WfromTop": ["inclusive", "200to300", "300to400", "400to500", "500toInf",
#                  "200to225", "225to250", "250to275", "275to300", "300to325", "325to350",
#                  "350to400", "400to500", "200to350"],
#     # "top": ["400to500"], "W": ["500to1200"], "WfromTop": ["225to250"],
# }

xTitles = {
    # 'PFMass':       'm_{#sum PF-Constituents}',
    # 'PFMass':       'm_{SD}^{*}',
    "mjet": "",
}

response_xTitles = {
    "response_mjet_mgenparticles": "m_{#sum PF-Constituents}/m_{#sum Gen-Constitutents}",
    "response_mjet_mgensubjets": "m_{#sum PF-Constituents}/m_{#sum subjets, gen}",
    # "response_mjet_SD_mgensubjets": 'm_{SD}/m_{#sum subjets, gen}',
    "response_mjet_SD_mgensubjets": "m_{SD, reco}/m_{SD, gen}",
    "response_mtopjet_mgentopjet": "m_{AK8, reco}/m_{AK8, gen}",
    "response_msubjets_mgensubjets": "m_{#sum subjets, reco}/m_{#sum subjets, gen}",
    "response_mjet_mjet_SD": "m_{#sum PF-Constituents}/m_{SD}",
}

selection_text = {
    "W": "V+jets Selection",
    "WfromTop": "W Selection",
    "top": "t#bar{t}\rightarrow #mu + jets Selection",
}

scaleQCD = False
reverse_stacking = False

# plotter.legend_on_extern_canvas = True
plotter.draw_extra_text = True
plotter.private_work = True


def plot_mass(
    selection,
    mass_name,
    hist_file_path,
    output_dir="test",
    binning="CMS",
    logY=False,
    signal_mc=[],
    scaleQCD=True,
):

    f_hists = ROOT.TFile(hist_file_path, "READ")
    ROOT.TH1.AddDirectory(0)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pt_bins = plotter.pt_bins_dict[selection]["jms"]
    nbjet_bins = plotter.nbjet_bins_dict[selection]

    for ipt_bin in range(len(pt_bins)):
        for inbjet_bin in range(len(nbjet_bins)):
            pt_bin = pt_bins[ipt_bin]
            nbjet_bin = nbjet_bins[inbjet_bin]

            regions = (
                ["pass", "passW", "fail"]
                if selection.lower() == "top"
                else ["pass", "fail"]
            )

            for region in regions:
                legend_entries = []
                hist_dir = (
                    selection
                    + "_%s__"
                    + mass_name
                    + "_"
                    + pt_bin
                    + ("" if nbjet_bin == "" else "_" + nbjet_bin)
                    + ("_" + region if len(region) > 0 else "")
                )
                h_data, mc_hists = plotter.get_hists(
                    f_hists,
                    plotter.mc_samples[selection],
                    hist_dir,
                    selection,
                    scaleQCD=scaleQCD,
                )
                h_data.GetXaxis().SetTitle("m_{SD}")
                bkg_stack = ROOT.THStack()
                for sample in plotter.mc_samples[selection]:
                    mc_hists[sample].SetFillColorAlpha(plotter.colors.get(sample), 0.8)
                    mc_hists[sample].SetLineColor(plotter.colors.get(sample))
                    # mc_hists[sample].SetLineColor(plotter.colors[sample])
                    # mc_hists[sample].SetFillColor(plotter.colors[sample])
                    legend_entries.append(
                        (mc_hists[sample], plotter.legend_names[sample], "f")
                    )
                    bkg_stack.Add(mc_hists[sample], "Hist")

                legend_entries.append((h_data, "Data", "p"))

                plotter.logY = logY
                bininfo_tex = plotter.pt_bins_tex_dict[pt_bin] + (
                    ""
                    if nbjet_bin == ""
                    else "\\" + plotter.nbjet_bins_tex_dict[nbjet_bin]
                )
                additional_text = (
                    " "
                    + plotter.selection_tex[selection]
                    + "\\"
                    + bininfo_tex
                    + "\\ %s #bf{%s}"
                    % (plotter.region_tex[selection][region], "prefit")
                )

                # plotter.plot_data_mc(h_data,bkg_stack,hist_dir.replace("_%s",""),out_dir=output_dir,legend_entries=legend_entries,additional_text=plotter.pt_bins_tex_dict[pt_bin])
                plotter.plot_data_mc(
                    h_data,
                    bkg_stack,
                    hist_dir.replace("_%s", ""),
                    out_dir=output_dir,
                    legend_entries=legend_entries,
                    additional_text=additional_text,
                )


if __name__ == "__main__":

    # hist_file = "../macros/Histograms_oneScale_JER_0p005_MIDDB_0p7.root"
    # hist_file = "../macros/Histograms_TopPtReweighting_NBjets.root"
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
    # hist_file = "../macros/Histograms_%s.root"%args.year
    hist_file = args.input
    print(hist_file)

    # for selection in ['top','W']:
    # for selection in ['top','W','Zbb']:
    # for selection in ['Zbb']:
    output_dir = (
        (
            "../Plots/%s_coffea%s/"
            % (
                args.input.replace(".root", ""),
                "_originalQCDScale" if args.origQCDScale else "",
            )
        )
        if args.output == ""
        else args.output
    )
    for selection in ["W", "top"]:
        print(selection, hist_file)
        plot_mass(
            selection,
            "mjet",
            hist_file,
            output_dir,
            binning="CMS",
            scaleQCD=not args.origQCDScale,
        )
    # for selection in ['W']:
    #     plot_mass(selection,'mjet',hist_file,'../Plots/templates',binning="CMS")
    # hist_file = "../macros/Histograms.root"
    # for selection in ['top']:
    #     plot_mass(selection,'mjet',hist_file,'../Plots/templates',binning="CMS")
