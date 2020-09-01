import ROOT
from ROOT import gROOT,gStyle
import argparse, os, json
import numpy as np
import sys
import cms_style
cms_style.cms_style()


# gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)


leftMargin = 0.14



ultrawide = False
legend_on_extern_canvas = False


new_binning_dict = {'top':np.linspace(0,350,71),'W':np.linspace(40,160,13),'WfromTop':np.linspace(0,200,40)}
# new_binning_dict = {'top':np.linspace(-10,350,73),'W':np.linspace(-10,160,18),'WfromTop':np.linspace(-10,200,42)}
#myBinning
# new_binning_dict = {'top':np.linspace(0,350,71),'W':np.linspace(0,500,101),'WfromTop':np.linspace(0,300,61)}

binning_dict ={
    "ATLAS": {'top':np.linspace(0,350,71),'W':np.linspace(40,160,13),'WfromTop':np.linspace(0,200,40)},
    # "CMS":{'top':np.linspace(0,350,71),'W':np.linspace(0,500,101),'WfromTop':np.linspace(0,300,61)}
    "CMS":{'top':np.linspace(50,200,16),'W':np.linspace(0,500,101),'WfromTop':np.linspace(0,300,61)}
}

rebin = True
scaleQCD=True

colors = {
    'QCD':867,
    'QCD_lowPt':867,
    'QCD_highPt':867,
    'WJets':413,
    'WUnmatched':419,
    'WMatched':413,
    'WJetsUnmatched':419,
    'WJetsMatched':413,
    'ZJets':800,
    'DYJets':800,
    'ZJetsUnmatched':797,
    'ZJetsMatched':800,
    'ZUnmatched':797,
    'ZMatched':800,
    'SingleTop':800,
    'ST':800,
    'ST_tWch':41,
    'ST_tch':43,
    'ST_tW_top':41,
    'ST_tW_antitop':42,
    'ST_tWch_top':41,
    'ST_tWch_antitop':42,
    'ST_tch_top':43,
    'ST_tch_antitop':44,
    'ST_sch':40,    
    'ttbar':810,
    'TTbar':810,
    'TTbar_Hadronic':810,
    'TTbar_SemiLeptonic':804,
    'TTToHadronic':810,
    'TTToSemiLeptonic':804,
    'TTbar_had':810,
    'TTbar_semilep':804,
    'TTbar_dilep':803,
    'other':921,
    'other_vjets':921,
    'other_ttbar':921,
}


legend_names = {
    'QCD':"QCD",
    'WJets':"W+Jets",
    'WUnmatched':"W+Jets (Unmatched)",
    'WMatched':"W+Jets (Matched)",
    'WJetsUnmatched':"W+Jets (Unmatched)",
    'WJetsMatched':"W+Jets (Matched)",
    'ZJets':"Z+Jets",
    'DYJets':"DY+Jets",
    'ZJetsUnmatched':"Z+Jets (Unmatched)",
    'ZJetsMatched':"Z+Jets (Matched)",
    'ZUnmatched':"Z+Jets (Unmatched)",
    'ZMatched':"Z+Jets (Matched)",
    'ST':"Single top",
    'SingleTop':"Single top",
    'ST_tWch':"Single anti-/top (tW-ch.)",
    'ST_tch':"Single anti-/top (t-ch.)",
    'ST_tW_top':"Single top (tW-ch.)",
    'ST_tW_antitop':"Single antitop (tW-ch.)",
    'ST_tWch_top':"Single top (tW-ch.)",
    'ST_tWch_antitop':"Single antitop (tW-ch.)",
    'ST_tch_top':"Single top (t-ch.)",
    'ST_tch_antitop':"Single antitop (t-ch.)",
    'ST_sch':"Single anti-/top (s-ch.)",    
    'ttbar':"t#bar{t}",
    'TTbar':"t#bar{t}",
    'TTbar_Hadronic':"t#bar{t} (hadronic)",
    'TTbar_SemiLeptonic':"t#bar{t} (semi-leptonic)",
    'TTToHadronic':"t#bar{t} (hadronic)",
    'TTToSemiLeptonic':"t#bar{t} (semi-leptonic)",
    'TTbar_had':"t#bar{t} (hadronic)",
    'TTbar_semilep':"t#bar{t} (semi-leptonic)",
    'TTbar_dilep':"t#bar{t} (di-leptonic)",
    'other':"other",
    'other_ttbar':"other",
    'other_vjets':"other",
    'others':"other"
}


samples = {
    "W":["QCD","ST_tW_top","ST_tW_antitop","TTToSemiLeptonic","TTToHadronic", "ZJetsUnmatched","ZJetsMatched", "WJetsUnmatched","WJetsMatched"],
    "WfromTop":["QCD", "ST_sch", "ST_tch","ST_tWch", "DYJets", "WJets", "TTbar_had","TTbar_dilep", "TTbar_semilep"],
    "top":["QCD", "ST_sch", "ST_tch","ST_tWch", "DYJets", "WJets", "TTbar_had","TTbar_dilep", "TTbar_semilep"],
}

mc_samples = {
    "W":["other_vjets","QCD", "ZJetsUnmatched","ZJetsMatched", "WJetsUnmatched","WJetsMatched"],
    "WfromTop":["other_ttbar", "WJets", "ST", "TTbar"],
    "top":["other_ttbar", "WJets", "ST", "TTbar"],
}

merged_hists={
    "TTbar":["TTbar_had","TTbar_dilep", "TTbar_semilep"],
    "ST":["ST_sch", "ST_tch","ST_tWch"],
    "other_ttbar":["QCD","DYJets"],
    "other_vjets":["ST_tW_top","ST_tW_antitop","TTToSemiLeptonic","TTToHadronic"]
    
}


pt_bins_dict = {
    "top":["inclusive", "300to400", "400to500", "500toInf","300to500","300to350","350to400","400to450","450to500","500to550","550to600","600toInf"],
    "W":["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf",
         "500to600", "600to750", "750to900", "900to1200", "500to1200"],    
    "WfromTop":["inclusive", "200to300", "300to400", "400to500", "500toInf",
                "200to225", "225to250", "250to275", "275to300", "300to325", "325to350", "350to400", "400to500","200to350"],
}

pt_bins_tex_dict = {}
for selection, bins in pt_bins_dict.items():
    for pt_bin in bins:
        pt_bins_tex_dict[pt_bin] = "" if pt_bin is "inclusive" else " %s GeV #leq p_{T} < %s GeV"%(pt_bin.split('to')[0],pt_bin.split('to')[1])



obs_draw_option = 'PE1X0'
draw_option = 'H'
ratio_draw_options='PE1X0'

obs_line_color=1
obs_marker_style=8

extra_text="Preliminary"

ratio_plot = True
ratio_hist_yTitle = "#frac{Data}{MC}"
xLabelSize=18.
yLabelSize=18.
xTitleSize=20.
yTitleSize=22.
xTitleOffset=2.8
yTitleOffset=1.5

logX=False
logY=False

yplot=0.7
yratio=0.3
ymax=1.0
xmax=1.0
xmin=0.0

y_range = [None,None]
x_range = [None,None]
YRangeUser = all(list(map(lambda a: a is not None ,y_range)))
XRangeUser = all(list(map(lambda a: a is not None ,x_range)))                


def get_hists(f_hists, samples, hist_dir='_%s_mjet_inclusive_pass',selection='top',pseudo_data=False):
    hists = []
    hist_dir = selection+hist_dir if hist_dir[0] is '_' else hist_dir

    new_binning = None
    if(rebin):
        new_binning = binning_dict['CMS'].get(selection,None)
        
    h_data = None
    if(not pseudo_data):
        for data_name in ['Data','data_obs']:
            try:
                h_data = f_hists.Get(hist_dir%data_name)
                h_data.GetYaxis().SetTitle("Events / %i GeV"%(int(new_binning[1]-new_binning[0])))
                break
            except:
                print('tried getting data hist with', data_name,'(which failed miserably)')
                
        if(new_binning is not None):
            h_data = h_data.Rebin(len(new_binning)-1,"",new_binning)
        
        h_qcd_from_data = h_data.Clone()
        
    mc_hists = {}
    for sample in samples:
        this_hist = None
        if(sample in merged_hists):
            for subsample in merged_hists[sample]:
                if(this_hist is None):
                    this_hist = f_hists.Get(hist_dir%subsample).Clone()
                    this_hist.SetTitle(this_hist.GetTitle().replace(subsample,sample))
                    this_hist.SetName(this_hist.GetName().replace(subsample,sample))
                else:
                    this_hist.Add(f_hists.Get(hist_dir%subsample).Clone())
        else:                
            this_hist = f_hists.Get(hist_dir%sample).Clone()
                
        if(new_binning is not None):
            this_hist = this_hist.Rebin(len(new_binning)-1,"",new_binning)

        if('QCD' not in sample and not pseudo_data):
            h_qcd_from_data.Add(this_hist,-1.0)

        mc_hists.update({sample:this_hist})

    if(pseudo_data):
        for sample_name,h in mc_hists.items():
            if(h_data is None):
                h_data = h.Clone()
                h_data.SetName(h_data.GetName().replace(sample_name,'PseudoData'))
            else:
                h_data.Add(h,1.0)
    else:
        if(scaleQCD and selection == "W"):
            norm = mc_hists['QCD'].Integral()
            mc_hists['QCD'].Scale((h_qcd_from_data.Integral()/norm) if norm > 0 else 1.0)

    return (h_data,mc_hists)
        
        
def plot_data_mc(h_data=None,h_mc=None,plot_title="",out_dir=None,legend_entries=[],additional_text="",signal_mc=[],additional_hists=[]):
    if(h_data is None and h_mc is None):
        raise ValueError("No Histograms were provided!")
           
    c_width = 1800 if ultrawide else 600 
            
    c = ROOT.TCanvas(plot_title,plot_title,c_width,600)            
    legend = ROOT.TLegend(0,0,1,1) if legend_on_extern_canvas else ROOT.TLegend(0.60,0.6,0.9,0.9)
    legend.SetFillStyle(0)
    cms_style.ratio_plot = True
    cms_style.additional_pad = len(signal_mc)>0
    plotpad,ratiopad,additional_pad = cms_style.setup_pads(c,logY=logY )
            
    cms_style.extra_text = extra_text
    cms_style.extra_text_rel_X = 0.14
    cms_style.font_size_modifier = 0.8
    cms_style.text_padding = 0.4
    cms_style.draw_lumi(plotpad, lumi = 41.8,do_extra_text=True,out_of_frame = True,cms_text=True)
    plotpad.cd()

    bkg_err = None
    bkg_stack = None
    
    if(type(h_mc) is ROOT.THStack):
        bkg_stack = h_mc
        for h in h_mc:
            if(bkg_err is None):
                bkg_err = h.Clone()
            else:
                bkg_err.Add(h)
    elif(type(h_mc) is ROOT.TH1F):
        bkg_stack = h_mc.Clone()
        bkg_err = h_mc.Clone()
    elif(type(h_mc) is list or type(h_mc) is dict):
        bkg_stack = ROOT.THStack()
        bkg_err = None
        h_mc_list = h_mc if type(h_mc) is list else h_mc.values()
        for h in h_mc_list:
            bkg_stack.Add(h.Clone(),"Hist")
            if(bkg_err is None):
                bkg_err = h.Clone()
            else:
                bkg_err.Add(h)
    elif(h_mc is None):
        pass
    else:
        raise NotImplementedError("Type conversion for MC Hist(s) is not implementd (provided type: "+str(type(h_mc))+")")    
        

    for entry in reversed(legend_entries):
        legend.AddEntry(*entry)
    if(bkg_err is not None):
        bkg_err.SetFillStyle(3204)
        bkg_err.SetFillColor(922)
        bkg_err.SetLineWidth(0)
        bkg_err.SetMarkerSize(0)
        legend.AddEntry(bkg_err,"MC stat. Unc.","f")

                # for entry in reversed(legend_mc_entries):
            #     legend.AddEntry(*entry)

    max_val = h_data.GetMaximum()
    if(bkg_err is not None):
        max_val = max(max_val ,bkg_err.GetMaximum())
    # elif(len(additional_hists)>0):
    #     max_val = max(max_val ,additional_hists[0].GetMaximum())
        
    if(h_data is not None):        
        h_data.SetLineColor(obs_line_color)
        h_data.SetMarkerStyle(obs_marker_style)
        h_data.SetMarkerSize(0.5)
        h_data.SetTitle(plot_title)
        h_data.Draw(obs_draw_option)
        h_data.GetYaxis().SetRangeUser(0.9*pow(10,round(np.log10(max_val))-4),1.4*max_val)
        cms_style.setup_hist(h_data)
        if(bkg_stack is not None):
            bkg_stack.Draw(draw_option + 'SAME')
    else:
        if(bkg_stack is not None):
            bkg_stack.Draw(draw_option)
            bkg_stack.GetYaxis().SetRangeUser(0.9*pow(10,round(np.log10(max_val))-4),1.4*max_val)
            cms_style.setup_hist(bkg_stack)
                   
    if(bkg_err is not None):
        bkg_err.Draw("E2SAME")
    for h in additional_hists:
        h.Draw(draw_option+'SAME')

    if(h_data is not None):
        h_data.Draw(obs_draw_option+'SAME')
        
    if(legend_on_extern_canvas):
        legend_canvas = ROOT.TCanvas('legend','legend',100,150)
        legend_canvas.cd()
        legend.Draw()
        legend_canvas.SaveAs(out_dir+'/legend.pdf')
        ext_legend_canvas_done=True
        plotpad.cd()
    else:
        legend.Draw('SAME')

    if(additional_pad is not None):

        bkg_names = [h.GetName() for h in bkg_stack]
        bkg_dict = dict(zip(bkg_names,[h for h in bkg_stack]))
        signal_dict = {}

        for signal_name in signal_mc:
            for bkg_name in bkg_names:
                if signal_name in bkg_name:
                    signal_dict.update({signal_name:bkg_dict[bkg_name]})
                    bkg_dict.pop(bkg_name)


        h_data_minus_bg = h_data.Clone()
        for h in bkg_dict.values():
            h_data_minus_bg.Add(h,-1.0)

        h_signal = ROOT.THStack()
        for signal_name in signal_mc:
            h_signal.Add(signal_dict[signal_name],"HIST")        

        additional_pad.cd()

        cms_style.setup_add_hist(h_data_minus_bg)
        h_data_minus_bg.Draw(obs_draw_option)
        h_data_minus_bg.GetYaxis().SetTitle("Data-Bkg")
        h_data_minus_bg.GetYaxis().SetRangeUser(-1,1.1*h_signal.GetMaximum())
        h_signal.Draw(draw_option+'SAME')
        h_data_minus_bg.Draw(obs_draw_option+'SAME')
        
    ratiopad.cd()
    ratio_hist = []
    if(bkg_err is not None):
        ratio_hist.append(h_data.Clone())    
        ratio_hist[0].Divide(bkg_err)
    
        ratio_hist[0].SetLineColor(1)
        ratio_hist[0].SetMarkerColor(1)
        ratio_hist[0].SetLineWidth(2)
        ratio_hist[0].SetMarkerStyle(8)
        ratio_hist[0].SetMarkerSize(0.5)
        cms_style.setup_ratio_hist(ratio_hist[0])
        ratio_hist[0].GetYaxis().SetRangeUser(0.3,1.7)
        ratio_hist[0].Draw('PE1X0')

        ratio_stat_err = bkg_err.Clone()
        ratio_stat_err.Divide(bkg_err)
        ratio_stat_err.SetFillStyle(3204)
        ratio_stat_err.SetFillColor(922)
        ratio_stat_err.SetLineColor(922)
        ratio_stat_err.SetMarkerSize(0)
        ratio_stat_err.SetLineStyle(ROOT.kDashed)

        ratio_stat_err.Draw('E2SAME')
        ratio_hist[0].Draw('PE1X0SAME')
    else:
        for h in additional_hists:
            ratio_hist.append(h.Clone())
            ratio_hist[-1].Divide(h_data)

        cms_style.setup_ratio_hist(ratio_hist[0])
        ratio_hist[0].GetYaxis().SetRangeUser(0.3,1.7)
        ratio_hist[0].Draw(ratio_draw_options)
        for h in ratio_hist[1:]:
            h.Draw(ratio_draw_options+'SAME')
            

    for h in ratio_hist:
        h.GetXaxis().SetTitle(h_data.GetXaxis().GetTitle())
        h.GetYaxis().SetTitle(ratio_hist_yTitle)
        h.GetYaxis().SetRangeUser(0.5,1.5)

    ratioXMin = ratio_hist[0].GetXaxis().GetXmin()
    ratioXMax = ratio_hist[0].GetXaxis().GetXmax()
    
    zeropercent = ROOT.TLine(ratioXMin,1,ratioXMax,1)
    plus10percent = ROOT.TLine(ratioXMin,1.1,ratioXMax,1.1)
    minus10percent = ROOT.TLine(ratioXMin,0.9,ratioXMax,0.9)
    plus10percent.SetLineStyle(ROOT.kDashed)
    minus10percent.SetLineStyle(ROOT.kDashed)
    
    zeropercent.Draw()
    plus10percent.Draw()
    minus10percent.Draw()

                
    plotpad.cd()
    additional_latex = ROOT.TLatex()
    additional_latex.SetNDC(1)
    additional_latex.SetTextFont(42)
    additional_latex.SetTextSize(cms_style.additional_text_size)
    # additional_latex.SetTextAlign(31)
    additional_latex.SetTextAlign(11)
    for i in range(len(additional_text.split("\\"))):
        additional_latex.DrawLatex(plotpad.GetLeftMargin()*(1+cms_style.text_padding),cms_style.additional_text_ypos-1.1*i*cms_style.additional_text_size,additional_text.split("\\")[i])


    plotpad.RedrawAxis()
    c.RedrawAxis()
    if(out_dir is not None):
        c.SaveAs(out_dir+'/'+plot_title.replace(" ","_")+'.pdf')
    return (c,(plotpad,ratiopad,additional_pad))


    
if(__name__ == '__main__'):
    pass
