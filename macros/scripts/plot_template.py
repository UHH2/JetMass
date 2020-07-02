import ROOT
from ROOT import gROOT,gStyle
import argparse, os, json
import numpy as np
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)

gStyle.SetTextFont(43)

gStyle.SetTitleOffset(0.86,"X")
gStyle.SetTitleOffset(1.8,"Y")
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
min_MSD = 0
max_MSD = 250
n_MSD = 25

ultrawide = False
legend_on_extern_canvas = True

new_binning = np.linspace(min_MSD,max_MSD,n_MSD+1) if n_MSD > 0 else None
rebin_MSD = n_MSD > 0
print(new_binning)
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
    'ST_tWch':41,
    'ST_tch':43,
    'ST_tW_top':41,
    'ST_tW_antitop':42,
    'ST_tWch_top':41,
    'ST_tWch_antitop':42,
    'ST_tch_top':43,
    'ST_tch_antitop':44,
    'ST_sch':40,    
    'TTbar_Hadronic':810,
    'TTbar_SemiLeptonic':804,
    'TTToHadronic':810,
    'TTToSemiLeptonic':804,
    'TTbar_had':810,
    'TTbar_semilep':804,
    'TTbar_dilep':803,
    'other':867
}


samples = {
    # "W":["ST_tW_top","ST_tW_antitop","TTbar_SemiLeptonic","TTbar_Hadronic", "ZJetsUnmatched","ZJetsMatched", "WJetsUnmatched","WJetsMatched", "QCD"],
    "W":["ST_tW_top","ST_tW_antitop","TTToSemiLeptonic","TTToHadronic", "ZJetsUnmatched","ZJetsMatched", "WJetsUnmatched","WJetsMatched", "QCD"],
    "WfromTop":["QCD", "ST_sch", "ST_tch","ST_tWch", "DYJets", "WJets", "TTbar_had","TTbar_dilep", "TTbar_semilep"],
}
pt_bins_dict = {
    "W":["inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf"],
    "WfromTop":["inclusive", "200to300", "300to400", "400to500", "500toInf"],
}


obs_draw_option = 'PE1'
draw_option = 'H'

ratio_plot = True
xLabelSize=18.
yLabelSize=18.
xTitleSize=20.
yTitleSize=22.
xTitleOffset=2.8
yTitleOffset=1.5

logX=False
logY=True

yplot=0.7
yratio=0.3
ymax=1.0
xmax=1.0
xmin=0.0

y_range = [None,None]
# y_range = [0,2*10**4]
x_range = [None,None]
YRangeUser = all(list(map(lambda a: a is not None ,y_range)))
XRangeUser = all(list(map(lambda a: a is not None ,x_range)))

def setup_hist(hist):
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(yTitleSize)
    hist.GetYaxis().SetTitleOffset(yTitleOffset)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(yLabelSize)
    if(ratio_plot):
        hist.GetXaxis().SetTitleSize(0.0)
        hist.GetXaxis().SetLabelSize(0.0)
    else:
        hist.GetXaxis().SetTitle(xTitle)
        hist.GetXaxis().SetTitleFont(43)
        hist.GetXaxis().SetTitleSize(xTitleSize)
        hist.GetXaxis().SetTitleOffset(xTitleOffset)
        hist.GetXaxis().SetLabelFont(43)
        hist.GetXaxis().SetLabelSize(xLabelSize)

    if(YRangeUser):
        hist.GetYaxis().SetRangeUser(y_range[0],y_range[1])
    if(XRangeUser):
        hist.GetXaxis().SetRangeUser(x_range[0],x_range[1])
                
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
    if(ratio_plot):
        plotpad = ROOT.TPad("plotpad","Plot",xmin,ymax-yplot,xmax,ymax)
        ratiopad = ROOT.TPad("ratiopad","Ratio",xmin,ymax-yplot-yratio,xmax,ymax-yplot)
    else:
        plotpad = ROOT.TPad("plotpad","Plot",xmin,ymax-yplot-yratio,xmax,ymax)
        ratiopad = None

    plotpad.SetTopMargin(0.08)
    plotpad.SetLeftMargin(leftMargin)
    plotpad.SetRightMargin(0.05)
    plotpad.SetTicks()
    plotpad.Draw()
    
    if(ratio_plot):
        plotpad.SetBottomMargin(0.016)
        ratiopad.SetTopMargin(0.016)
        ratiopad.SetBottomMargin(0.35)
        ratiopad.SetLeftMargin(leftMargin)
        ratiopad.SetRightMargin(0.05)
        ratiopad.SetTicks()
        ratiopad.Draw()
    else:
        plotpad.SetBottomMargin(0.1)
        
    if(logY):
        plotpad.SetLogy()
        c.SetLogy()
    if(logX):
        plotpad.SetLogx()
        if(ratio_plot):
            ratiopad.SetLogx()
        c.SetLogx()
    return plotpad,ratiopad

def plot_template(selection = 'W',hist_file_path = "../Histograms.root",out_dir='../../Plots/Templates',pseudo_data=""):
    f_hists = ROOT.TFile(hist_file_path,"READ")
    ROOT.TH1.AddDirectory(0)

    if(not os.path.exists(out_dir)):
        os.makedirs(out_dir)

    pt_bins = pt_bins_dict[selection]
    pt_bins_str = [(pt_bin.replace('to',' #geq p_{T} > ')+('' if 'in' in pt_bin.lower() else ' GeV')).replace(' > Inf','') for pt_bin in pt_bins]
    
    mc_samples = samples[selection]

    ext_legend_canvas_done = False
    
    for ipt_bin in range(len(pt_bins)):
        pt_bin = pt_bins[ipt_bin]
        
        for region in ['pass','fail']:
           
            plot_title = '%s %s'%(pt_bins_str[ipt_bin],region)

            c_width = 1800 if ultrawide else 600 
            
            c = ROOT.TCanvas(plot_title,plot_title,c_width,600)            
            legend = ROOT.TLegend(0,0,1,1) if legend_on_extern_canvas else ROOT.TLegend(0.70,0.6,0.9,0.9)
            legend.SetFillStyle(0)
            plotpad,ratiopad = setup_pads(c)              
            plotpad.cd()
            
            hist_dir = selection+'_%s__mjet_'+pt_bin+'_'+region
            if(pseudo_data==""):
                h_data = f_hists.Get(hist_dir%"Data")
            else:
                if("nominal" in pseudo_data):
                    hist_dir_pseudo = hist_dir
                else:
                    hist_dir_pseudo = selection+'_%s__mjet_'+pseudo_data[0]+pt_bin+'_'+region+'__'+pseudo_data[1] 
                h_data = f_hists.Get(hist_dir_pseudo%mc_samples[0])
                print(hist_dir_pseudo%mc_samples[0])
                for sample in mc_samples[1:]:
                    print(hist_dir_pseudo%sample)
                    h_data.Add(f_hists.Get(hist_dir_pseudo%sample))
            if(rebin_MSD):
                h_data = h_data.Rebin(n_MSD,"",new_binning)            
            legend.AddEntry(h_data,'Data','p')

            h_qcd_from_data = h_data.Clone()
            
            bkg_stack = ROOT.THStack()
            bkg_Err = None

            for sample in mc_samples:
                if('QCD' in sample):
                    continue
                this_hist = f_hists.Get(hist_dir%sample)
                if(rebin_MSD):
                    this_hist = this_hist.Rebin(n_MSD,"",new_binning)                
                    h_qcd_from_data.Add(this_hist,-1)
                    
            for sample in mc_samples:
                this_hist = None
                # if('QCD' in sample):                    
                #     this_hist = f_hists.Get(hist_dir%'QCD_lowPt')
                #     this_hist.Add(f_hists.Get(hist_dir%'QCD_highPt'))
                # else:
                this_hist = f_hists.Get(hist_dir%sample)
                if(rebin_MSD):
                    this_hist = this_hist.Rebin(n_MSD,"",new_binning)

                #scale QCD to match Data-W-Z (Integral)
                if('QCD' in sample):
                    norm = this_hist.Integral()
                    QCD_scale = (h_qcd_from_data.Integral()/norm) if norm > 0 else 1.0
                    this_hist.Scale(QCD_scale)
                
                this_hist.SetLineColor(colors[sample])
                this_hist.SetFillColor(colors[sample])
                legend.AddEntry(this_hist,sample,'f')
                if(bkg_Err is None):
                    bkg_Err = this_hist.Clone()
                else:
                    bkg_Err.Add(this_hist)
                bkg_stack.Add(this_hist,'Hist')

            bkg_Err.SetFillStyle(3204);
            bkg_Err.SetFillColor(922);
            bkg_Err.SetLineColor(1);
            legend.AddEntry(bkg_Err,"MC stat. Unc.","f");

            h_data.SetLineColor(1)
            h_data.SetMarkerStyle(8)
            h_data.SetMarkerSize(0.5)
            h_data.SetTitle(plot_title)
            h_data.Draw(obs_draw_option)
            
            h_data.GetYaxis().SetRangeUser(0.9*pow(10,round(np.log10(h_data.GetMaximum()))-4),1.4*h_data.GetMaximum())
            setup_hist(h_data)
            
            bkg_stack.Draw(draw_option + 'SAME')
            bkg_Err.Draw("E2SAME")
            h_data.Draw(obs_draw_option+'SAME')

            if(legend_on_extern_canvas):
                if(not ext_legend_canvas_done):
                    legend_canvas = ROOT.TCanvas('legend','legend',100,150)
                    legend_canvas.cd()
                    legend.Draw()
                    legend_canvas.SaveAs(out_dir+'/legend.pdf')
                    ext_legend_canvas_done=True
                    plotpad.cd()
            else:
                legend.Draw('SAME')

            ratiopad.cd()

            ratio_hist = h_data.Clone()
            ratio_hist.GetYaxis().SetTitle("#frac{Data}{MC}")
            ratio_hist.GetYaxis().SetRangeUser(0.5,1.5)
            ratio_hist.Divide(bkg_Err)
            ratio_hist.SetLineColor(1)
            ratio_hist.SetMarkerColor(1)
            ratio_hist.SetLineWidth(2)
            ratio_hist.SetMarkerStyle(1)
            setup_ratio_hist(ratio_hist)
            ratio_hist.GetYaxis().SetRangeUser(0.3,1.7)
                

            ratio_hist.Draw('EP')

            ratioXMin = ratio_hist.GetXaxis().GetXmin()
            ratioXMax = ratio_hist.GetXaxis().GetXmax()

            zeropercent = ROOT.TLine(ratioXMin,1,ratioXMax,1)
            plus10percent = ROOT.TLine(ratioXMin,1.1,ratioXMax,1.1)
            minus10percent = ROOT.TLine(ratioXMin,0.9,ratioXMax,0.9)
            plus10percent.SetLineStyle(ROOT.kDashed)
            minus10percent.SetLineStyle(ROOT.kDashed)
            
            zeropercent.Draw()
            plus10percent.Draw()
            minus10percent.Draw()

                
            plotpad.cd()
            latex = ROOT.TLatex()
            latex.SetNDC(1)
            latex.SetTextSize(24)
            if(pt_bin is not "inclusive"):
                print(pt_bin)
                latex.DrawLatex(0.20,0.95," %s GeV #leq p_{T} < %s GeV"%(pt_bin.split('to')[0],pt_bin.split('to')[1]))


            c.RedrawAxis()
            c.SaveAs(out_dir+'/'+hist_dir.replace("%s","")+'.pdf')


if(__name__ == '__main__'):
    # plot_template("W","../Histograms.root",'../../Plots/Vqq_jets_templates/')
    plot_template("W","../scaleStudy/Histograms_oneScale.root",'../../Plots/oneScale_0p01_all_down/',pseudo_data=["0_0_all_","down"])
