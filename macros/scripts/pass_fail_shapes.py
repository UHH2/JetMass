import ROOT
import os
from ROOT import gROOT,gStyle
import numpy as np
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)

gStyle.SetTextFont(43)

gStyle.SetTitleOffset(0.86,"X")
gStyle.SetTitleOffset(1.6,"Y")
# gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadLeftMargin(0.15)
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

ratio_plot = True

min_MSD = 0
max_MSD = 496
# n_MSD = 62
n_MSD = -1
new_binning = np.linspace(min_MSD,max_MSD,n_MSD+1) if n_MSD > 0 else None
rebin_MSD = n_MSD > 0


leftMargin = 0.14
logY=False
logX=False
yplot=0.7
yratio=0.3
ymax=1.0
xmax=1.0
xmin=0.0
xLabelSize=18.
yLabelSize=18.
xTitleSize=20.
yTitleSize=22.
xTitleOffset=2.8
yTitleOffset=1.5
y_range = [None,None]
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

  # ratioHist.GetXaxis().SetTitle("m_{SD} [GeV]")
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

def plot_var(hist_file,dir_suffix,variable,useMC):


    f_hists = ROOT.TFile(hist_file,'read')

    # var_str = {'msd':['Mass_central','m_{SD}'],'pt':['Pt_central_jet','p_{T}'],'rho':['Rho_central','#rho']}
    # var_str = {'msd':['Mass_central_jet','m_{SD}'],'pt':['Pt_central_jet','p_{T}'],'rho':['Rho_central','#rho']}
    var_str = {'msd':['mjet','m_{SD}'],'pt':['Pt_central_jet','p_{T}'],'rho':['Rho_central','#rho']}
    ptbins =[500,550,600,675,800,1200,99999]
    for i in range(len(ptbins)-1):
        tagger = 'N2DDT'
        # tagger = 'DAK8DDT'
        pt_bin_str = ('%ito%i'%(ptbins[i],ptbins[i+1])).replace('99999','Inf')
        hist_name='W_%s__'+var_str[variable][0]+'_'+pt_bin_str+'_%s'
        if(useMC): 
            h_pass = f_hists.Get(hist_name%('QCD','pass'))
            h_fail = f_hists.Get(hist_name%('QCD','fail'))       
        else:
            h_pass = f_hists.Get(hist_name%('Data','pass'))
            h_pass.Add(f_hists.Get(hist_name%('WMatched','pass')),-1)
            h_pass.Add(f_hists.Get(hist_name%('ZMatched','pass')),-1)
            h_pass.Add(f_hists.Get(hist_name%('WUnmatched','pass')),-1)
            h_pass.Add(f_hists.Get(hist_name%('ZUnmatched','pass')),-1)
            h_pass.Add(f_hists.Get(hist_name%('TTbar','pass')),-1)
            h_fail = f_hists.Get(hist_name%('Data','fail'))
            h_fail.Add(f_hists.Get(hist_name%('WMatched','fail')),-1)
            h_fail.Add(f_hists.Get(hist_name%('ZMatched','fail')),-1)
            h_fail.Add(f_hists.Get(hist_name%('WUnmatched','fail')),-1)
            h_fail.Add(f_hists.Get(hist_name%('ZUnmatched','fail')),-1)
            h_fail.Add(f_hists.Get(hist_name%('TTbar','fail')),-1)

        if(rebin_MSD):
            h_pass = h_pass.Rebin(n_MSD,"",new_binning)
            h_fail = h_fail.Rebin(n_MSD,"",new_binning)
            
        norm_pass = h_pass.Integral()
        norm_fail = h_fail.Integral()
        h_pass.Scale((1./norm_pass) if norm_pass>0. else 1)
        h_fail.Scale((1./norm_fail) if norm_fail>0. else 1)

        c = ROOT.TCanvas('c%i'%i,'c%i'%i,700,600)
        legend = ROOT.TLegend(0.60,0.78,0.83,0.90)

        plotpad,ratiopad = setup_pads(c)              
        plotpad.cd()

        setup_hist(h_fail)
        h_pass.SetTitle('QCD shape comparison %i GeV < p_{T} < %i GeV'%(ptbins[i],ptbins[i+1]))
        h_pass.SetLineColor(8)
        h_pass.SetLineWidth(2)
        h_fail.SetLineColor(46)
        h_fail.SetLineWidth(2)
        h_fail.GetYaxis().SetTitle('#DeltaN/N')
        y_range_max = max(0.03,1.1*h_pass.GetMaximum())
        h_fail.GetYaxis().SetRangeUser(0,y_range_max)
        legend.AddEntry(h_pass,'%s pass'%('QCD' if useMC else '(Data-W-Z-TTbar)'),'l')
        legend.AddEntry(h_fail,'%s fail'%('QCD' if useMC else '(Data-W-Z-TTbar)'),'l')
        
        h_fail.Draw('H')
        h_pass.Draw('HSAME')
        latex = ROOT.TLatex()
        latex.SetNDC(1)
        latex.SetTextSize(20)
        latex_str = ("%.0f GeV #leq p_{T} < %.0f GeV"%(ptbins[i],ptbins[i+1])).replace('99999','Inf')
        latex.DrawLatex(0.20,0.80,latex_str)
        
        legend.Draw('SAME')

        ratiopad.cd()
        ratio_hist = h_pass.Clone()
        ratio_hist.Divide(h_fail)
        setup_ratio_hist(ratio_hist)
        ratio_hist.GetYaxis().SetRangeUser(0.3,1.7)
        ratio_hist.GetYaxis().SetTitle("pass/fail")
        ratio_hist.GetXaxis().SetTitle(var_str[variable][1])
        ratio_hist.SetLineColor(1)
        ratio_hist.Draw("EP")

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

        if not os.path.isdir('../../Plots/'+dir_suffix):
            os.makedirs('../../Plots/'+dir_suffix)

        c.SaveAs('../../Plots/'+dir_suffix+'/%s_shape_comp_%s%s.pdf'%(variable,pt_bin_str,"" if useMC else "_fromData"))

def plot_dir(hist_file,dir_name):

    for useMC in [True,False]:
        # for variable in ["msd","pt",'rho']:
        for variable in ["msd"]:
            plot_var(hist_file,dir_name,variable,useMC)
        
if(__name__ == '__main__'):

    plot_dir('../Histograms.root','pass_fail_shapes')
