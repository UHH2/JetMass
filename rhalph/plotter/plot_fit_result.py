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


colors = {
    'qcd':867,
    'QCD':867,
    'WJets':413,
    'WJetsUnmatched':419,
    'WJetsMatched':413,
    'ZJets':800,
    'DYJets':800,
    'ZJetsUnmatched':797,
    'ZJetsMatched':800,
    'WUnmatched':419,
    'WMatched':413,
    'ZUnmatched':797,
    'ZMatched':800,
    'SingleTop':800,
    'ST_tWch':41,
    'ST_tch':43,
    'ST_tWch_top':41,
    'ST_tWch_antitop':42,
    'ST_tW_top':41,
    'ST_tW_antitop':42,
    'ST_tch_top':43,
    'ST_tch_antitop':44,
    'ST_sch':40,    
    'TTbar':810,
    'TTbar_Hadronic':810,
    'TTbar_SemiLeptonic':804,
    'TTToHadronic':810,
    'TTToSemiLeptonic':804,
    'TTbar_had':810,
    'TTbar_semilep':804,
    'TTbar_dilep':803,
    'WJets':413,
    'other':867
}

signal_scale = -1

obs_draw_option = 'PE1'
fit_draw_option = 'H'

xLabelSize=18.
yLabelSize=18.
xTitleSize=20.
yTitleSize=22.
xTitleOffset=2.8
yTitleOffset=1.5

logX=False
logY=False

ratio_plot = True
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

def plot_fit_result(config={'ModelName':'WMassModel'},logY=False):
    f_shapes = ROOT.TFile(config['ModelName']+"/fit_shapes.root","READ")
    ROOT.TH1.AddDirectory(0)

    out_dir = config['ModelName']+'/plots'+('/logY/' if logY else '')
    if(not os.path.exists(out_dir)):
        os.makedirs(out_dir)

    for channel_str, channel in config['channels'].items():
        print('plotting', channel_str)
        signals = channel['signal']
        backgrounds = [bg for bg in channel['samples'] if bg not in signals]
        backgrounds = list(map(lambda bg: 'qcd' if ('QCD' in bg and "QcdEstimation" in channel and channel["QcdEstimation"]=="True") else bg ,backgrounds))
        regions = channel['regions'] if 'regions' in channel else [""]
        for region in regions:
            for suffix in ['prefit','postfit']:
                plot_title = '%s %s %s'%(channel_str,region,suffix)
                c = ROOT.TCanvas(plot_title,plot_title,600,600)
                legend = ROOT.TLegend(0.60,0.65,0.9,0.9)
                legend.SetFillStyle(0)
                plotpad,ratiopad = setup_pads(c)              
                plotpad.cd()

                hist_dir = channel_str + region + '_' + suffix 
                h_obs = f_shapes.Get(hist_dir+'/data_obs')
                total_b = h_obs.Clone()
                total_b.Reset()
                total_sb = h_obs.Clone()
                total_sb.Reset()

                shapes_signal = [f_shapes.Get(hist_dir+'/'+signal) for signal in signals]
                minimum_signal_maximum=1e10
                # for s_sig in shapes_signal:
                #     minimum_signal_maximum = s_sig.GetMaximum() if s_sig.GetMaximum() > minimum_signal_maximum else minimum_signal_maximum

                #     print('max:',s_sig.GetMaximum())
                shapes_background = [f_shapes.Get(hist_dir+'/'+background) for background in backgrounds]

                minima = [hist.GetMinimum() for hist in shapes_signal]+[hist.GetMinimum() for hist in shapes_background]
                while 0.0 in minima:
                    minima.remove(0.0)
                if(len(minima)==0):
                    minima.append(1.0)
                Y_minimum = min(minima)

                h_obs.SetLineColor(1)
                h_obs.SetMarkerStyle(8)
                h_obs.SetTitle(plot_title)
                h_obs.Draw(obs_draw_option)
                h_obs.GetYaxis().SetRangeUser(Y_minimum*0.9 if logY else 0,(10.0 if logY else 1.1)*h_obs.GetMaximum())
                setup_hist(h_obs)
                h_obs.Draw(obs_draw_option+'SAME')
                if(logY):
                    plotpad.SetLogy()
                legend.AddEntry(h_obs,channel['obs'],'p')

                
                for ibackground in range(len(backgrounds)):
                    shapes_background[ibackground].SetLineColor(colors[backgrounds[ibackground]])
                    shapes_background[ibackground].SetLineWidth(2)                    
                    total_b.Add(shapes_background[ibackground],1)
                    total_sb.Add(shapes_background[ibackground],1)
                    shapes_background[ibackground].SetLineStyle(2)                    
                    shapes_background[ibackground].Draw(fit_draw_option+'SAME')
                    legend.AddEntry(shapes_background[ibackground],backgrounds[ibackground],'l')
                total_b.SetLineColor(46)
                total_b.SetLineWidth(2)
                total_b.Draw(fit_draw_option + 'SAME')
                legend.AddEntry(total_b,'Fit b','l')

                for isignal in range(len(signals)):
                    shapes_signal[isignal].SetLineColor(colors[signals[isignal]])
                    shapes_signal[isignal].SetLineWidth(2)
                    total_sb.Add(shapes_signal[isignal])
                    if(signal_scale>0):
                        shapes_signal[isignal].Scale(signal_scale)
                    shapes_signal[isignal].Draw(fit_draw_option+'SAME')
                    legend.AddEntry(shapes_signal[isignal],signals[isignal]+('*%.1f'%signal_scale if signal_scale>0 else ''),'l')
                total_sb.SetLineColor(32)
                total_sb.SetLineWidth(2)
                total_sb.Draw(fit_draw_option+'SAME')
                legend.AddEntry(total_sb,'Fit s+b','l')
                h_obs.Draw(obs_draw_option+'SAME')
                legend.Draw('SAME')

                normal_ratio = True
                ratiopad.cd()
                if(normal_ratio):
                    ratio_hist = h_obs.Clone()
                    ratio_hist.GetYaxis().SetTitle("#frac{obs.}{fit}")
                    ratio_hist.GetYaxis().SetRangeUser(0.5,1.5)
                    ratio_hist.Divide(total_sb)
                    ratio_hist.SetLineColor(1)
                    ratio_hist.SetLineWidth(2)
                    ratio_hist.SetMarkerColor(1)
                    ratio_hist.SetMarkerStyle(33)
                    setup_ratio_hist(ratio_hist)
                    ratio_hist.GetYaxis().SetRangeUser(0.3,1.7)

                    ratio_hist_1 = h_obs.Clone()
                    ratio_hist_1.GetYaxis().SetTitle("#frac{obs.}{total b}")
                    ratio_hist_1.Divide(total_b)
                    ratio_hist_1.SetLineColor(46)
                    ratio_hist_1.SetLineWidth(2)
                    ratio_hist_1.SetMarkerStyle(1)

                    ratio_hist.Draw('EP1')
                    # ratio_hist_1.Draw('EPSAME')

                else:
                    #constructing data-BG/sigma
                    ratio_hist = h_obs.Clone()
                    ratio_hist.GetYaxis().SetTitle("#frac{X - total B}{#sigma_{Data}}")
                    ratio_hist.Add(total_b,-1.)

                    obs_sigma = ratio_hist.Clone()
                    for i in range(obs_sigma.GetNbinsX()+1):
                        obs_sigma.SetBinContent(i,ratio_hist.GetBinError(i))
                        obs_sigma.SetBinError(i,1)
                    ratio_hist.Divide(obs_sigma)
                    
                    #constructing fit-BG/sigma
                    ratio_fit = total_sb.Clone()
                    ratio_fit.Add(total_b,-1)

                    fit_sigma = ratio_fit.Clone()
                    for i in range(fit_sigma.GetNbinsX()):
                        fit_sigma.SetBinContent(i,ratio_fit.GetBinError(i))
                    ratio_fit.Divide(obs_sigma)

                    signal_ratios = []
                    for signal_shape in shapes_signal:
                        signal_sigma = signal_shape.Clone()
                        for i in range(signal_shape.GetNbinsX()):
                            signal_sigma.SetBinContent(i,signal_shape.GetBinError(i))
                            signal_sigma.SetBinError(i,1)

                        new_signal_ratio = signal_shape.Clone()
                        new_signal_ratio.Divide(obs_sigma)
                        signal_ratios.append(new_signal_ratio)

                    ratio_hist.GetXaxis().SetTitleSize(xTitleSize)
                    setup_ratio_hist(ratio_hist)
                    # ratioHist.GetYaxis().SetRangeUser(-1,1.7)
                    ratio_hist.Draw('EP')
                    ratio_fit.Draw('HistSAME')
                    for signal_ratio in signal_ratios:
                        signal_ratio.Draw('HistSAME')
                ratioXMin = ratio_hist.GetXaxis().GetXmin()
                ratioXMax = ratio_hist.GetXaxis().GetXmax()

                zeropercent = ROOT.TLine(ratioXMin,1,ratioXMax,1)
                plus10percent = ROOT.TLine(ratioXMin,1.1,ratioXMax,1.1)
                minus10percent = ROOT.TLine(ratioXMin,0.9,ratioXMax,0.9)
                plus10percent.SetLineStyle(ROOT.kDashed)
                minus10percent.SetLineStyle(ROOT.kDashed)
                if(normal_ratio):
                    zeropercent.Draw()
                    plus10percent.Draw()
                    minus10percent.Draw()

                plotpad.cd()
                latex = ROOT.TLatex()
                latex.SetNDC(1)
                latex.SetTextSize(24)
                if('w' in channel_str.lower()):
                    pt_bin = tuple(channel['pt_bin'].split('to'))
                    latex.DrawLatex(0.20,0.95,"%s-region  - %s GeV #leq p_{T} < %s GeV"%(region,*pt_bin))


                c.RedrawAxis()
                # c.SaveAs(out_dir+'/'+hist_dir+('_logY' if logY else '')+'.pdf')
                c.SaveAs(out_dir+'/'+hist_dir+('_logY' if logY else '')+'.pdf')
