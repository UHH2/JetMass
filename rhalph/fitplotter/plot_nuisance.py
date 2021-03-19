import ROOT
from ROOT import gStyle,gROOT
import sys,os
sys.path.append(os.getcwd()+'/rhalphalib/')
import rhalphalib as rl
rl.util.install_roofit_helpers()
import numpy as np

def _RooFitResult_nameArray(self):
    return np.array([p.GetName() for p in self.floatParsFinal()])

ROOT.RooFitResult.nameArray = _RooFitResult_nameArray

def _RooFitResult_massScales(self):
    result = []
    for p in self.floatParsFinal():
        if("massScale" in p.GetName()):
             result.append([p.GetName(),p.getVal(),p.getErrorHi(),p.getErrorLo()])
    return np.array(result,dtype=object)

ROOT.RooFitResult.massScales = _RooFitResult_massScales

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
gStyle.SetTitleSize(0.02, "XYZ")
gStyle.SetLabelSize(0.04, "XYZ")
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

def update_nuisance(name,fitargs,nuisances):
    nuisance = fitargs.find(name)
    effect = nuisances[name]['effect']
    nuisances[name]['central'] = 1 + nuisance.getVal() * effect
    nuisances[name]['down'] = 1 + (nuisance.getVal() - (nuisance.getAsymErrorLo() if nuisance.getAsymErrorLo() != 0 else nuisance.getErrorLo()) ) * effect
    nuisances[name]['up'] = 1 + (nuisance.getVal() + (nuisance.getAsymErrorHi() if nuisance.getAsymErrorHi() != 0 else nuisance.getErrorHi()) ) * effect
    
def plot_nuisances(config):

    f_r = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    def dummy(effect=0.01):
        return {'effect':effect,'central':1,'up':1.1,'down':0.9}
        
    nuisances = {'wqq_norm':dummy(0.00089),'zqq_norm':dummy(0.00072),'CMS_lumi':dummy(0.027)}

    fitargs = f_r.Get('fit_s').floatParsFinal()

    for nu in nuisances.keys():
        update_nuisance(nu,fitargs,nuisances)
    print(nuisances)

def print_nuisance(config):
    f_r = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    print(dir(f_r.Get('fit_s')))
    print(type(f_r.Get('fit_s').covQual()))
    c = ROOT.TCanvas("correlationHist","correlationHist",1000,1000)

    f_r.Get('fit_s').correlationHist().Draw("colz2")

    c.SaveAs('correlationHist.png')
    c.SaveAs('correlationHist.pdf')

def plot_mass_scale_nuisances(config):
    f_ = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    f_r_ = f_.Get('fit_s')
    ms_par = f_r_.massScales()

    variation = 0.01
    central = 1.+ms_par[:,1]*variation
    err_up = abs(ms_par[:,2])*variation
    err_down =  abs(ms_par[:,3])*variation

    Npoints = len(ms_par)
    g = ROOT.TGraphAsymmErrors(Npoints)
    for i in range(Npoints):
        g.SetPoint(i, i+0.5, central[i])        
        g.SetPointError(i, 0.0, 0.0, err_down[i], err_up[i])

    frame = ROOT.TH1F("frame", " ", Npoints, 0, Npoints)
    frame.GetYaxis().SetRangeUser(0.98, 1.02)
    for i in range(Npoints):
        frame.GetXaxis().SetBinLabel(i+1, ms_par[i,0])

    width = 200 + Npoints*75
    c = ROOT.TCanvas("c", "c", 2000, 600)
    ROOT.gPad.SetBottomMargin(.2)
    ROOT.gPad.SetLeftMargin(.08)
    frame.GetYaxis().SetTitle("JMS-SF")
    frame.GetYaxis().CenterTitle()

    frame.SetFillColor(ROOT.kWhite)
    frame.SetLineColor(ROOT.kWhite)
    frame.Draw()
    line = ROOT.TLine(0,1.0,Npoints,1.0)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineStyle(ROOT.kDashed)
    line.Draw("SAME")
    g.SetMarkerColor(ROOT.kBlack)
    g.SetMarkerStyle(8)
    g.SetMarkerSize(1)
    g.Draw("P SAME")
    c.SaveAs(config['ModelName']+'/plots/JetMassScale_Nuisance.pdf')



    
if __name__ == '__main__':
    import json
    config={'ModelName':'/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/TTbarSelectionNoWTopSplit'}
    a=plot_mass_scale_nuisances(config)    
    # plot_nuisances(config)
    #print_nuisance(config)
