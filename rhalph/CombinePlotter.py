#!/usr/bin/env python
from __future__ import print_function
import os,ROOT,argparse
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

import sys
sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/python')
# import fitplotter
import cms_style
cms_style.extra_text="Preliminary Simulation"
cms_style.cms_style()
cms_logo = False

class CombinePlotter:
    def __init__(self):
        # print('Nothing to do here. This class is just a wrapper for some worflow methods using combine')
        self.methods = [func for func in dir(CombinePlotter) if callable(getattr(CombinePlotter, func))]        
        self.methods.remove('__init__')
        

    def plot_gof_result(self,gof_base_file="",toys_result_file="",algo="saturated",lumi=41.8,bernstein_order="2x2"):
        algo_x_label = {
            "saturated":"#Pi_{i} exp(-(d_{i}-f_{i})^{2}/(2#sigma_{i}^{2}))",
            "KS":"",
            "AD":""
        }

        observed_file = ROOT.TFile(gof_base_file,"READ")
        toys_file = ROOT.TFile(toys_result_file,"READ")
        
        
        observed = [entry.limit for entry in observed_file.Get("limit")][0]
        toys = [entry.limit for entry in toys_file.Get("limit")]
        
        minimum_test_statistic = 0.8*min(observed,min(toys))
        maximum_test_statistic = 1.2*max(observed,max(toys))
        
        hist_title = "GOF - BernsteinOrders: %s"%(bernstein_order)
        toys_hist = ROOT.TH1F(hist_title ,hist_title,70,minimum_test_statistic,maximum_test_statistic)
        for val in toys:
            toys_hist.Fill(val,1)
            
            
        toys_hist_norm = toys_hist.Clone()
        test_stat_norm = toys_hist_norm.Integral()
        if test_stat_norm <= 0:
            raise ValueError("Integral of test statistic from toys is <=0. This is very wrong...")
        toys_hist_norm.Scale(1./test_stat_norm)
    
        canv = ROOT.TCanvas("GoodnessOfFit","GoodnessOfFit",600,600)
        legend = ROOT.TLegend(0.6,0.6,0.9,0.9)
        canv.SetLeftMargin(0.15)
        toys_hist.SetMarkerStyle(8)
        toys_hist.SetMarkerColor(1)
        toys_hist.SetLineColor(1)
        toys_hist.Draw("PE1")
        toys_hist.GetYaxis().SetTitle("toys")
        toys_hist.GetXaxis().SetTitle("-2 ln #lambda_{%s}"%algo)
        toys_hist.GetXaxis().SetTitleSize(0.045)
        toys_hist.GetYaxis().SetRangeUser(0,toys_hist.GetMaximum()*1.5)
        legend.AddEntry(toys_hist,"toys (N=%i)"%len(toys),"pe1")
            
        observed_arrow = ROOT.TArrow(observed,0,observed,toys_hist.GetMaximum()/2,0.05,"<|");
        observed_arrow.SetAngle(40)
            
        observed_arrow.SetLineColor(46)
        observed_arrow.SetFillColor(46)
        observed_arrow.SetLineWidth(2)
        
        observed_arrow.Draw()
        legend.AddEntry(observed_arrow,"observed = %.2f"%observed,"l")
        legend.AddEntry(None,"p = %.2f"%toys_hist_norm.Integral(toys_hist_norm.FindBin(observed),toys_hist_norm.GetNbinsX()),"")
            
        legend.Draw("SAME")
        cms_style.extra_text_rel_X=0.12
        cms_style.draw_lumi(canv, lumi,do_extra_text=cms_logo,out_of_frame = True,do_cms_text=cms_logo,private_work=False)
            
        canv.SaveAs("goodness_of_fit.pdf")

    def plot_fit_diagnostics_results(self,model_name=""):
        import json
        import fitplotter
        print("plotting ",model_name)
        config = json.load(open(model_name+'/config.json'))

        qcd_only_fit = (not os.path.isfile(model_name+"/fit_shapes.root") and os.path.isfile(model_name+"/qcdmodel/fit_shapes.root") )        
        if(qcd_only_fit):
            config['ModelName'] += '/qcdmodel'
            for ch_name,ch_config in config['channels'].items():
                ch_config['samples'] = ["QCD"]
                ch_config['signal'] = []
        fitplotter.plot_fit_result(config,plot_total_sig_bkg = False,use_config_samples=True)
        fitplotter.plot_fit_result(config,logY=True,plot_total_sig_bkg = False,use_config_samples=True)

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--method',type=str,choices=globals()['CombinePlotter']().methods,required=True)
    parser.add_argument('--parameter',type=str,help="pass parameters to plotting function. delimiter: ';'")
    args = parser.parse_args()

    method = getattr(globals()['CombinePlotter'](),args.method)
    print(method.__defaults__)
    parameter_types = [type(a) for a in method.__defaults__]
    print(parameter_types)
    print(args.parameter.split(';'))
    parameters = [parameter_types[i](args.parameter.split(";")[i]) for i in range(len(args.parameter.split(';')))]
    method(*tuple(parameters))


