#!/nfs/dust/cms/user/albrechs/miniconda3/envs/coffea/bin/python
# from __future__ import print_function
import numpy as np 
import sys,os
import json
# sys.path.append(os.getcwd()+'/rhalphalib/')
# import ROOT as _ROOT
#import uproot4 as uproot
import uproot
# import rhalphalib as rl
# rl.util.install_roofit_helpers()

# _ROOT.gROOT.SetBatch(True)
# def _RooFitResult_nameArray(self):
#     return np.array([p.GetName() for p in self.floatParsFinal()])

# _ROOT.RooFitResult.nameArray = _RooFitResult_nameArray

# def _RooFitResult_massScales(self):
#     result = []
#     for p in self.floatParsFinal():
#         if("massScale" in p.GetName()):
#              result.append([p.GetName(),p.getVal(),p.getErrorHi(),p.getErrorLo()])
#     return np.array(result,dtype=object)

# _ROOT.RooFitResult.massScales = _RooFitResult_massScales

# def _RooFitResult_filterParameters(self,substr=""):
#     result = []
#     substr_list = substr if isinstance(substr, list) else [substr]
#     for p in self.floatParsFinal():
#         if(all([s.lower() in p.GetName().lower() for s in substr_list])):
#              result.append([p.GetName(),p.getVal(),p.getErrorHi(),p.getErrorLo()])
#     return np.array(result,dtype=object)

# _ROOT.RooFitResult.filterParameters = _RooFitResult_filterParameters


import glob,os,sys
# sys.path.append("../python/")
# import cms_style
# cms_style.cms_style()

def iterative_gaus(hist, peak_to_extract = 'W',outname='test.pdf',boundaries=(50,120,80)):
    print('h.name',hist.GetName(),'peak',peak_to_extract)
    fitOpt = "RSQ"
    
    h_norm = hist.Clone()
    norm = h_norm.Integral()
    h_norm.Scale(1./norm if norm>0 else 1.0)

    c = _ROOT.TCanvas("c","c",600,600)
    cms_style.ratio_plot = False
    plotpad,_,_ = cms_style.setup_pads(c)
    plotpad.cd()
    cms_style.setup_hist(h_norm)
    h_norm.Draw("HIST")

    legend = _ROOT.TLegend(0.60,0.8,0.9,0.9)
    print(boundaries)
    gaus1_tf1 = _ROOT.TF1('firstGaus','gaus(0)',*boundaries[:-1])
    gaus1_tf1.SetParameter(0,1.)
    gaus1_tf1.SetParameter(1,boundaries[2])
    # gaus1_tf1.SetParameter(2,h_norm.GetRMS())
    gaus1_tf1.SetLineColor(32)

    fit1 = h_norm.Fit('firstGaus',fitOpt)
    gaus1_norm = gaus1_tf1.GetParameter(0)
    gaus1_mean = gaus1_tf1.GetParameter(1)
    gaus1_sigma = gaus1_tf1.GetParameter(2)
    gaus2_tf1 = _ROOT.TF1('my-gaus2','gaus(0)',gaus1_mean-1.5*gaus1_sigma,gaus1_mean+1.5*gaus1_sigma)
    gaus2_tf1.SetParameter(0,gaus1_norm)
    gaus2_tf1.SetParameter(1,gaus1_mean)
    gaus2_tf1.SetParameter(2,gaus1_sigma)
    gaus2_tf1.SetLineColor(46)

    fit2 = h_norm.Fit('my-gaus2',fitOpt)

    fit1_status,fit2_status = None,None
    try:
        fit1_status = fit1.Status()
    except:
        print("fit 1 failed very much!")
    try:
        fit2_status = fit2.Status()
    except:
        print("fit 2 failed very much!")
    
    print("Fit status:")
    print("first fit:\t",fit1_status)
    print("second fit:\t",fit2_status)    
    
    
    if(fit1_status is not None):
        gaus1_tf1.Draw("SAME")
    if(fit2_status is not None):
        gaus2_tf1.Draw("SAME")

    legend.AddEntry(gaus1_tf1,"first fit","l")
    legend.AddEntry(gaus2_tf1,"second fit","l")
    legend.Draw("SAME")
    
    c.SaveAs(outname)
    try:
        if(fit2_status == 0):
            fitresults = gaus2_tf1.GetParameter(1)
        elif(fit1_status == 0):
            fitresults = gaus1_tf1.GetParameter(1)
        else:            
            fitresults =  -1.0
    except:
        print("fit failed. setting mass manually")
        fitresults = -1.0
        # fitresults = h_norm.GetMean()
    if(fitresults == boundaries[2]):
        fitresults = -1.0
    return fitresults


    return peak_positions

class FitResults:
    def __init__(self,config):
        self._config = config
        self.name = self._config['ModelName']
        # diagnostics_file_path = self.name +'/fitDiagnostics.root'
        # if(os.path.isfile(diagnostics_file_path)):
        #     f_ = _ROOT.TFile(diagnostics_file_path,'READ')    
        #     self.fit_result = f_.Get('fit_s')
        # elif(os.path.isfile(diagnostics_file_path.replace("fitDiagnostics","fitDiagnosticsTest"))):
        #     diagnostics_file_path = diagnostics_file_path.replace("fitDiagnostics","fitDiagnosticsTest")
        #     f_ = _ROOT.TFile(diagnostics_file_path,'READ')    
        #     self.fit_result = f_.Get('fit_s')
        # else:
        #     self.fit_result = None
        # fit_result_file = np.load(f"{self.name}/{self.name}fitResult.npy",allow_pickle=True,encoding='bytes')
        # self.fit_result = fit_result_file.item()
        self.fit_result = json.load(open(f"{self.name}/{self.name}fitResult.json","r"))
        
        shapes_file_path = self.name +'/fit_shapes.root'
        if(os.path.isfile(shapes_file_path)):
            # self.fit_shapes = _ROOT.TFile(shapes_file_path,'READ')
            self.fit_shapes_uproot = uproot.open(shapes_file_path)
        else:
            self.fit_shapes = None
        self.results = {}
        # self.hists=[]

    def get_tagging_sf(self, prefix = "top"):
        pass
        # pt_edges = np.array(self._config.get('pt_edges',[200,300,400,500,550,600,10000]),dtype=float)

        # sf_infos = ["nominal","up","down"]
        # hists = {d:_ROOT.TH1F("{prefix}_tagging_sf_{d}".format(prefix=prefix,d=d),"{prefix}_tagging_sf_{d}".format(prefix=prefix,d=d),len(pt_edges)-1,pt_edges) for d in sf_infos}

        # if(self.fit_result is None):
        #     return {}
        # tagging_sfs = {p[0]:list(p[1:]) for p in self.fit_result.filterParameters(prefix+'_tag_eff_sf')}
        # for ipt,pt in enumerate(pt_edges[:-1]):
        #     this_par = None
        #     for parname in tagging_sfs.keys():
        #         if(str(int(pt)) in parname):
        #             this_par = parname
        #             break
        #     for i,n in enumerate(sf_infos):
        #         if(i==0):
        #             hists[n].SetBinContent(ipt+1,tagging_sfs[this_par][i])
        #         else:
        #             hists[n].SetBinContent(ipt+1,tagging_sfs[this_par][0]+tagging_sfs[this_par][i])

        # self.results.update({prefix+"tagging_sf":tagging_sfs})
        # for h in hists.values():
        #     self.hists.append(h)
        # return tagging_sfs


    # parameters = self.fit_result.filterParameters(prefix+'_tag_eff_sf')
        # print(parameters)
    
    # def get_massScales(self):    
    #     massScales = {}
    #     print("getting massscales")
    #     if(self.fit_result is None):
    #         print("No fit result. skipping")
    #         return {}
    #     massScales_obj = self.fit_result.massScales()
    #     # print("parameters",self.fit_result.nameArray())
    #     # print(massScales_obj)
    #     for param in massScales_obj:
    #         #print(param)
    #         pt_edges = None
    #         for channel_name,channel in self._config['channels'].items():
    #             pt_bin = channel_name.lower()
    #             pt_bin = pt_bin.split("pt")[-1]
    #             if(pt_bin in param[0]):
    #                 pt_edges = channel['pt_bin'].split("to")
    #         factor = 0.5 if 'Unscaled' in self.name else 1.0
    #         central = (100+param[1]*factor)/100
    #         error_up = (param[2]*factor)/100
    #         error_down = (param[3]*factor)/100
    #         massScales.update({param[0]:{'edges':pt_edges,'vals':[central,error_up,error_down]}})
    #     self.results.update({"jms":massScales})
    #     return massScales
    def get_massScales(self):    
        massScales = {}
        print("getting massscales")
        if(self.fit_result is None):
            print("No fit result. skipping")
            return {}
        #print(self.fit_result)
        massScales_parameters = [[name,*self.fit_result[name]] for name in self.fit_result.keys() if 'massScale' in name]
        print(massScales_parameters)
        # print("parameters",self.fit_result.nameArray())
        # print(massScales_obj)
        for param in massScales_parameters:
            #print(param)
            pt_edges = None
            for channel_name,channel in self._config['channels'].items():
                pt_bin = channel_name.lower()
                pt_bin = pt_bin.split("pt")[-1]
                if(pt_bin in param[0]):
                    pt_edges = channel['pt_bin'].split("to")
            factor = 0.5 if 'Unscaled' in self.name else 1.0
            central = (100+param[1]*factor)/100
            error_up = (param[2]*factor)/100
            error_down = (param[3]*factor)/100
            massScales.update({param[0]:{'edges':pt_edges,'vals':[central,error_up,error_down]}})
        self.results.update({"jms":massScales})
        return massScales

    def get_peak_positions_new(self):
        selections = set([self._config['channels'][c]['selection'] for c in self._config['channels']])
        fit_setup_info = {
            'top':{'samples':['TTToSemiLeptonic_mergedTop','TTToSemiLeptonic_mergedW']},
            'W':  {'samples':['WMatched','ZMatched']},
            'Zbb':{'samples':['WMatched','ZMatched']},
        }

        import numpy as np
        from scipy.optimize import curve_fit


        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        import mplhep as hep
        plt.style.use(hep.style.CMS)

        font_size = 20
        mpl.rcParams['axes.labelsize'] = font_size
        mpl.rcParams['axes.labelsize'] = font_size

        def iterative_gaussian_fit(th1,name='',outDir='./',title=''):
            hist, bin_edges = th1.to_numpy()
            hist = hist/hist.sum()
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
                
            mean = np.average(bin_centres, weights=hist)
    
            def gauss(x, *p):
                A, mu, sigma = p
                return A*np.exp(-(x-mu)**2/(2.*sigma**2))

            p0 = [1., mean, 10.]
            
            coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
            #sometimes curve_fit returns sigma as negative
            coeff[2]=abs(coeff[2])
            
            xmin = coeff[1]-2.5*coeff[2]
            xmax = coeff[1]+2.5*coeff[2]
            bin_min = np.searchsorted(bin_centres , xmin, side="left")
            bin_max = np.searchsorted(bin_centres , xmax, side="right")

            if(bin_max-bin_min < 3): # if number of bins is smaller than 3 increase by 1 in both directions (curve_fit will run into problems otherwise)
                bin_min -= 1
                bin_max += 1

            coeff1, var_matrix1 = curve_fit(gauss, bin_centres[bin_min:bin_max], hist[bin_min:bin_max], p0=coeff)
            #sometimes curve_fit returns sigma as negative
            coeff1[2]=abs(coeff1[2])
            
            #getting drawable arrays from fits
            msd = np.linspace(bin_edges[0],bin_edges[-1],1000)
            hist_fit = gauss(msd, *coeff)
            msd1 = np.linspace(bin_edges[bin_min],bin_edges[bin_max],1000)
            hist_fit1 = gauss(msd1,*coeff1)
            
            #setting up canvas
            fig,ax = plt.subplots(figsize = (10,8))
            hep.cms.text(",work in progress",ax=ax, fontsize=font_size)
            # ax.plot([],[], label=name,color='white')
            for subtitle in title.split("//"):
                ax.plot([],[], label="$\mathbf{%s}$"%subtitle,color='white')
            
            #plotting hist
            hep.histplot((hist,bin_edges),ax=ax)
            #plotting first fit
            ax.plot(msd, hist_fit, label='first fit: $\mu = %.2f, \sigma=%.2f$'%(coeff[1],coeff[2]))
            #plotting second fit
            ax.plot(msd1, hist_fit1, label='second fit: $\mu = %.2f, \sigma=%.2f$'%(coeff1[1],coeff1[2]))
            
            plt.xticks(fontsize = font_size-2)
            plt.yticks(fontsize = font_size-2)
            ax.set_xlabel("$m_{SD}$ [GeV]")
            ax.set_ylabel("$\Delta N / N$")
            #ax.legend(fontsize = font_size-2,loc='upper left',bbox_to_anchor=(1, 1))
            ax.legend(fontsize = font_size-2,loc='upper right')#,bbox_to_anchor=(1, 1))
            
            # print('Fitted mean = ', coeff[1],coeff1[1])
            # print('Fitted standard deviation = ', coeff[2],coeff1[2])
            if(not os.path.isdir(outDir)):
                os.makedirs(outDir)

            outname = outDir+name+".pdf"
            plt.savefig(outname)
            return coeff1[1],coeff1[2]

        
        peak_positions = {}
        for ch in self._config['channels']:
            selection = self._config['channels'][ch]['selection']
            regions = self._config['channels'][ch]['regions']
            regions.remove("fail")
            peak_positions.update({ch:{}})
            for region in regions:
                peak_positions[ch].update({region:{}})
                for sample in self._config['channels'][ch]['signal']:
                    peak_positions[ch][region].update({sample:{}})
                    for fit in ['prefit','postfit']:
                        hist_dir = ch+region+"_"+fit+"/"+sample
                        th1 = self.fit_shapes_uproot[hist_dir]
                        name = self.name + "_" +selection +"_"+ hist_dir.replace("/","_")
                        title = selection + "~" + region + "~" + sample.replace("_","~") + "//" +ch+" "+fit
                        outDir = 'iterative_gaussians/'
                        mu,sigma = iterative_gaussian_fit(th1,name,outDir,title)
                        pt_edges = [float(pt) for pt in self._config['channels'][ch]['pt_bin'].split("to")]
                        peak_positions[ch][region][sample].update({fit:{'values':[mu,sigma],'pt_edges':pt_edges}})
        self.results.update({'peak_positions':peak_positions})
        

    # def get_peak_positions(self, peak_to_extract = 'W',method="fit"):
    #     peak_positions={}
    #     if(self.fit_shapes is None):
    #         return peak_positions
    #     #print(self.fit_shapes)
    #     for ch_name,ch in self._config['channels'].items():
    #         if(ch['selection'] == 'W'):
    #             peak_to_extract_ = 'W'
    #         else:
    #             peak_to_extract_ = peak_to_extract

    #         regions = ch['regions']
    #         if('fail' in regions):
    #             regions.remove('fail')
    #         for region in regions:
    #             region_name = ch_name+region
    #             peak_positions.update({region_name:{}})
    #             for peak in (['W','top'] if peak_to_extract_=='both'  else [peak_to_extract_]):
    #                 peak_process = None
    #                 for p in ch['signal']:
    #                     if peak.lower() in p.lower():
    #                         peak_process = p
    #                         break
    #                 processes = [peak_process,'data_obs']
    #                 peak_positions[region_name].update({peak_process:{},'data_obs_'+peak_process:{}})
    #                 for process in processes:
    #                     for fit_suffix in ['prefit','postfit']:
    #                         hist_dir = region_name+'_'+fit_suffix+'/%s'
    #                         shape = self.fit_shapes.Get(hist_dir%process)
    #                         shape.SetName(region_name + " " + process+ " " +fit_suffix)
    #                         shape.SetTitle(region_name + " " + process+ " " +fit_suffix)
    #                         shape.GetXaxis().SetTitle("m_{SD} [GeV]")
    #                         shape.GetYaxis().SetTitle("#Delta N/ N")
    #                         if('data' in process):
    #                             total_procs = self.fit_shapes.Get(hist_dir%'TotalProcs')
    #                             total_procs.Add(self.fit_shapes.Get(hist_dir%peak_process),-1.0)
    #                             shape.Add(total_procs,-1.0)

    #                         outdir = config["ModelName"]+'/peak_fits/'
    #                         if(not os.path.isdir(outdir)):
    #                             os.makedirs(outdir)
                            
    #                         boundaries = {'W':(50,100 if ch["selection"]=="W" else 120,80),'top':(50,220,173)}
    #                         if(method == 'fit'):
    #                             mean = iterative_gaus(shape, peak,outdir + '/'+region_name+process+fit_suffix+'.pdf',boundaries = boundaries[peak])
    #                         else:
    #                             mean = shape.GetMean()
    #                         # print(mean)
    #                         peak_positions[region_name][process+("_"+peak_process if "data" in process else "")].update({fit_suffix:mean})
    #                     # c.SaveAs()
    #                     # del c
    #     self.results.update({"peak_positions":peak_positions})
        
    # def write_to_root_file(self,filename="fitResults.root"):
    #     f_ = _ROOT.TFile(filename,"UPDATE")
    #     f_.mkdir(self.name)
    #     f_.cd(self.name)
    #     for h in self.hists:
    #         h.Write()
    #     f_.Close()

    def write_to_root_file(self,filename="fitResults.root"):
        pass
    
if __name__ == '__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--fits',nargs='+', default=[],help='List of fits (separated by spaces).',required=True)
    parser.add_argument('--fit-peaks',action="store_true")
    parser.add_argument('-o','--output',default="fitResults.json")
    # configs = glob.glob("*/config.json")
    args = parser.parse_args()
    configs = [( fit_dir if "config.json" else f'{fit_dir}/config.json' ) for fit_dir in args.fits]
    
    all_results = {}
    # rf_name = "fitResults.root"
    # rf = _ROOT.TFile(rf_name,"RECREATE")
    # rf.Close()
    for config_path in configs:
        config = json.load(open(config_path))
        name =config["ModelName"]
        print(name)
        FR = FitResults(config)
        sf = FR.get_massScales()
        # if(config.get('TTbarTaggingEff','False')=='True'):
        #     top_tag_sf = FR.get_tagging_sf("top")
        #     W_tag_sf = FR.get_tagging_sf("W")
            
        #JMS.update({name:{'jms':sf,'peak_positions':None}})
        
        #FR.get_peak_positions(peak_to_extract='both',method='mean')

        if args.fit_peaks:
            FR.get_peak_positions_new()
        
        all_results.update({name:FR.results})

        # FR.write_to_root_file(filename=rf_name)
        
    open(args.output,"w").write(json.dumps(all_results,indent=2))
