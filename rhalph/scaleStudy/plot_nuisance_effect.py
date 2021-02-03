from __future__ import print_function
import ROOT,json,os
import numpy as np
import sys
# sys.path.append('../python/')
# import cms_style
# cms_style.cms_style()
#selection = 'WJetsL4fb'
# selection = 'WJetsL4fbTF1x1TF1x1NormUnc0p05Bins5GeV'
# selection = 'WJets10GeVBins'
# selection = 'Top'               
# selection = 'WfromTop'
# prefix="/nfs/dust/cms/user/albrechs/JetMassCalibration/scaleStudy/WJets_CHFHighBin_scaleStudy/"

#prefix = "TopPassAndFail/"
#prefix = "TopFail/"
# workdir_varied = prefix+"workdir_scaleStudy_"+selection
# workdir_unvaried = prefix+"workdir_nominal_"+selection

# rhalphdir = os.getcwd()
rhalphdir = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/'

ROOT.gROOT.SetBatch(True)

def collect_and_plot_scaleStudy(nominal_diagnostics,up_diagnostics,down_diagnostics,outdir,categoryName,categoryConfig):
    if(not os.path.isdir(outdir)):
        os.makedirs(outdir)
        os.makedirs(outdir+'/workdir_nominal/%s'%categoryName)
        os.makedirs(outdir+'/workdir_variations/%s00allup'%categoryName)
        os.makedirs(outdir+'/workdir_variations/%s00alldown'%categoryName)
        os.system('cp '+nominal_diagnostics+' ' +outdir+'/workdir_nominal/%s/fitDiagnostics.root'%categoryName)
        os.system('cp '+up_diagnostics+' ' +outdir+'/workdir_variations/%s00allup/fitDiagnostics.root'%categoryName)
        os.system('cp '+down_diagnostics+' ' +outdir+'/workdir_variations/%s00alldown/fitDiagnostics.root'%categoryName)

    # workdir_unvaried+'/'+name.replace('_','')+'/fitDiagnostics.root'
    # model_name = (name+'_'+category+direction).replace('_','')

    # workdir_varied+'/'+model_name+'/fitDiagnostics.root'
    plot_nuisance_effect_grid(categoryName,categoryConfig,cycle_suffix='',workdir_prefix=outdir)
    
def plot_nuisance_effect_grid(name,config,cycle_suffix='',workdir_prefix="/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/scaleStudy/WJets_CHFLowBin_scaleStudy/"
):
    workdir_varied = workdir_prefix+"/workdir_variations"
    workdir_unvaried = workdir_prefix+"/workdir_nominal"
    print(workdir_unvaried)
    unvaried_nuisances = []
    unvaried_grid_file = ROOT.TFile(rhalphdir+"/../Histograms/grid_"+name+".root")
    h_grid = unvaried_grid_file.Get("grid")
    h_cat = unvaried_grid_file.Get("categories")
    categories = [h_cat.GetXaxis().GetBinLabel(i) for i in range(1,h_cat.GetXaxis().GetNbins()+1)]
    Npt = h_grid.GetXaxis().GetNbins()
    Neta = h_grid.GetYaxis().GetNbins()
    Ncat = len(categories)

    variation = 0.1
    parname_unvaried = []
    par_unvaried = []
    error_unvaried = []
    
    matrix={'up':np.zeros((Ncat,Ncat)),'down':np.zeros((Ncat,Ncat))}
    print('using',workdir_unvaried+cycle_suffix+'/'+name.replace('_','')+'/fitDiagnostics.root')
    if(not os.path.isfile(workdir_unvaried+cycle_suffix+'/'+name.replace('_','')+'/fitDiagnostics.root')):
        return matrix
    f_r_unvaried = ROOT.TFile(workdir_unvaried+cycle_suffix+'/'+name.replace('_','')+'/fitDiagnostics.root','READ')
    if('fit_s' not in f_r_unvaried.GetListOfKeys()):
        return matrix
    fitargs_unvaried = f_r_unvaried.Get('fit_s').floatParsFinal()
    
    for i in range(Npt):
        for j in range(Neta):
            for k in range(Ncat):
                parameter_name = "massScale_"
                parameter_name += "pt" + str(i)
                parameter_name += "_eta" + str(j)
                parameter_name += "_" + categories[k]
                params = fitargs_unvaried.find(parameter_name)
                central = params.getValV()    
                parname_unvaried.append(parameter_name)
                par_unvaried.append(central)
                error_unvaried.append(params.getError())                

    Npoints = len(parname_unvaried)

    h_unvaried = ROOT.TH1F('central values','central',Npoints,0,Npoints)
    
    for i in range(Npoints):
        h_unvaried.GetXaxis().SetBinLabel(i+1,parname_unvaried[i])
        h_unvaried.SetBinContent(i+1,par_unvaried[i])
        h_unvaried.SetBinError(i+1,error_unvaried[i])

    c_unvaried = ROOT.TCanvas("c_unvaried","c_unvaried",900,600)
    h_unvaried.GetYaxis().SetRangeUser(-1.5,2.5)
    h_unvaried.SetMarkerStyle(8)
    h_unvaried.SetLineColor(1)
    h_unvaried.Draw('PE0X0')

    zero  = ROOT.TLine(0,0,Npoints,0)
    zero.SetLineWidth(2)
    zero.SetLineStyle(ROOT.kDashed)
    zero.Draw("SAME")
    
    # outdir = 'nuisance_plots_'+selection+'_'+prefix.replace('_','')
    outdir = workdir_prefix+'/nuisance_plots'
    print(outdir)
    if(not os.path.isdir(outdir)):
        os.makedirs(outdir)

    c_unvaried.SaveAs(outdir+('/central_fit_'+name+cycle_suffix+'.pdf').replace('_',''))
    
    
    pars_varied={'up':{},'down':{}}
    for ptbin in range(len(config['ptbins'])-1):
        for etabin in range(len(config['etabins'])-1):
            for pfflavour in config['pfflavours']:
                for direction in ['up','down']:
                    # if 'oneScale' in name and direction == 'down':
                        # continue
                    variation_name = "massScale_"
                    variation_name += "pt" + str(ptbin)
                    variation_name += "_eta" + str(etabin)
                    variation_name += "_" + pfflavour+'_'+direction
                    category = '%i_%i_%s_'%(ptbin,etabin,pfflavour)                        
                    # variation_name = name+'_'+category+direction
                    model_name = (name+'_'+category+direction).replace('_','')

                    parname = []
                    par = []
                    if(not os.path.isfile(workdir_varied+cycle_suffix+'/'+model_name+'/fitDiagnostics.root')):
                        print(workdir_varied+cycle_suffix+'/'+model_name+'/fitDiagnostics.root does not exist!')
                        pars_varied[direction].update({variation_name:[-10.0]*Npt*Neta*Ncat})
                        continue
                    f_r = ROOT.TFile(workdir_varied+cycle_suffix+'/'+model_name+'/fitDiagnostics.root','READ')

                    if('fit_s' not in f_r.GetListOfKeys()):
                        pars_varied[direction].update({variation_name:[-10.0]*Npt*Neta*Ncat})
                        continue
                    
                    try:
                        fitargs = f_r.Get('fit_s').floatParsFinal()
                    except:
                        pars_varied[direction].update({variation_name:[-10.0]*Npt*Neta*Ncat})
                        continue

                    for i in range(Npt):
                        for j in range(Neta):
                            for k in range(Ncat):
                                parameter_name = "massScale_"
                                parameter_name += "pt" + str(i)
                                parameter_name += "_eta" + str(j)
                                parameter_name += "_" + categories[k]
                                params = fitargs.find(parameter_name)
                                central = params.getValV()    
                                parname.append(parameter_name)
                                par.append(central)
                    print(variation_name)
                    pars_varied[direction].update({variation_name:par})

    # print(parname_unvaried,par_unvaried)

    # print(pars_varied)
    # n_massscales=len(pars_varied['up'].keys())
    n_massscales=len(par_unvaried)
    max_val = 0
    print(n_massscales)
    comp_plot = ROOT.TH2F("nuisance_comp_"+name,"nuisance_comp_"+name,n_massscales,0,n_massscales,2*n_massscales+1,0,2*n_massscales+1)
    for i in range(n_massscales):
        comp_plot.GetXaxis().SetBinLabel(i+1,parname_unvaried[i].replace('massScale_',''))
        comp_plot.GetYaxis().SetBinLabel(i+1,parname_unvaried[i].replace('massScale_','')+'_down')
        comp_plot.GetYaxis().SetBinLabel(i+1+n_massscales+1,parname_unvaried[i].replace('massScale_','')+'_up')
        print(parname_unvaried[i].replace('massScale_',''))
        for j in range(len(par_unvaried)):
            # if 'oneScale' not in name:
            print(parname_unvaried[j])
            print('central',par_unvaried[j])
            print('down',pars_varied['down'][parname_unvaried[i]+'_down'][j])
            print('down',pars_varied['up'][parname_unvaried[i]+'_up'][j])
            down_value = -par_unvaried[j]+pars_varied['down'][parname_unvaried[i]+'_down'][j]
            comp_plot.SetBinContent((j+1),i+1,down_value)
            max_val = max(abs(down_value),max_val)
            matrix['down'][j][i] = down_value
            # else:
                # print(par_unvaried[j],pars_varied['up'][parname_unvaried[i]+'_up'][j])
            up_value = -par_unvaried[j]+pars_varied['up'][parname_unvaried[i]+'_up'][j]
            comp_plot.SetBinContent((j+1),i+1+n_massscales+1,up_value)
            max_val = max(abs(up_value),max_val)
            matrix['up'][j][i] = up_value
    
            
    c = ROOT.TCanvas('nuisance_comp_'+name,'nuisance_comp_'+name,900,600)
    
    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPalette(ROOT.kBlackBody)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    ROOT.gStyle.SetNumberContours(250)

    # Red    = np.array([ 1.00, 1.0 ,0.00])
    # Green  = np.array([ 0.00, 1.0 ,0.00])
    # Blue   = np.array([ 0.00, 0.8 ,1.00])
    # Length = np.array([ 0.00, 0.5,1.00])
    # nb=100
    # ROOT.TColor.CreateGradientColorTable(3,Length,Red,Green,Blue,nb,1.0)
    # comp_plot.SetContour(nb)


    max_val*=1.1
    comp_plot.GetZaxis().SetRangeUser(-max_val,max_val)
    comp_plot.GetZaxis().SetTitle("value from varied fit - value from unvaried fit ")
    ROOT.gPad.SetLeftMargin(.2)
    ROOT.gPad.SetRightMargin(.2)
    ROOT.gPad.SetBottomMargin(.15)


    comp_plot.Draw('textcolz')

    c.SaveAs(outdir+('/nuisance_comp_'+name+cycle_suffix+'.pdf').replace('_',''))
    # c.SaveAs(outdir+('/nuisance_comp_'+name+cycle_suffix+'.C').replace('_',''))

    print('up variation:',matrix['up'])
    print('-> determinant:',np.linalg.det(matrix['up']))
    print('down variation:',matrix['down'])
    print('-> determinant:',np.linalg.det(matrix['down']))

    bad_example = np.identity(4)
    print('-> determinant:',np.linalg.det(bad_example))
    return matrix

if(__name__ == "__main__"):
    grids={
        "oneScale":{
            "ptbins":[0.,100000.],
            "etabins":[0.,9.],
            "pfflavours":["all"]
        },
        # "PF_flavours":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # },
        # "2ptbins":{
        #     "ptbins":[0.,10.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["all"]
        # },
        # "2etabins":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0., 1.479, 9.],
        #     "pfflavours":["all"]
        # },
        # "PF_flavours_2ptbins":{
        #     "ptbins":[0.,10.,100000.],
        #     "etabins":[0., 9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # },
        # "PF_flavours_2etabins":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0., 1.479, 9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # }
    }
    # config_template = json.load(open('WJets_scaleStudy.json'))
    # submit_variations(grids,config_template)

    # pt_bins = ['_inclusive','_300To500','_300To400','_400To500','_500ToInf',"_300To350","_350To400","_400To450","_450To500","_500To550","_550To600","_600ToInf"]
    #cycles = ["_300To500","_300To350","_350To400","_400To450","_450To500","_500To550","_550To600","_600ToInf"]
    cycles = [""]
    # cycles = ['TFPt%iRho%iTFPt%iRho%i'%(qcd_pt_order_MCTF,qcd_rho_order_MCTF,qcd_pt_order,qcd_rho_order) for qcd_pt_order_MCTF in range(1,6) for qcd_rho_order_MCTF in range(1,6) for qcd_pt_order in range(1,6) for qcd_rho_order in range(1,6)]
    print(cycles)
    # if(len(cycles)>0):
    #     matrices = {}
    #     for suffix in cycles:
    #         for name,config in grids.items():
    #             matrices.update({name+suffix:plot_nuisance_effect_grid(name, config, suffix)})

    #     h_det_up = ROOT.TH1F('determinants_up','determinants_up',len(matrices.keys()),0,len(matrices.keys()))
    #     h_det_down = ROOT.TH1F('determinants_down','determinants_down',len(matrices.keys()),0,len(matrices.keys()))

    #     for i,pt_bin in enumerate(cycles):
    #         for name,config in grids.items():
    #             h_det_up.GetXaxis().SetBinLabel(i+1,pt_bin.replace('_',''))
    #             h_det_up.GetYaxis().SetTitle('determinant')
    #             h_det_up.GetYaxis().SetRangeUser(-1.5,2.5)
        
    #             # h_det_up.SetBinContent(i+1,abs(1-np.linalg.det(matrices[name+pt_bin]['up'])))
    #             # h_det_down.SetBinContent(i+1,abs(1-np.linalg.det(matrices[name+pt_bin]['down'])))
    #             h_det_up.SetBinContent(i+1,(np.linalg.det(matrices[name+pt_bin]['up'])))
    #             h_det_down.SetBinContent(i+1,(np.linalg.det(matrices[name+pt_bin]['down'])))
                
    #             h_det_up.SetLineColor(ROOT.kBlue+2)
    #             h_det_down.SetLineColor(ROOT.kRed+2)

    #     c = ROOT.TCanvas('canv_determinant','canv_determinant',600,600)
    #     leg = ROOT.TLegend(0.1,0.7,0.9,0.9)

    #     leg.AddEntry(h_det_up,'determinant up-variations','l')
    #     leg.AddEntry(h_det_down,'determinant down-variations','l')

    #     h_det_up.Draw('H')
    #     h_det_down.Draw('HSAME')

    #     leg.Draw('SAME')
    
    #     c.SaveAs('nuisance_plots_'+selection+'_'+prefix.replace('_','')+'/determinants.pdf')
    # else:
    # for name,config in grids.items():
    #     plot_nuisance_effect_grid(name, config)

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--nominal',type=str)
    parser.add_argument('--up',type=str)
    parser.add_argument('--down',type=str)
    parser.add_argument('--out',type=str)
    parser.add_argument('--catName',type=str)
    args = parser.parse_args()
    collect_and_plot_scaleStudy(args.nominal,args.up,args.down,args.out,args.catName,grids[args.catName])
