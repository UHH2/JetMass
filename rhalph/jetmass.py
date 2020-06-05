import sys,os
sys.path.append(os.getcwd()+'/rhalphalib/')
import rhalphalib as rl
import numpy as np
import ROOT
from ROOT import gROOT
rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False


def jet_mass_producer(configs=None):
    """
    configs: configuration dict including:
    ModelName,gridHistFileName,QcdEstimation,channels
      -> channels: dict with dict for each channels:
        -> includes histLocation,histDir,samples,NormUnc,signal,obs,regions    
    """
    rebin_msd = True
    # min_msd, max_msd = (50,210)
    # binwidth = 16
    # nbins = int(np.floor((max_msd-min_msd)/binwidth))
    # msd_bins = np.linspace(min_msd, nbins*binwidth+min_msd, nbins+1)
    min_msd, max_msd = (50,190)
    binwidth = 4
    nbins = int(np.floor((max_msd-min_msd)/binwidth))
    print(nbins)
    msd_bins = np.linspace(min_msd, nbins*binwidth+min_msd, nbins+1)
    print(msd_bins)
    
    #channels for combined fit
    channels = configs['channels']
    qcd_estimation_channels = {k:v for k,v in channels.items() if "QcdEstimation" in v and v["QcdEstimation"]=="True"}
    
    print('channels:',channels.keys())
    
    #getting path of dir with root file from config
    hist_file = ROOT.TFile(configs['histLocation'])

    do_qcd_estimation =  len(qcd_estimation_channels)>0
    
    #specify if QCD estimation (using Bernstein-polynomial as TF) should be used
    ################
    #QCD Estimation#
    ################
    # derive pt bins from channel names for the pt,rho grid for the Bernstein-Polynomial
    if(do_qcd_estimation):
        # qcd_eff = get_qcd_efficiency(configs['histLocation'], w_channels)
        qcd_model = rl.Model('qcd_helper')
        qcd_pass, qcd_fail = 0.,0.
        for channel_name, config in qcd_estimation_channels.items():
            fail_ch = rl.Channel(channel_name + 'fail')
            pass_ch = rl.Channel(channel_name + 'pass')
            qcd_model.addChannel(fail_ch)
            qcd_model.addChannel(pass_ch)
            fail_hist = hist_file.Get('W_QCD__mjet_'+config['pt_bin']+'_fail')
            pass_hist = hist_file.Get('W_QCD__mjet_'+config['pt_bin']+'_pass')
            fail_ch.setObservation(fail_hist)
            pass_ch.setObservation(pass_hist)
            qcd_fail += fail_ch.getObservation().sum()
            qcd_pass += pass_ch.getObservation().sum()
        qcd_eff = qcd_pass / qcd_fail

        #get all lower edges from channel names
        pt_edges = [float(channel.split('Pt')[-1]) for channel in qcd_estimation_channels]
        #get last upper edge from name of last channel
        pt_edges.append(float(channels[list(qcd_estimation_channels.keys())[-1].split('Pt')[0]+'Pt%i'%pt_edges[-1]]['pt_bin'].split('to')[-1]))
        pt_bins = np.array(pt_edges)
        # pt_bins = np.array([500, 550, 600, 675, 800, 1200])
        n_pt = len(pt_bins) - 1
        msd = rl.Observable('msd',msd_bins)

        # here we derive these all at once with 2D array
        ptpts, msdpts = np.meshgrid(pt_bins[:-1] + 0.3 * np.diff(pt_bins), msd_bins[:-1] + 0.5 * np.diff(msd_bins), indexing='ij')
        rhopts = 2*np.log(msdpts/ptpts)
        ptscaled = (ptpts - 500.) / (1200. - 500.)
        rhoscaled = (rhopts - (-6)) / ((-2.1) - (-6))
        validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
        rhoscaled[~validbins] = 1  # we will mask these out later

    #get name from config, or fall back to default
    if('ModelName' in configs):
        model_name = configs['ModelName']
    else:
        model_name = 'Jet_Mass_Model'

    #Reading categories of consituent-variations for nuisance paramters from gridHist    
    grid_hist_file_name = configs['gridHistFileName']

    print('reading grid for nuisance parameter:')
    grid_hist_file = ROOT.TFile(grid_hist_file_name,'READ')
    grid_hist = grid_hist_file.Get('grid')
    grid_axes = dict(item.strip().split("=") for item in grid_hist.GetTitle().split(","))
    x_bins=range(grid_hist.GetNbinsX())
    y_bins=range(grid_hist.GetNbinsY())

    categories_hist=grid_hist_file.Get('categories')
    particle_categories=[]
    for i in range(1,categories_hist.GetNbinsX()+1):
        particle_categories.append(categories_hist.GetXaxis().GetBinLabel(i))

    grid_hist_file.Close()
    print('used variation categories:',particle_categories)
    print('X: %s , %i bins'%(grid_axes['x'],len(x_bins)))
    print('Y: %s , %i bins'%(grid_axes['y'],len(y_bins)))

    
    #setting up rhalphalib roofit model
    model = rl.Model(model_name)

    #setting up nuisances correspondig to consituent-variation according to categories from grid
    grid_nuisances = []
    print('adding nuisance paramters:')
    for category in particle_categories:
        for x_bin in x_bins:
            for y_bin in y_bins:
                print('massScale_%s%i_%s%i_%s'%(grid_axes['x'],x_bin,grid_axes['y'],y_bin,category),'shape')
                grid_nuisances.append([rl.NuisanceParameter('massScale_%s%i_%s%i_%s'%(grid_axes['x'],x_bin,grid_axes['y'],y_bin,category),'shape'),x_bin,y_bin,category])

    #setting up nuisances for systematic uncertainties
    print('CMS_lumi', 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    lumi_effect = 1.027 
    
    for channel_name, config in channels.items():
        print('setting up channel:',channel_name)
        
        #using hists with /variable/ in their name (default: Mass, if defined get from config) 
        variable = 'mjet' if 'variable' not in config else config['variable']        
        #getting list of samples from config
        samples = config['samples']
        #for WMass fit there are multiple regions per sample
        regions = [''] if 'regions' not in config else config['regions']

        print('getting template of variable:',variable)
        print('samples:',samples)
        print('regions:',regions)

        for region in regions:
            region_suffix = '_'+region if len(region)>0 else ''            
            hist_dir = config['selection']+'_%s__'+variable+'_%s'+config['pt_bin'] + region_suffix
            print('hist_dir:',hist_dir)
            #setting up channel for fit (name must be unique and can't include any '_')
            region_name = channel_name + region
            ch = rl.Channel(region_name)
            model.addChannel(ch)
            print('rl.Channel:',ch)

            for sample_name in samples:
                #do not include QCD template here, but rather use qcd estimation below
                if(('QcdEstimation' in config and config['QcdEstimation']=='True') and  'qcd' in sample_name.lower()):
                    continue

                #specify if sample is signal or background type
                sample_type = rl.Sample.SIGNAL if sample_name in config['signal'] else rl.Sample.BACKGROUND
                sample_hist = hist_file.Get(hist_dir%(sample_name,""))
                sample_hist.SetName('msd')

                #rebin hist
                if(rebin_msd > 0):
                    sample_hist = sample_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                
                #setup actual rhalphalib sample
                sample = rl.TemplateSample(ch.name + '_' + sample_name, sample_type, sample_hist)

                #setting effects of constituent variation nuisances (up/down)
                for grid_nuisance, x, y, category in grid_nuisances:
                    hist_up = hist_file.Get(hist_dir%(sample_name,str(x) + '_' + str(y) + '_' + category +'_')+ '__up')
                    hist_down = hist_file.Get(hist_dir%(sample_name,str(x) + '_' + str(y) + '_' + category +'_')+ '__down')

                    #rebin hists
                    if(rebin_msd > 0):
                        hist_up = hist_up.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                        hist_down = hist_down.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                        # hist_up = hist_up.Rebin(len(msd_bins)-1, hist_up.GetName() + "_" + channel_name, msd_bins)
                        # hist_down = hist_down.Rebin(len(msd_bins)-1, hist_down.GetName() + "_" + channel_name, msd_bins)

                    
                    sample.setParamEffect(grid_nuisance, hist_up, hist_down)
                sample.setParamEffect(lumi, lumi_effect)

                ch.addSample(sample)


            if 'varyPseudoLike' in config and config['obs'] == "Pseudo": #TODO: fix this for new HistFileFormat
                hist_path = config['varyPseudoLike']
                if '/' not in hist_path:
                    hist_path = hist_dir+'/'+hist_path
                data_hist = data_file.Get(hist_path)
                data_hist.SetName('Mass_central')
            else:
                #necessary for 2016 JetHT UHH Ntuples because they are missing JetConstituents
                # if('W' in channel_name and config['obs'] == "Data"):
                #     hist_path = hist_dir.split('pt')[0]+'noConst_pt'+hist_dir.split('pt')[1]
                # else:
                data_hist=hist_file.Get(hist_dir%("Data",""))

            if(rebin_msd > 0):
                data_hist = data_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            data_hist.SetName('msd')
            ch.setObservation(data_hist)
            if('QcdEstimation' in config and config['QcdEstimation']=='True'):
                mask = validbins[np.where(pt_bins==float(channel_name.split('Pt')[-1]))[0][0]]
                dropped_events = np.sum(ch.getObservation().astype(float)[~mask])
                percentage = dropped_events/np.sum(ch.getObservation().astype(float))
                print('dropping due to mask: %.2f events (out of %.2f -> %.2f%%)'%(dropped_events,np.sum(ch.getObservation().astype(float)),percentage*100))
                ch.mask = mask
                
    if(do_qcd_estimation):
        #QCD TF
        tf_params = rl.BernsteinPoly('tf_params', (2,2), ['pt','rho'], limits = (0,10))
        print('Using QCD efficiency (N2-ddt) of %.2f%% to scale initial QCD in pass region'%(qcd_eff*100))
        tf_params = qcd_eff * tf_params(ptscaled,rhoscaled)
        
        for channel_name, config in channels.items():
            if('QcdEstimation' not in config or config['QcdEstimation']=="False"):
                continue
            print(channel_name,'qcd estimation')
            fail_ch = model[channel_name + 'fail']
            pass_ch = model[channel_name + 'pass']
            ptbin = np.where(pt_bins==float(channel_name.split('Pt')[-1]))[0][0]
            qcd_params = np.array( [rl.IndependentParameter('qcdparam_ptbin%i_msdbin%i'%(ptbin,i),0) for i in range(msd.nbins)] )
            initial_qcd = fail_ch.getObservation().astype(float)
            for sample in fail_ch:
                initial_qcd -= sample.getExpectation(nominal=True)
            if np.any(initial_qcd<0.):
                raise ValueError('inital qcd (fail qcd from data - mc) negative at least one bin')
            sigmascale = 10.
            scaledparams = initial_qcd * ( 1 + sigmascale / np.maximum(1., np.sqrt(initial_qcd)))**qcd_params
            fail_qcd = rl.ParametericSample('%sfail_qcd' %channel_name, rl.Sample.BACKGROUND, msd, scaledparams)
            fail_ch.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('%spass_qcd'% channel_name, rl.Sample.BACKGROUND, tf_params[ptbin,:], fail_qcd)
            pass_ch.addSample(pass_qcd)
            
    model.renderCombine(model_name)


if(__name__ == "__main__"):
    import json,argparse
    import plotter
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with config")
    parser.add_argument('--justplots',action='store_true', help='just make plots.')
    parser.add_argument('--build', action='store_true', help='just build workspace for combine')
    args = parser.parse_args()
    
    try:
        configs = json.load(open(args.config))
    except IndexError:
        print("You must specify a configuration JSON!")
        sys.exit(0)
    if(not args.justplots):
        jet_mass_producer(configs)
        if(args.build):
            exit(0) 
        from runFit import runFits
        runFits([configs['ModelName']])

    plotter.plot_fit_result(configs)
    plotter.plot_fit_result(configs,logY=True)

    qcd_estimation_channels = {k:v for k,v in configs['channels'].items() if "QcdEstimation" in v and v["QcdEstimation"]=="True"}
    if(len(qcd_estimation_channels)>0):
        plotter.plot_qcd_bernstein(configs,(2,2))
