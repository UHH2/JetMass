import sys,os
sys.path.append(os.getcwd()+'/rhalphalib/')
import rhalphalib as rl
import numpy as np
import ROOT


def jet_mass_producer(configs=None):
    """
    configs: configuration dict including:
    ModelName,gridHistFileName,pathCMSSW,QcdEstimation,channels
      -> channels: dict with dict for each channels:
        -> includes histLocation,histDir,samples,NormUnc,signal,obs,regions    
    """
    #specify how many bins you want the softdropmass hists to be rebinned to. (-1: no rebinning)
    rebin_msd = 26
    n_msd_bins = rebin_msd if rebin_msd > 0 else 24
    # msd_bins = np.linspace(50,170,n_msd_bins)
    msd_bins = np.linspace(50,200,n_msd_bins)

    #channels for combined fit
    channels = configs['channels']
    w_channels = {k:v for k,v in channels.items() if "WMass" in k}

    print('channels:',channels.keys())

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
        
        #getting path of dir with root file from config
        hist_location = config['histLocation']
        #using hists with /variable/ in their name (default: Mass, if defined get from config) 
        variable = 'Mass' if 'variable' not in config else config['variable']        
        #getting list of samples from config
        samples = config['samples']
        #for WMass fit there are multiple regions per sample
        regions = [''] if 'regions' not in config else config['regions']

        print('getting template of variable:',variable)
        print('from dir:',hist_location)
        print('samples:',samples)
        print('regions:',regions)

        for region in regions:
            region_suffix = '_'+region if len(region)>0 else ''
            hist_dir = config['histDir'] + region_suffix
            print('hist_dir:',hist_dir)
            #setting up channel for fit (name must be unique and can't include any '_')
            region_name = channel_name + region
            ch = rl.Channel(region_name)
            model.addChannel(ch)
            print('rl.Channel:',ch)

            for sample_name in samples:
                #do not include QCD template here, but rather use qcd estimation below
                if('qcd' in sample_name.lower()):
                    continue

                #specify if sample is signal or background type
                sample_type = rl.Sample.SIGNAL if sample_name == config['signal'] else rl.Sample.BACKGROUND
                sample_file = ROOT.TFile(hist_location + '/' + sample_name + '.root','READ')
                sample_hist = sample_file.Get(hist_dir + '/' + variable + '_central')

                #rebin hist
                if(rebin_msd > 0):
                    sample_hist = sample_hist.Rebin(len(msd_bins)-1, sample_hist.GetName() + "_" + channel_name, msd_bins)
                
                #setup actual rhalphalib sample
                sample = rl.TemplateSample(ch.name + '_' + sample_name, sample_type, sample_hist)

                #setting effects of constituent variation nuisances (up/down)
                for grid_nuisance, x, y, category in grid_nuisances:
                    hist_up = sample_file.Get(hist_dir + '/' + variable + '_' + str(x) + '_' + str(y) + '_' + category + '_up')
                    hist_down = sample_file.Get(hist_dir + '/' + variable + '_' + str(x) + '_' + str(y) + '_' + category + '_down')

                    #rebin hists
                    if(rebin_msd > 0):
                        hist_up = hist_up.Rebin(len(msd_bins)-1, hist_up.GetName() + "_" + channel_name, msd_bins)
                        hist_down = hist_down.Rebin(len(msd_bins)-1, hist_down.GetName() + "_" + channel_name, msd_bins)

                    
                    sample.setParamEffect(grid_nuisance, hist_up, hist_down)
                sample.setParamEffect(lumi, lumi_effect)

                ch.addSample(sample)

            data_file = ROOT.TFile(hist_location+ '/' +config['obs'] +'.root')

            if 'varyPseudoLike' in config and config['obs'] == "Pseudo":
                hist_path = config['varyPseudoLike']
                if '/' not in hist_path:
                    hist_path = hist_dir+'/'+hist_path
                data_hist = data_file.Get(hist_path)
                data_hist.SetName('Mass_central')
            else:
                #necessary for 2016 JetHT UHH Ntuples because they are missing JetConstituents
                if('W' in channel_name and config['obs'] == "Data"):
                    hist_path = hist_dir.split('pt')[0]+'noConst_pt'+hist_dir.split('pt')[1]
                else:
                    hist_path = hist_dir
                data_hist=data_file.Get(hist_path+'/'+'Mass_central')

            if(rebin_msd > 0):
                data_hist = data_hist.Rebin(len(msd_bins)-1, data_hist.GetName() + "_" + channel_name, msd_bins)

            ch.setObservation(data_hist)

        model.renderCombine(model_name)


if(__name__ == "__main__"):
    import json,sys

    try:
        configs = json.load(open(sys.argv[1]))
    except IndexError:
        print("You must specify a configuration JSON!")
        sys.exit(0)

    jet_mass_producer(configs)

    from runFit import runFits
    # if('pathCMSSW' in configs):
    #     runFits([configs['ModelName']],configs['pathCMSSW'])
    # else:
    runFits([configs['ModelName']])
    
