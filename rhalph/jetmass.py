import sys,os
sys.path.append(os.getcwd()+'/rhalphalib/')
import rhalphalib as rl
import numpy as np
import ROOT
from ROOT import gROOT
rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def scale_lumi(hist,lumiscale):
    if(lumiscale == 1.0):
        return hist
    hist_old = hist.Clone()
    hist.Reset()
    for i in range(1,hist.GetNbinsX()+1):
        hist.SetBinContent(i,hist_old.GetBinContent(i)*lumiscale)
        hist.SetBinError(i,hist_old.GetBinError(i)*np.sqrt(lumiscale))
    return hist

def build_pseudo(samples,hist_file,hist_dir,vary_pseudo_like,MINIMAL_MODEL=False,seed=123456):
    lumiScale=1.0
    if(len(vary_pseudo_like)>0 and 'lumiScale' in vary_pseudo_like[0]):
        lumiScale = float(vary_pseudo_like[0].split(':')[-1])
    print('LumiScale:',lumiScale)
    if(len(vary_pseudo_like)>1):
        print('building pseudo data from ',samples,' (varied like',vary_pseudo_like,')')
        first_hist_dir = (hist_dir+'__%s')%(samples[0],*tuple(vary_pseudo_like[1:]))
        print(first_hist_dir)
        print(hist_file)
        pseudo_data_hist = hist_file.Get(first_hist_dir)
        for sample in samples[1:]:
            this_hist_dir = (hist_dir+'__%s')%(sample,*tuple(vary_pseudo_like[1:]))
            pseudo_data_hist.Add(hist_file.Get(this_hist_dir))
        pseudo_data_hist = scale_lumi(pseudo_data_hist,lumiScale)
        return pseudo_data_hist
    else:        
        print('building pseudo data from ',samples,' (unvaried)')
        toys=False
        if(len(vary_pseudo_like)>0 and 'toys' in vary_pseudo_like[0]):
            toys=True
        first_hist_dir = hist_dir%(samples[0],"")
        pseudo_data_hist = hist_file.Get(first_hist_dir)
        for sample in samples[1:]:
            if ('qcd' not in sample.lower() or "WJetsMatched" not in sample.lower()) and MINIMAL_MODEL:
                continue
            this_hist_dir = hist_dir%(sample,"")
            print(this_hist_dir)
            pseudo_data_hist.Add(hist_file.Get(this_hist_dir))

        if(toys):
            N_entries_toys = pseudo_data_hist.Integral()*lumiScale
            toy_pseudo_data = pseudo_data_hist.Clone()
            toy_pseudo_data.Reset()
            n_toys=100
            for i in range(n_toys):
                ROOT.gRandom.SetSeed(seed+i)
                toy_pseudo_data.FillRandom(pseudo_data_hist,int(N_entries_toys/n_toys))
            toy_pseudo_data.Scale(N_entries_toys/toy_pseudo_data.Integral())
            return toy_pseudo_data
        else:
            pseudo_data_hist = scale_lumi(pseudo_data_hist,lumiScale) 
            return pseudo_data_hist
        
def build_mass_scale_variations(grid_hist_file_name):
    
    # print('reading grid for nuisance parameter:')
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
    # print('used variation categories:',particle_categories)
    # print('X: %s , %i bins'%(grid_axes['x'],len(x_bins)))
    # print('Y: %s , %i bins'%(grid_axes['y'],len(y_bins)))

    #setting up nuisances correspondig to consituent-variation according to categories from grid
    grid_nuisances = []
    mass_scale_names = []
    # print('adding nuisance paramters:')
    for category in particle_categories:
        for x_bin in x_bins:
            for y_bin in y_bins:
                scale_name = 'massScale_%s%i_%s%i_%s'%(grid_axes['x'],x_bin,grid_axes['y'],y_bin,category)
                # print(scale_name,'shape')
                mass_scale_names.append(scale_name)
                grid_nuisances.append([rl.NuisanceParameter('massScale_%s%i_%s%i_%s'%(grid_axes['x'],x_bin,grid_axes['y'],y_bin,category),'shape',0,-10,10),x_bin,y_bin,category])
    return grid_nuisances, mass_scale_names

def jet_mass_producer(args,configs=None,MINIMAL_MODEL=False,includeMassScales=True):
    """
    configs: configuration dict including:
    ModelName,gridHistFileName,channels,histLocation
      -> channels: dict with dict for each channels:
        -> includes histDir,samples,NormUnc,signal,regions,QcdEstimation
    """
    rebin_msd = True
    binnings = {"W":np.linspace(50,300,26),"top":np.linspace(50,300,26)}
    binning_from_config = configs.get('binning',{})
    for selection,bin_info in binning_from_config.items():
        min_msd, max_msd = (bin_info[0],bin_info[1])
        binwidth = bin_info[2]
        nbins = int(np.floor((max_msd-min_msd)/binwidth))
        msd_bins = np.linspace(min_msd, nbins*binwidth+min_msd, nbins+1)
        binnings[selection] = msd_bins
        
    #channels for combined fit
    channels = configs['channels']
    qcd_estimation_channels = {k:v for k,v in channels.items() if "QcdEstimation" in v and v["QcdEstimation"]=="True"}
    
    print('channels:',channels.keys())
    
    #getting path of dir with root file from config
    hist_file = ROOT.TFile(configs['histLocation'])

    do_qcd_estimation =  len(qcd_estimation_channels)>0
    do_initial_qcd_fit = (configs.get("InitialQCDFit","False") == "True")
    qcd_fail_region_constant = (configs.get("QCDFailConstant","False") == "True")

    lumi_scale = 1.
    if('Pseudo' in configs and len(configs['Pseudo'])>0 and  'lumiScale' in configs['Pseudo'][0]):
        lumi_scale = float(configs['Pseudo'][0].split(':')[-1])
    model_name=  configs.get('ModelName','Jet_Mass_Model')  #get name from config, or fall back to default

    #specify if QCD estimation (using Bernstein-polynomial as TF) should be used
    ################
    #QCD Estimation#
    ################
    # derive pt bins from channel names for the pt,rho grid for the Bernstein-Polynomial
    if(do_qcd_estimation):
        print('Doing some preparations for data driven QCD Estimate (Bernstein TF)')
        bernstein_orders = tuple(configs.get('BernsteinOrders',[2,2]))
        qcd_model = rl.Model('qcdmodel')
        qcd_pass, qcd_fail = 0.,0.
        qcd_estimation_relevant_selection = 'W'
        for channel_name, config in qcd_estimation_channels.items():
            qcd_estimation_relevant_selection = config['selection']
            msd_bins = binnings[qcd_estimation_relevant_selection]
            fail_ch = rl.Channel(channel_name + 'fail')
            pass_ch = rl.Channel(channel_name + 'pass')
            qcd_model.addChannel(fail_ch)
            qcd_model.addChannel(pass_ch)
            additional_bin = config.get('additional_bin','')
            fail_hist = hist_file.Get('W_QCD__mjet_'+config['pt_bin']+ additional_bin+'_fail')
            pass_hist = hist_file.Get('W_QCD__mjet_'+config['pt_bin']+ additional_bin+'_pass')            
            if(rebin_msd > 0):
                fail_hist = fail_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                pass_hist = pass_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            if(lumi_scale != 1.0):
                fail_hist = scale_lumi(fail_hist,lumi_scale) 
                pass_hist = scale_lumi(pass_hist,lumi_scale)

            empty_hist = fail_hist.Clone()
            empty_hist.Reset()
            signal_fail = rl.TemplateSample(channel_name + 'fail' + '_' + 'Signal', rl.Sample.SIGNAL, empty_hist)
            fail_ch.addSample(signal_fail)
            signal_pass = rl.TemplateSample(channel_name + 'pass' + '_' + 'Signal', rl.Sample.SIGNAL, empty_hist)
            pass_ch.addSample(signal_pass)

            fail_ch.setObservation(fail_hist)
            pass_ch.setObservation(pass_hist)
            qcd_fail += fail_ch.getObservation().sum()
            qcd_pass += pass_ch.getObservation().sum()
        qcd_eff = qcd_pass / qcd_fail
        #get all lower edges from channel names
        # pt_edges = [float(channel.split('Pt')[-1]) for channel in qcd_estimation_channels]
        # #get last upper edge from name of last channel
        # pt_edges.append(float(channels[list(qcd_estimation_channels.keys())[-1].split('Pt')[0]+'Pt%i'%pt_edges[-1]]['pt_bin'].split('to')[-1]))
        pt_edges = configs.get('pt_edges',[500,550,600,675,800,1200])
        pt_bins = np.array(pt_edges)
        # pt_bins = np.array([500, 550, 600, 675, 800, 1200])
        n_pt = len(pt_bins) - 1
        msd_bins = binnings[qcd_estimation_relevant_selection]
        msd = rl.Observable('msd',msd_bins)

        # here we derive these all at once with 2D array
        ptpts, msdpts = np.meshgrid(pt_bins[:-1] + 0.3 * np.diff(pt_bins), msd_bins[:-1] + 0.5 * np.diff(msd_bins), indexing='ij')
        rhopts = 2*np.log(msdpts/ptpts)
        ptscaled = (ptpts - 500.) / (1200. - 500.)
        rhoscaled = (rhopts - (-6)) / ((-2.1) - (-6))
        validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
        rhoscaled[~validbins] = 1  # we will mask these out later

        TF_suffix = configs.get('TFSuffix',"")
        
        if(do_initial_qcd_fit):
            initial_qcd_fit_orders = tuple(configs.get('InitialQCDFitOrders',[2,2]))        
            if not os.path.exists(model_name):
                os.makedirs(model_name)
            print('QCD eff:',qcd_eff)
            # tf_MCtempl = rl.BernsteinPoly("tf_MCtempl", initial_qcd_fit_orders, ['pt', 'rho'], init_params = np.ones((initial_qcd_fit_orders[0]+1,initial_qcd_fit_orders[1]+1)), limits=(-1,10))
            tf_MCtempl = rl.BernsteinPoly("tf_MCtempl_"+model_name+TF_suffix, initial_qcd_fit_orders, ['pt', 'rho'], init_params = np.ones((initial_qcd_fit_orders[0]+1,initial_qcd_fit_orders[1]+1)), limits=(-50,50))
            tf_MCtempl_params = qcd_eff * tf_MCtempl(ptscaled, rhoscaled)
            for channel_name, config in channels.items():
                # ptbin = np.where(pt_bins==float(channel_name.split('Pt')[-1]))[0][0]
                ptbin = np.where(pt_bins==float(config['pt_bin'].split('to')[0]))[0][0]
                failCh = qcd_model[channel_name + 'fail']
                passCh = qcd_model[channel_name + 'pass']
                failObs = failCh.getObservation()
                if(qcd_fail_region_constant):
                    print("Setting QCD parameters in fail region constant")
                qcdparams = np.array([rl.IndependentParameter('qcdparam_'+model_name+TF_suffix+'_ptbin%d_msdbin%d' % (ptbin, i), 0, constant=qcd_fail_region_constant) for i in range(msd.nbins)])

                sigmascale = 10.
                scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
                fail_qcd = rl.ParametericSample('%sfail_qcd' %channel_name, rl.Sample.BACKGROUND, msd, scaledparams)
                failCh.addSample(fail_qcd)
                pass_qcd = rl.TransferFactorSample('%spass_qcd' %channel_name, rl.Sample.BACKGROUND, tf_MCtempl_params[ptbin, :], fail_qcd)
                passCh.addSample(pass_qcd)
                
                failCh.mask = validbins[ptbin]
                passCh.mask = validbins[ptbin]

            qcd_model.renderCombine(model_name+"/qcdmodel")

            qcdfit_ws = ROOT.RooWorkspace('w')
            simpdf, obs = qcd_model.renderRoofit(qcdfit_ws)
            ROOT.Math.MinimizerOptions.SetDefaultPrecision(1e-18)
            # ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
            # ROOT.Math.MinimizerOptions.SetDefaultTolerance(0.0001)
            # ROOT.Math.MinimizerOptions.SetDefaultPrecision(-1.0)
            qcdfit = simpdf.fitTo(obs,
                                  ROOT.RooFit.Extended(True),
                                  ROOT.RooFit.SumW2Error(True),
                                  ROOT.RooFit.Strategy(1),
                                  ROOT.RooFit.Save(),
                                  ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                                  # ROOT.RooFit.PrintLevel(-1),                              
                                  ROOT.RooFit.PrintLevel(1),
                                  ROOT.RooFit.Minos(0)
            )

            qcdfit_ws.add(qcdfit)
            if "pytest" not in sys.modules:
                qcdfit_ws.writeToFile(model_name+ '/qcdfit_'+model_name+TF_suffix+'.root')
            if qcdfit.status() != 0:
                raise RuntimeError('Could not fit qcd')

            qcd_model.readRooFitResult(qcdfit)
            
            param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
            decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', qcdfit, param_names)
            tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
            tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)
            tf_dataResidual = rl.BernsteinPoly("tf_dataResidual_"+model_name+TF_suffix, bernstein_orders, ['pt', 'rho'], limits=(-50,50))
            # tf_dataResidual = rl.BernsteinPoly("tf_dataResidual", bernstein_orders, ['pt', 'rho'], limits=(0,10))
            tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
            tf_params = qcd_eff * tf_MCtempl_params_final * tf_dataResidual_params
        else:
            tf_params = None # define later
        

    #Reading categories of consituent-variations for nuisance paramters from gridHist    

    grid_nuisances, _ =build_mass_scale_variations(configs['gridHistFileName'])
    
    #setting up rhalphalib roofit model
    model = rl.Model(model_name)

    #setting up nuisances for systematic uncertainties
    print('CMS_lumi', 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    lumi_effect = 1.027 

    norm_nuisances = {}
    for channel_name in channels.keys():
        if( MINIMAL_MODEL ):
            break
        for i, sample in enumerate(channels[channel_name]['samples']):
            if 'NormUnc' not in channels[channel_name]:
                continue
            norm_uncertainties = channels[channel_name]['NormUnc']
            for name,norm_unc in norm_uncertainties.items():
                nuisance_par = [rl.NuisanceParameter(name+'_normUnc','lnN'),norm_unc]
                for k,v in norm_nuisances.items():
                    if name in v[0].name:
                        nuisance_par = v
                if norm_unc>0 and name in sample and  sample not in norm_nuisances:
                    norm_nuisances.update({sample:nuisance_par})
                
    
    for channel_name, config in channels.items():
        print('setting up channel:',channel_name)
        
        #using hists with /variable/ in their name (default: Mass, if defined get from config) 
        variable = 'mjet' if 'variable' not in config else config['variable']        
        #getting list of samples from config
        if MINIMAL_MODEL:
            config['samples'] = ['QCD','WJetsMatched'] 
        samples =  config['samples']
        #for WMass fit there are multiple regions per sample
        regions = [''] if 'regions' not in config else config['regions']

        print('getting template of variable:',variable)
        print('samples:',samples)
        print('regions:',regions)
        msd_bins = binnings[config['selection']]

        for region in regions:
            additional_bin = config.get('additional_bin','')
            region_suffix = '_'+region if len(region)>0 else ''            
            hist_dir = config['selection']+'_%s__'+variable+'_%s'+config['pt_bin'] + additional_bin + region_suffix
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
                print(hist_dir%(sample_name,""))
                sample_hist.SetName('msd')

                #rebin hist
                if(rebin_msd > 0):
                    sample_hist = sample_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                if(lumi_scale != 1.0):
                    sample_hist = scale_lumi(sample_hist,lumi_scale) 

                #setup actual rhalphalib sample
                sample = rl.TemplateSample(ch.name + '_' + sample_name, sample_type, sample_hist)                
                #sample.autoMCStats()
                #setting effects of constituent variation nuisances (up/down)
                for grid_nuisance, x, y, category in grid_nuisances:
                    hist_up = hist_file.Get(hist_dir%(sample_name,str(x) + '_' + str(y) + '_' + category +'_')+ '__up')
                    hist_down = hist_file.Get(hist_dir%(sample_name,str(x) + '_' + str(y) + '_' + category +'_')+ '__down')

                    #rebin hists
                    if(rebin_msd > 0):
                        hist_up = hist_up.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                        hist_down = hist_down.Rebin(len(msd_bins)-1, 'msd', msd_bins)
                    if(lumi_scale != 1.0):
                        hist_up = scale_lumi(hist_up,lumi_scale) 
                        hist_down = scale_lumi(hist_down,lumi_scale) 

                    if(includeMassScales):
                        sample.setParamEffect(grid_nuisance, hist_up, hist_down)                    
                sample.setParamEffect(lumi, lumi_effect)
                if sample_name in norm_nuisances.keys():
                    sample.setParamEffect(norm_nuisances[sample_name][0],norm_nuisances[sample_name][1])
                    
                ch.addSample(sample)

            PseudoData = 'Pseudo' in configs and len(configs['Pseudo'])>0
            if PseudoData:
                data_hist = build_pseudo(samples,hist_file,hist_dir,configs['Pseudo'],MINIMAL_MODEL) 
            else:
                print('using data!!!!!')
                data_hist=hist_file.Get(hist_dir%("Data",""))

            if(rebin_msd > 0):
                data_hist = data_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            data_hist.SetName('msd')
            ch.setObservation(data_hist,read_sumw2=PseudoData)
            if('QcdEstimation' in config and config['QcdEstimation']=='True'):
                mask = validbins[np.where(pt_bins==float(config['pt_bin'].split('to')[0]))[0][0]]
                # dropped_events = np.sum(ch.getObservation().astype(float)[~mask])
                # percentage = dropped_events/np.sum(ch.getObservation().astype(float))
                # print('dropping due to mask: %.2f events (out of %.2f -> %.2f%%)'%(dropped_events,np.sum(ch.getObservation().astype(float)),percentage*100))
                ch.mask = mask
                
    if(do_qcd_estimation):
        #QCD TF
        if(not do_initial_qcd_fit):
            tf_params = rl.BernsteinPoly('tf_params_'+model_name+TF_suffix, bernstein_orders, ['pt','rho'], limits = (-50,50))
            print('Using QCD efficiency (N2-ddt) of %.2f%% to scale initial QCD in pass region'%(qcd_eff*100))
            tf_params = qcd_eff * tf_params(ptscaled,rhoscaled)
        
        for channel_name, config in channels.items():
            if('QcdEstimation' not in config or config['QcdEstimation']=="False"):
                continue
            print(channel_name,'qcd estimation')
            fail_ch = model[channel_name + 'fail']
            pass_ch = model[channel_name + 'pass']
            ptbin = np.where(pt_bins==float(config['pt_bin'].split('to')[0]))[0][0]
            if(qcd_fail_region_constant):
                print("Setting QCD parameters in fail region constant")
            qcd_params = np.array( [rl.IndependentParameter('qcdparam_'+model_name+TF_suffix+'_ptbin%i_msdbin%i'%(ptbin,i),0, constant=qcd_fail_region_constant) for i in range(msd.nbins)] )

            initial_qcd = fail_ch.getObservation()[0].astype(float) if isinstance(fail_ch.getObservation(), tuple) else fail_ch.getObservation().astype(float)
            for sample in fail_ch:
                initial_qcd -= sample.getExpectation(nominal=True)
            if np.any(initial_qcd<0.):
                initial_qcd  = np.where(initial_qcd<=0.,0,initial_qcd)
                print('negative bins in initial_qcd in ',channel_name)
                # continue
                minimum = np.amin(initial_qcd)
                initial_qcd = np.where(initial_qcd == 0,minimum,initial_qcd)
                initial_qcd += abs(minimum)
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
    import fitplotter
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with config")
    parser.add_argument('--justplots',action='store_true', help='just make plots.')
    parser.add_argument('--build', action='store_true', help='just build workspace for combine')
    parser.add_argument('--job_index', default="", type=str)
    parser.add_argument('--minimalModel',action='store_true')
    parser.add_argument('--noMassScales',action='store_true')
    parser.add_argument('--defaultPOI',action='store_true')
    parser.add_argument('--skipTemplatePlots',action='store_true')
    parser.add_argument('--customCombineWrapper',action='store_true')
    parser.add_argument('--combineOptions',type=str,help="string with additional cli-options passed to combine",default="")
    # parser.add_argument('--QCDOnly',action="store_true",help="perform full TF fit to QCD MC - for BernsteinOptimization only")
    args = parser.parse_args()
    
    try:
        configs = json.load(open(args.config))
    except IndexError:
        print("You must specify a configuration JSON!")
        sys.exit(0)
    configs['ModelName'] = configs['ModelName']+str(args.job_index)

    if(not args.justplots):
        jet_mass_producer(args,configs,args.minimalModel,(not args.noMassScales))
        open(configs['ModelName']+'/config.json','w').write(json.dumps(configs,sort_keys=False,indent=2))
        use_r_poi = args.defaultPOI or args.noMassScales
        if not args.customCombineWrapper:
            _,mass_scale_names = build_mass_scale_variations(configs['gridHistFileName'])
            from CombineWorkflows import CombineWorkflows
            cw = CombineWorkflows()
            cw.workspace = configs['ModelName']+'/'+configs['ModelName']+'_combined.root'
            cw.extraOptions = "--freezeParameters r --preFitValue 0 " + args.combineOptions
            cw.POI = "r" if use_r_poi else mass_scale_names
            cw.method = "diagnostics"
            cw.write_wrapper()
        else:
            if(not os.path.isfile(configs['ModelName']+'/wrapper.sh')):
                import warnings
                warnings.warn("\033[93mYou used the option --CustomCombineWrapper, but no wrapper can be found in the modeldir!\033[0m",RuntimeWarning)
            #write_wrapper([configs['ModelName']],"r" if use_r_poi else mass_scale_names,additional_options=args.combineOptions,combineWorkflow=args.combineWorkflow)
        if(args.build):
            exit(0)     
        from runFit import runFits        
        runFits([configs['ModelName']])

    if(args.customCombineWrapper):
        exit(0)
        
    do_postfit = True
    try:
        fit_diagnostics = ROOT.TFile(configs['ModelName']+'/fitDiagnostics.root',"READ")
        fit_result = fit_diagnostics.Get("fit_s")
        massScales = []
        fitargs = fit_diagnostics.Get("fit_s").floatParsFinal()
        for name in build_mass_scale_variations(configs['gridHistFileName'])[1]:
            param = fitargs.find(name)
            center = param.getValV()
            error_up = abs(param.getErrorHi())
            error_down = abs(param.getErrorLo())
            
            massScales.append([center,error_up,-error_down])
        np.save(configs['ModelName']+"/"+configs['ModelName']+'MassScales.npy',np.array(massScales))
        do_postfit = fit_result.status() <= 1
    except:
        print("fit failed. only plotting prefit distributions from fitDiangnostics (beware weird shape uncertainties)")
        
        do_postfit = False
    if(not args.skipTemplatePlots):
        fitplotter.plot_fit_result(configs,plot_total_sig_bkg = False,do_postfit=do_postfit)
        fitplotter.plot_fit_result(configs,logY=True,plot_total_sig_bkg = False,do_postfit=do_postfit)
    if(do_postfit and not args.noMassScales):
        fitplotter.plot_mass_scale_nuisances(configs)
    
    qcd_estimation_channels = {k:v for k,v in configs['channels'].items() if "QcdEstimation" in v and v["QcdEstimation"]=="True"}
    if(len(qcd_estimation_channels)>0 and do_postfit):
        fitplotter.plot_qcd_bernstein(configs,do_3d_plot = False)
        if(configs.get('QCDFailConstant','False') == 'False'):
            fitplotter.plot_qcd_fail_parameters(configs)
