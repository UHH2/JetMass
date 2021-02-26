import os,sys
import ROOT
sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/python/')
import plotter,cms_style
cms_style.cms_style()

def plot_fit_result(config={'ModelName':'WMassModel'},logY=False, fit_shapes_root="fit_shapes.root",plot_total_sig_bkg=True,do_postfit=True,use_config_samples=False): 
    print('opening file',config['ModelName']+"/"+fit_shapes_root)
    fit_shapes_root = "fit_shapes.root" if do_postfit else "fitDiagnostics.root"
    f_shapes = ROOT.TFile(config['ModelName']+"/"+fit_shapes_root,"READ")
    ROOT.TH1.AddDirectory(0)
    PrePostFit = ['prefit','postfit'] if do_postfit else ['prefit']
    out_dir = config['ModelName']+'/plots/'+fit_shapes_root.replace('.root','/')+('/logY/' if logY else '')
    if(not os.path.exists(out_dir)):
        os.makedirs(out_dir)

    fit_hist_names={
        'fit_shapes':{},
        'fitDiagnostics':{},
    }
        
    plotter.xTitle = "m_{SD} [GeV]"
    plotter.logY = logY
    plotter.additional_text_size_modifier = 1.3
    plotter.draw_extra_text = False
    plotter.luminosity = config.get('Luminosity',41.8)

    pseudo_data_info = config.get('Pseudo',[])
    if(len(pseudo_data_info)>0 and 'lumiScale' in pseudo_data_info[0]):
            lumiScale = float(pseudo_data_info[0].split(':')[-1])
            plotter.luminosity *= lumiScale
    for channel_str, channel in config['channels'].items():
        binning = config.get('binning',[50,300,26])
        plotter.yTitle = "Events / %i GeV"%binning[2]
        # plotter.yTitle = "Events / %i GeV"%config['binning'][channel['selection']][2]
        print('plotting', channel_str)
        backgrounds = list(map(lambda bg: 'qcd' if ('QCD' in bg and "QcdEstimation" in channel and channel["QcdEstimation"]=="True") else bg ,plotter.mc_samples[channel['selection']]))
        if(use_config_samples):
            backgrounds = [bg.replace("QCD","qcd") for bg in channel['samples']]
        regions = channel['regions'] if 'regions' in channel else [""]
        for region in regions:
            for suffix in PrePostFit:
                hist_dir = (channel_str + region + '_' + suffix +'/%s') if do_postfit else ('shapes_prefit/'+channel_str+region+'/%s')
                plotter.rebin = False               
                h_obs,fit_shapes = plotter.get_hists(f_shapes,backgrounds,hist_dir)
                print(fit_shapes)
                # h_obs.GetXaxis().SetTitle('m_{SD} [GeV]')
                legend_entries=[]
                for sample, h in fit_shapes.items():
                    h.SetFillColorAlpha(plotter.colors.get(sample),.8)
                    # h.SetLineColorAlpha(plotter.colors.get(sample),.8)
                    h.SetLineColor(plotter.colors.get(sample))
                    legend_entries.append((h,plotter.legend_names.get(sample),'f'))
                
                additional_hists=[]
                if(plot_total_sig_bkg):
                    h_total_b = f_shapes.Get(hist_dir%('TotalBkg' if do_postfit else 'total_background'))
                    h_total_b.SetLineColor(46)
                    h_total_b.SetLineWidth(2)
                    legend_entries.append((h_total_b,'Total Bkg.','l'))
                    additional_hists.append(h_total_b)

                    h_total_s = f_shapes.Get(hist_dir%('TotalSig' if do_postfit else 'total_signal'))
                    h_total_s.SetLineColor(32)
                    h_total_s.SetLineWidth(2)
                    legend_entries.append((h_total_s,'Total Signal','l'))
                    additional_hists.append(h_total_s)

                legend_entries.append((h_obs,'Data','pe1x0'))
                additional_text = ' '+plotter.selection_tex[channel['selection']]+"\\"+plotter.pt_bins_tex_dict[channel['pt_bin']]+"\\ %s #bf{%s}"%(plotter.region_tex[channel['selection']][region],suffix)
                
                plotter.plot_data_mc(h_obs,fit_shapes,channel_str+'_'+region+'_'+suffix,out_dir,legend_entries=legend_entries,additional_hists=additional_hists,additional_text=additional_text)
                
                
