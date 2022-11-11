import os,sys
import ROOT
sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/')
import plotter,cms_style
cms_style.cms_style()

genbin_colors_hex=[
    #diverging color-scheme
    # "#a50026",
    # "#d73027",
    # "#f46d43",
    # "#fdae61",
    # "#fee090",
    
    
    # "#e0f3f8",
    # "#abd9e9",
    # "#74add1",
    # "#4575b4",
    # "#313695",

    # "#40004b",
    # "#762a83",
    # "#9970ab",
    # "#c2a5cf",
    # "#e7d4e8",

    # "#d9f0d3",
    # "#a6dba0",
    # "#5aae61",
    # "#1b7837",
    # "#00441b",


    #sequential colo-scheme
    '#fee6ce',
    '#fdd0a2',
    '#fdae6b',
    '#fd8d3c',
    '#f16913',
    '#d94801',
    '#a63603',
    '#7f2704',

    '#efedf5',
    '#dadaeb',
    '#bcbddc',
    '#9e9ac8',
    '#807dba',
    '#6a51a3',
    '#54278f',
    '#3f007d',

    '#e5f5e0',
    '#c7e9c0',
    '#a1d99b',
    '#74c476',
    '#41ab5d',
    '#238b45',
    '#006d2c',
    '#00441b',
]
genbin_TColors = [ROOT.TColor.GetColor(c) for c in genbin_colors_hex]

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
    plotter.private_work = True
    plotter.logY = logY
    plotter.additional_text_size_modifier = 1.5
    plotter.draw_extra_text = True
    plotter.luminosity = config.get('Luminosity',41.8)
    plotter.legend_bbox = (0.60,0.2,0.9,0.6)

    pseudo_data_info = config.get('Pseudo',[])
    if(len(pseudo_data_info)>0 and 'lumiScale' in pseudo_data_info[0]):
            lumiScale = float(pseudo_data_info[0].split(':')[-1])
            plotter.luminosity *= lumiScale
    for channel_str, channel in config['channels'].items():
        binning = config.get('binning',[50,300,26])
        plotter.yTitle = "Events / %i GeV"%binning[2]
        # plotter.yTitle = "Events / %i GeV"%config['binning'][channel['selection']][2]
        #print('plotting', channel_str)
        backgrounds = list(map(lambda bg: 'qcd' if ('QCD' in bg and "QcdEstimation" in channel and channel["QcdEstimation"]=="True") else bg ,plotter.mc_samples[channel['selection']]))
        if(use_config_samples):
            backgrounds = [bg.replace("QCD","qcd") for bg in channel['samples']]
        regions = channel['regions'] if 'regions' in channel else [""]
        for region in regions:
            for suffix in PrePostFit:
                hist_dir = (channel_str + region + '_' + suffix +'/%s') if do_postfit else ('shapes_prefit/'+channel_str+region+'/%s')
                plotter.rebin = False               
                h_obs,fit_shapes = plotter.get_hists(f_shapes,backgrounds,hist_dir, scaleQCD=False, selection=channel['selection'])
                #print(fit_shapes)
                # h_obs.GetXaxis().SetTitle('m_{SD} [GeV]')
                legend_entries=[]
                sample_counts = {}
                for sample, h in fit_shapes.items():
                    sample_name = sample
                    legend_suffix = ''
                    import re
                    genbin_pattern = re.compile('(?P<sample_name>[A-Za-z_]*)_ptgen(?P<ptgenbin>[0-9]+)_msdgen(?P<msdgenbin>[0-9]+)')
                    genbins = None
                    if genbin_pattern.match(sample):
                        regex_result = genbin_pattern.match(sample).groupdict()
                        sample_name = regex_result['sample_name']
                        genbins = (int(regex_result['ptgenbin']),int(regex_result['msdgenbin']))
                    if 'onegenbin' in sample:
                        sample_name = sample.replace("_onegenbin","")
                    if sample_name not in sample_counts:
                        sample_counts[sample_name] = 0
                    else:
                        sample_counts[sample_name] += 1
                    h.SetLineColor(plotter.colors.get(sample_name))
                    h.SetFillColorAlpha(plotter.colors.get(sample_name),.8)
                    if(genbins is not None):
                        # h.SetFillColorAlpha(plotter.colors.get(sample_name) + sample_counts[sample_name],.8)
                        h.SetFillColorAlpha(genbin_TColors[sample_counts[sample_name]],.8)
                        # h.SetLineColor(plotter.colors.get(sample_name) + sample_counts[sample_name])
                        legend_suffix = ' (ptgen%i_msdgen%i)'%genbins
                    # h.SetLineColorAlpha(plotter.colors.get(sample),.8)
                    legend_entries.append((h,plotter.legend_names.get(sample_name)+legend_suffix,'f'))
                
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
                
                
