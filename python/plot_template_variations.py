from __future__ import print_function
import ROOT, os
import UHH2.JetMass.plotter as plotter

def plot_template_variations(c_name,c_grid,file_name='../macros/Histograms.root',selection='top',output_dir='template_variations'):
    f_hists = ROOT.TFile(file_name,'READ')
    if(not os.path.exists(output_dir)):
        os.makedirs(output_dir)

    plotter.obs_draw_option='Hist'
    plotter.ratio_draw_options='H'
    plotter.ratio_hist_yTitle="#frac{variation}{nominal}"
    plotter.extra_text = "Simulation"
    
    for ipt_bin in range(len(plotter.pt_bins_dict[selection])):
        pt_bin = plotter.pt_bins_dict[selection][ipt_bin]
        for region in ['_pass','_fail']:

            region_out_dir = output_dir+'/'+pt_bin+'_'+region
            if(not os.path.exists(region_out_dir)):
                os.makedirs(region_out_dir)
                
            for i in range(len(c_grid['ptbins'])-1):
                for j in range(len(c_grid['etabins'])-1):
                    for pfflavour in c_grid['pfflavours']:
                        hist_dir = selection+'_%s__mjet_'+pt_bin+region
                        h_pseudo,h_mc = plotter.get_hists(f_hists,plotter.mc_samples[selection],hist_dir,pseudo_data=True)
                        varied_hists = {'up':(),'down':()}
                        for direction in ['up','down']:
                            category = '%i_%i_%s'%(i,j,pfflavour)
                            hist_dir_varied = selection+'_%s__mjet_'+category+'_'+pt_bin+region+'__'+direction
                            varied_hists[direction] = plotter.get_hists(f_hists,plotter.mc_samples[selection],hist_dir_varied,pseudo_data=True)
                        
                        for sample,hist in h_mc.items():
                            variations = [varied_hists[direction][1][sample] for direction in ['up','down']]
                            variations[0].SetLineColor(ROOT.kBlue+2)
                            variations[1].SetLineColor(ROOT.kRed+2)

                            legend_entries=[(hist,'nominal','l')]
                            legend_entries.append((variations[0],'+1% variation','l'))
                            legend_entries.append((variations[1],'-1% variation','l'))
                            plotter.plot_data_mc(h_data=hist, h_mc=None ,plot_title=sample + ' '+pt_bin+' '+region+' '+ category,out_dir=region_out_dir,additional_hists=variations,legend_entries=legend_entries,additional_text=sample+'\\'+plotter.pt_bins_tex_dict[pt_bin]+'\\'+category)
                        
                        variations = [varied_hists[direction][0] for direction in ['up','down']]
                        variations[0].SetLineColor(ROOT.kBlue+2)
                        variations[1].SetLineColor(ROOT.kRed+2)

                        legend_entries=[(h_pseudo,'nominal','l')]
                        legend_entries.append((variations[0],'+1% variation','l'))
                        legend_entries.append((variations[1],'-1% variation','l'))
                        plotter.plot_data_mc(h_data=h_pseudo, h_mc=None ,plot_title='PseudoData' + ' '+pt_bin+' '+region+' '+ category,out_dir=region_out_dir,additional_hists=variations,legend_entries=legend_entries,additional_text='PseudoData\\'+plotter.pt_bins_tex_dict[pt_bin]+'\\'+category)
                        
if(__name__ ==  '__main__'):
    
    selection = 'top'    

    category_grids={
        "PF_flavours":{
            "ptbins":[0.,100000.],
            "etabins":[0.,9.],
            "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        }
    }

    output_dir = "../Plots/template_variations"
    
    file_name='../macros/scaleStudy/Histograms_PF_flavours.root'

    for c_name,c_grid in category_grids.items():
        plot_template_variations(c_name,c_grid,file_name,selection,output_dir)
