from __future__ import print_function
import ROOT
import os
import UHH2.JetMass.plotter as plotter


def plot_template_variations(c_name,
                             c_grid,
                             file_name: str = '../macros/Histograms.root',
                             selection: str = 'top',
                             output_dir: str = 'template_variations',
                             variation: int = 1
                             ):
    f_hists = ROOT.TFile(file_name, 'READ')
    if (not os.path.exists(output_dir)):
        os.makedirs(output_dir)

    plotter.obs_draw_option = 'Hist'
    plotter.ratio_draw_options = 'H'
    plotter.ratio_hist_yTitle = "#frac{variation}{nominal}"
    plotter.extra_text = "private work"
    plotter.private_work = True
    plotter.legend_bbox = (0.60, 0.2, 0.9, 0.6)

    # colors = {"chargedH":ROOT.kGreen+2, "neutralH":ROOT.kRed-4, "gamma":ROOT.kOrange, "other":ROOT.kBlue+2}

    # colors = {"chargedH":ROOT.kGreen+2, "neutralH":ROOT.kRed+1, "gamma":ROOT.kBlue, "other":ROOT.kOrange}
    colors = {"chargedH": ROOT.kRed-4, "neutralH": ROOT.kAzure+7, "gamma": 798, "other": 13, "all": ROOT.kRed-4}

    for ipt_bin in range(len(plotter.pt_bins_dict[selection]['jms'])):
        pt_bin = plotter.pt_bins_dict[selection]['jms'][ipt_bin]
        plotter.lumi_text_padding = 0.4
        plotter.ratio_plot = True

        for region in plotter.regions[selection]:

            # region_out_dir = output_dir+'/'+pt_bin.replace('to','To')+'_'+region
            # region_out_dir = output_dir+'/'+pt_bin+'_'+region
            region_out_dir = output_dir
            if (not os.path.exists(region_out_dir)):
                os.makedirs(region_out_dir)
            for i in range(len(c_grid['ptbins'])-1):
                for j in range(len(c_grid['etabins'])-1):
                    all_mc = {}
                    all_pseudo = {}
                    for pfflavour in c_grid['pfflavours']:
                        hist_dir = selection+'_%s__mjet_'+pt_bin+region
                        h_pseudo, h_mc = plotter.get_hists(f_hists, plotter.mc_samples[selection],
                                                           hist_dir, pseudo_data=True)
                        varied_hists = {'up': (), 'down': ()}
                        for direction in ['up', 'down']:
                            category = '%i_%i_%s' % (i, j, pfflavour)
                            hist_dir_varied = selection+'_%s__mjet_'+category+'_'+pt_bin+region+'__'+direction
                            varied_hists[direction] = plotter.get_hists(f_hists, plotter.mc_samples[selection],
                                                                        hist_dir_varied, pseudo_data=True)

                        for sample, hist in h_mc.items():
                            if (sample not in all_mc.keys()):
                                all_mc[sample] = {}
                            variations = [varied_hists[direction][1][sample] for direction in ['up', 'down']]
                            variations[0].SetLineColor(ROOT.kBlue+2)
                            variations[1].SetLineColor(ROOT.kRed+2)
                            add_text = plotter.selection_tex[selection]+'\\' + \
                                plotter.legend_names[sample]+'\\'+plotter.pt_bins_tex_dict[pt_bin]

                            legend_entries = [(hist, 'nominal', 'l')]
                            legend_entries.append((variations[0], '+%.1f%% variation ' % variation, 'l'))
                            legend_entries.append((variations[1], '-%.1f%% variation ' % variation, 'l'))
                            # plotter.plot_data_mc(h_data=hist, h_mc=None,
                            #                      plot_title=sample + ' '+pt_bin + ' '+region+' ' + category,
                            #                      out_dir=region_out_dir, additional_hists=variations,
                            #                      legend_entries=legend_entries,
                            #                      additional_text=add_text+'\\'+category)
                            all_mc[sample][pfflavour] = variations
                            all_mc[sample]['nominal'] = hist

                        variations = [varied_hists[direction][0] for direction in ['up', 'down']]
                        variations[0].SetLineColor(ROOT.kBlue+2)
                        variations[1].SetLineColor(ROOT.kRed+2)
                        all_pseudo[pfflavour] = variations
                        all_pseudo['nominal'] = h_pseudo

                        legend_entries = [(h_pseudo, 'nominal', 'l')]
                        legend_entries.append((variations[0], '+%.1f%% variation ' % variation, 'l'))
                        legend_entries.append((variations[1], '-%.1f%% variation ' % variation, 'l'))
                        # plotter.plot_data_mc(h_data=h_pseudo, h_mc=None,
                        #                      plot_title='PseudoData' + ' '+pt_bin.replace('to', 'To')+' ' /
                        #                      + region+' ' + category,
                        #                      out_dir=region_out_dir,
                        #                      additional_hists=variations, legend_entries=legend_entries,
                        #                      additional_text=plotter.selection_tex[selection]+'\\PseudoData\\' /
                        #                      +plotter.pt_bins_tex_dict[pt_bin]+'\\'+category)

                    for sample, hists in all_mc.items():
                        add_text = plotter.selection_tex[selection]+'\\'+plotter.legend_names[sample]+'\\' + \
                            plotter.pt_bins_tex_dict[pt_bin]

                        nominal = hists['nominal']
                        nominal.SetFillStyle(0)
                        nominal.SetLineWidth(2)
                        frame = nominal.Clone()
                        frame.Reset()
                        legend_entries = [(nominal, 'nominal', 'l')]
                        all_variations = []
                        variation_diffs = []
                        up_var = ROOT.THStack()
                        up_var.Add(nominal, 'Hist')
                        down_var = ROOT.THStack()
                        down_var.Add(nominal, 'Hist')
                        # print(hists)
                        max_y = 0.
                        for pfflavour in c_grid['pfflavours']:
                            up = hists[pfflavour][0]
                            up.SetLineColorAlpha(colors[pfflavour], .99)
                            up.SetLineStyle(1)
                            up_var_pf = up.Clone()
                            up_var_pf.Add(nominal, -1.0)
                            max_y = max(max_y, max(abs(up_var_pf.GetMaximum()), abs(up_var_pf.GetMinimum())))
                            up_var.Add(up_var_pf, 'Hist')
                            down = hists[pfflavour][1]
                            down.SetLineColorAlpha(colors[pfflavour], .99)
                            down.SetLineStyle(2)
                            down_var_pf = down.Clone()
                            down_var_pf.Add(nominal, -1.0)
                            max_y = max(max_y, max(abs(down_var_pf.GetMaximum()), abs(down_var_pf.GetMinimum())))
                            down_var.Add(down_var_pf, 'Hist')
                            up.SetLineWidth(2)
                            down.SetLineWidth(2)
                            legend_entries.append((up, '+%.1f%% variation ' % variation+pfflavour, 'l'))
                            legend_entries.append((down, '-%.1f%% variation ' % variation+pfflavour, 'l'))
                            # up.SetLineWidth(1)
                            # down.SetLineWidth(1)
                            all_variations.append(up)
                            all_variations.append(down)
                            variation_diffs.append(up_var_pf)
                            variation_diffs.append(down_var_pf)
                        plotter.ratio_plot = False
                        # plotter.obs_draw_option = 'HIST'
                        plotter.draw_option = 'HIST'
                        plotter.lumi_text_padding = 0.1
                        legend = ROOT.TLegend(*plotter.legend_bbox)
                        legend.SetFillStyle(0)
                        for entry in reversed(legend_entries):
                            legend.AddEntry(*entry)
                        c, _, out_name = plotter.plot_data_mc(h_data=hists['nominal'],  h_mc=None,
                                                              plot_title=selection + '_'+sample + ' '+pt_bin +
                                                              ' '+region[1:]+' all',
                                                              out_dir=region_out_dir,
                                                              additional_hists=all_variations,
                                                              legend_entries=legend_entries,
                                                              additional_text=add_text)
                        legend.Draw('SAME')
                        # hists['nominal'].GetYaxis().SetRangeUser(0,43.5)
                        hists['nominal'].GetXaxis().SetRangeUser(0, 200)
                        c.SaveAs(region_out_dir+'/'+out_name)


def build_asymmerrors(hists):
    nbins = hists[0].GetNbinsX()
    # graph = ROOT.TGraphAsymmErrors(nbins)
    # for i in range(nbins):
    #     nominal = hists[0].GetBinContent(i+1)
    #     low = hists[1].GetBinContent(i+1)
    #     high = hists[2].GetBinContent(i+1)
    #     graph.SetPoint(i,hists[0].GetBinCenter(i+1),nominal)
    #     graph.SetPointError(i,0,0,abs(nominal-low),abs(high-nominal))
    # return graph
    # print(nbins)
    nominal = hists[0]
    low = hists[2]
    high = hists[1]
    graph = ROOT.TGraph(2*nbins)
    for i in range(0, nbins):
        ibin = i+1
        # print('p',i,'bin',ibin)
        # print(i,nominal.GetBinCenter(ibin),high.GetBinContent(ibin))
        graph.SetPoint(i, nominal.GetBinCenter(ibin), high.GetBinContent(ibin))
        # print(nbins+i,nbins-ibin,nominal.GetBinCenter(nbins-ibin+1),low.GetBinContent(nbins-ibin+1))
        graph.SetPoint(nbins+i, nominal.GetBinCenter(nbins-ibin+1), low.GetBinContent(nbins-ibin+1))
    return graph


def plot_variation_envelope(c_name,
                            c_grid,
                            file_name: str = '../macros/Histograms.root',
                            selection: str = 'top',
                            output_dir: str = 'template_variations_envelope'):
    f_hists = ROOT.TFile(file_name, 'READ')
    if (not os.path.exists(output_dir)):
        os.makedirs(output_dir)
    colors = {"chargedH": ROOT.kGreen+2, "neutralH": ROOT.kRed+2, "gamma": ROOT.kOrange, "other": ROOT.kBlue+2}

    plotter.obs_draw_option = 'Hist'
    plotter.ratio_plot = False
    # plotter.ratio_draw_options='H'
    # plotter.ratio_hist_yTitle="#frac{variation}{nominal}"
    plotter.extra_text = "Simulation"
    plotter.obs_draw_option = "H"
    plotter.obs_marker_size = 0
    plotter.draw_option = "F"
    fill_style = 4050
    fill_color_alpha = 0.35

    for ipt_bin in range(len(plotter.pt_bins_dict[selection])):
        pt_bin = plotter.pt_bins_dict[selection][ipt_bin]
        for region in ['_pass', '_passW', '_fail'] if selection == "top" else ['_pass', '_fail']:

            # region_out_dir = output_dir+'/'+pt_bin.replace('to','To')+'_'+region
            region_out_dir = output_dir
            if (not os.path.exists(region_out_dir)):
                os.makedirs(region_out_dir)

            for i in range(len(c_grid['ptbins'])-1):
                for j in range(len(c_grid['etabins'])-1):
                    hist_dir = selection+'_%s__mjet_'+pt_bin+region
                    h_pseudo, h_mc = plotter.get_hists(f_hists, plotter.mc_samples[selection],
                                                       hist_dir, pseudo_data=True)
                    variation_bands = {}
                    for pfflavour in c_grid['pfflavours']:
                        varied_hists = {'up': (), 'down': ()}
                        for direction in ['up', 'down']:
                            category = '%i_%i_%s' % (i, j, pfflavour)
                            hist_dir_varied = selection+'_%s__mjet_'+category+'_'+pt_bin+region+'__'+direction
                            varied_hists[direction] = plotter.get_hists(f_hists, plotter.mc_samples[selection],
                                                                        hist_dir_varied, pseudo_data=True)

                        for sample, hist in h_mc.items():

                            variations = [varied_hists[direction][1][sample] for direction in ['up', 'down']]
                            graph = build_asymmerrors([hist]+variations)

                            graph.SetFillColorAlpha(colors[pfflavour], fill_color_alpha)
                            graph.SetLineColor(colors[pfflavour])
                            graph.SetFillStyle(fill_style)
                            if (sample not in variation_bands.keys()):
                                variation_bands[sample] = {}
                            variation_bands[sample][pfflavour] = graph

                        variations = [varied_hists[direction][0] for direction in ['up', 'down']]
                        graph = build_asymmerrors([h_pseudo]+variations)
                        graph.SetFillColorAlpha(colors[pfflavour], fill_color_alpha)
                        graph.SetLineColor(colors[pfflavour])
                        graph.SetFillStyle(fill_style)
                        if ('pseudo' not in variation_bands.keys()):
                            variation_bands['pseudo'] = {}
                        variation_bands['pseudo'][pfflavour] = graph

                    for sample, hist in h_mc.items():
                        legend_entries = [(hist, 'nominal (%s)' % plotter.legend_names[sample], 'l')]
                        for flavour in c_grid['pfflavours']:
                            h = variation_bands[sample][flavour]

                            legend_entries.append((h, '#pm 1%% variation %s' % flavour, 'f'))
                        plotter.plot_data_mc(h_data=hist, h_mc=None,
                                             plot_title=selection+' '+sample + ' '+pt_bin+' '+region,
                                             out_dir=region_out_dir,
                                             additional_hists=variation_bands[sample].values(),
                                             legend_entries=list(reversed(legend_entries)),
                                             additional_text=plotter.selection_tex[selection]+'\\' +
                                             plotter.legend_names[sample]+'\\'+plotter.pt_bins_tex_dict[pt_bin]
                                             )

                    legend_entries = [(h_pseudo, 'nominal (PseudoData)', 'l')]
                    for flavour in c_grid['pfflavours']:
                        h = variation_bands['pseudo'][flavour]
                        legend_entries.append((h, '#pm 1%% variation %s' % flavour, 'f'))
                    plotter.plot_data_mc(h_data=h_pseudo, h_mc=None,
                                         plot_title=selection+' '+'PseudoData' + ' ' +
                                         pt_bin.replace('to', 'To')+' '+region,
                                         out_dir=region_out_dir,
                                         additional_hists=variation_bands['pseudo'].values(),
                                         legend_entries=list(reversed(legend_entries)),
                                         additional_text=plotter.selection_tex[selection]+'\\'+'PseudoData\\' +
                                         plotter.pt_bins_tex_dict[pt_bin]
                                         )


if (__name__ == '__main__'):

    # selection = 'top'

    category_grids = {
        "oneScale": {
            "ptbins": [0., 100000.],
            "etabins": [0., 9.],
            "pfflavours": ["all"]
        },
        # "PF_flavours":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # }
    }

    # output_dir = "../Plots/template_variations_passW"
    # output_dir = "../Plots/template_variations_PF_flavours"
    output_dir = "../Plots/template_variations_closure_0p005Variation_newPtBins"

    # file_name='../macros/Histograms_PF_flavours_noJEC.root'
    # file_name='../macros/0p01Variation/Histograms_oneScale_JER.root'
    file_name = '../macros/Histograms_oneScale_JER_0p005.root'
    variation = 0.5
    for selection in ['top', 'W']:
        for c_name, c_grid in category_grids.items():
            # plot_template_variations(c_name,c_grid,file_name,selection,'test_variation_newstyle')
            plot_template_variations(c_name, c_grid, file_name, selection, output_dir, variation=variation)
            # plot_variation_envelope(c_name,c_grid,file_name,selection,'test_envelope')
