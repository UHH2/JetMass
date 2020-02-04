import uproot
if tuple(map(int,uproot.__version__.split('.'))) < (3,8):
    raise BaseException(""""You'll need a newer version of uproot (currently %s). 
Try 
pip install --user --upgrade uproot 
or (el6)
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-slc6-gcc8-opt/setup.sh 
or (el7)
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh"""%uproot.__version__)

import ROOT
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters


def numpy_to_TH2(np_arr,name,th3):
    '''
    converts 2D-numpy array with info of original TH3 into a TH2D histogram.
    '''
    np_arr_clean = np.where(np_arr<0,0,np_arr)
    from root_numpy import array2hist
    # The x bins here are hard coded just to make them comparable to annas old plots.
    # th2 = ROOT.TH2D(name,name,90,th3.xlow,-1,th3.ynumbins,th3.ylow,th3.yhigh)
    # array2hist(np_arr_clean[:-10,:],th2)
    th2 = ROOT.TH2D(name,name,th3.xnumbins,th3.xlow,th3.xhigh,th3.ynumbins,th3.ylow,th3.yhigh)
    array2hist(np_arr_clean,th2)
    return th2

def plot_map_root(file_path, map_name):
    from ROOT import gStyle, gROOT, TFile, TCanvas
    import os
    out_dir = os.path.dirname(file_path)
    out_file_path = out_dir+("/" if len(out_dir)>0 else "")+map_name;
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetLegendBorderSize(0)
    
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetTitleOffset(0.86,"X")
    gStyle.SetTitleOffset(1.6,"Y")
    
    gStyle.SetPadLeftMargin(0.14)
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadTopMargin(0.08)
    gStyle.SetPadRightMargin(0.16)
    
    gStyle.SetMarkerSize(0.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetTitleSize(0.05, "XYZ")
    gStyle.SetLabelSize(0.04, "XYZ")
    gStyle.SetNdivisions(506, "XYZ")
    gStyle.SetLegendBorderSize(0)
    f = TFile(file_path)
    hist = f.Get(map_name)
  
    hist.GetXaxis().SetTitle("#rho")
    hist.GetYaxis().SetTitle("p_{T} [GeV]")
    hist.GetZaxis().SetTitle("N2^{DDT} 5/% quantile")
    
    hist.GetXaxis().SetRangeUser(-6.0,-2.1)
    hist.GetYaxis().SetRangeUser(200,1200)
    hist.GetZaxis().SetRangeUser(0.12,0.3)

    Font=43
    TitleSize=24.0
    TitleOffset=1.3
    LabelSize=18.0
    hist.GetYaxis().SetTitleFont(Font)
    hist.GetYaxis().SetTitleSize(TitleSize)
    hist.GetYaxis().SetTitleOffset(TitleOffset)
    hist.GetYaxis().SetLabelFont(Font)
    hist.GetYaxis().SetLabelSize(LabelSize)

    hist.GetXaxis().SetTitleFont(Font)
    hist.GetXaxis().SetTitleSize(TitleSize)
    hist.GetXaxis().SetTitleOffset(TitleOffset)
    hist.GetXaxis().SetLabelFont(Font)
    hist.GetXaxis().SetLabelSize(LabelSize)
    
    hist.GetZaxis().SetTitleFont(Font)
    hist.GetZaxis().SetTitleSize(TitleSize)
    hist.GetZaxis().SetTitleOffset(TitleOffset)
    hist.GetZaxis().SetLabelFont(Font)
    hist.GetZaxis().SetLabelSize(LabelSize)
    
    hist.GetZaxis().SetNdivisions(510)
    
    hist.SetTitle(map_name.replace('_',' '))
    hist.SetTitleFont(43)
    hist.SetTitleSize(20.0)
	
    c1 = TCanvas("c1","c1",700,600)
    c1.cd()
    hist.Draw("colz")
    c1.SaveAs(out_file_path+".pdf")
    c1.SaveAs(out_file_path+".png")
    del c1

def plot_map(x_arr, y_arr, contents_arr, name = 'ddt_map', levels = None):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title(name.replace('_',' '))
    if levels is None:
        levels = np.linspace(0,1,500)
    ax.contour(x_arr, y_arr, contents_arr, levels)
    plt.savefig(name+'.png')
    plt.savefig(name+'.pdf')
    print('saved map as',name,'.[pdf,png]')
    

def build_ddt_map(hqcd, qcd_eff, smooth = True, gaus_sigma=1.0):
    vals_var, edges = hqcd.numpy()
    cum_sum  = np.cumsum(vals_var, axis = 2)
    max_val = cum_sum[:,:,-1]
    norma = cum_sum / np.maximum(1e-10,max_val[:,:,np.newaxis])

    if('high' in args.reject or 'h' in args.reject):
        qcd_eff = 1 - qcd_eff

    res = np.apply_along_axis(lambda norma: norma.searchsorted(qcd_eff), axis = 2, arr = norma)
    mask_zero = np.apply_along_axis(lambda norma: norma.searchsorted(1e-10), axis = 2, arr = norma)
    res = np.where(mask_zero == norma.shape[2], 0, res)

    #when X% quantiles position has not been found set it to last bin of discriminator-hist
    res = np.where(res > norma.shape[2]-1, norma.shape[2]-1, res)

    def bineval(a):
        return edges[0][2][a]
    bin_func = np.vectorize(bineval)
    quantile_map = bin_func(res)
    if(not smooth):
        return quantile_map
    smooth_map = filters.gaussian_filter(quantile_map,gaus_sigma)
    return smooth_map

    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help='<Required> Path to ROOT file containing the TH3 histogram(s), which will be used to build the DDT-map',required = True)
    parser.add_argument('-n','--hist-names', nargs='+', help='<Required> List of names of the TH3-Histogram(s) in specified ROOT-File', required = True)
    parser.add_argument('-o','--output',default = 'ddt_maps.root', help='Path to ROOT file, in which the DDT-map(s) will be saved as TH2-Histogram(s)')
    parser.add_argument('-wp','--working-point',nargs='+',default=[0.05], type=float,help='List of working-point percentile(s), at which the ddt-map(s) should be created')
    
    parser.add_argument('--hist-dir', nargs='+', default=None, help = 'List of the directories, in which the respective histograms specified with -n/--hist-names can be found in the specified input ROOT-File')

    parser.add_argument('-r','--reject',choices=['high','low','l','h'], default='low', help='Specifies what kind of values of the discriminator reject the background. If low/l wp*100 percent quantile is derived. If high/h (1-wp)*100 percent quantile is derived.')
    parser.add_argument('-s','--smooth', action='store_true', help='specify if you want to postprocess the maps with an gaussian filter' )
    parser.add_argument('-g','--gaus_sigma', type=float, default=1.0, help='specify what sigma should be used by the smoothing gaus-kernel.')
    parser.add_argument('-p','--makeplots',action='store_true', help='make some ugly plots')
    
    args = parser.parse_args()
    th3_infile = uproot.open(args.input)
   
    if(args.hist_dir is None):
        hist_dirs = [None for k in args.hist_names]
    elif(len(args.hist_dir) == len(args.hist_names)):
        hist_dirs = args.hist_dir
    elif(len(args.hist_dir) == 1):
        hist_dirs = [args.hist_dir[0] for k in args.hist_names]
    else:
        raise ValueError('number of provided hist_names (%i) and hist_dirs (%i) are not compatible! You must provide either no hist_dir, one hist_dir in total or one hist_dir per hist_name!'%(len(args.hist_names),len(args.hist_dir)))

    hist_names = zip(hist_dirs,args.hist_names)

    if 'list' in args.hist_names:
        print('listing available hists:')
        if(hist_dirs[0] is None):
                print(th3_infile.keys())
        for hist_dir in hist_dirs:
            print(th3_infile[hist_dir].keys())
        exit(0)

    if(args.smooth):
        smooth_suffix = ("_smooth_gaus%.2fsigma"%args.gaus_sigma).replace('.','p')
        args.output = args.output.replace('.root','')+smooth_suffix +".root"
    else:
        smooth_suffix = ''
        
    th2_outfile = ROOT.TFile(args.output,'RECREATE')
    th2_outfile.cd()
    for hist_dir,hist_name in hist_names:
        if hist_dir is None:
            th3 = th3_infile[hist_name]
        else:
            th3 = th3_infile[hist_dir][hist_name]

        _,edges = th3.numpy()
        edges_x = edges[0][0][:-1]
        edges_y = edges[0][1][:-1]
        
        for wp in args.working_point:
            smooth_map = build_ddt_map(th3, float(wp), smooth = args.smooth, gaus_sigma=args.gaus_sigma)
            smooth_map_th2 = numpy_to_TH2(smooth_map,hist_name+'_'+str(wp).replace('.','p')+smooth_suffix+'_'+hist_dir,th3)
            # if args.makeplots:
            #     plot_map(edges_x, edges_y, smooth_map.T,'ddt_map_'+hist_name+'_'+str(wp).replace('.','p'))

            th2_outfile.Write()
    th2_outfile.Close()
    if args.makeplots:
        hist_names = zip(hist_dirs,args.hist_names)
        for hist_dir,hist_name in hist_names:
            for wp in args.working_point:
                plot_map_root(args.output,hist_name+'_'+str(wp).replace('.','p')+smooth_suffix+'_'+hist_dir)
