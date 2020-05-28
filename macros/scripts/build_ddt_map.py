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

def plot_map_root(file_path, map_name, wp):
    from ROOT import gStyle, gROOT, TFile, TCanvas, TF1, TLatex
    import os
    # rho_min = -6.0
    # rho_max = -2.1
    # pt_min = 200
    # pt_max = 1200
    rho_min = -10
    rho_max = 0
    pt_min = 0
    pt_max = 1500
    disc_min = 0.12
    disc_max = 0.3
    left_margin = 0.14
    right_margin = 0.16
    top_margin = 0.08
    bottom_margin = 0.12
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
    gStyle.SetPadLeftMargin(left_margin)
    gStyle.SetPadBottomMargin(bottom_margin)
    gStyle.SetPadTopMargin(top_margin)
    gStyle.SetPadRightMargin(right_margin)
    
    gStyle.SetMarkerSize(0.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetTitleSize(0.05, "XYZ")
    gStyle.SetLabelSize(0.04, "XYZ")
    gStyle.SetNdivisions(506, "XYZ")
    gStyle.SetNumberContours(25)
    gStyle.SetLegendBorderSize(0)
    f = TFile(file_path)
    hist = f.Get(map_name)
  
    hist.GetXaxis().SetTitle("#rho")
    hist.GetYaxis().SetTitle("p_{T} [GeV]")
    if('n2' in map_name.lower()):
        hist.GetZaxis().SetTitle("N2^{DDT} %.1f"%(wp*100)+"% quantile")
    else:
        hist.GetZaxis().SetTitle("DeepBoosted WvsQCD %.1f"%(wp*100)+"% quantile")
        disc_max=1
    hist.GetXaxis().SetRangeUser(rho_min, rho_max)
    hist.GetYaxis().SetRangeUser(pt_min, pt_max)
    if(disc_max>hist.GetMaximum()):
        hist.GetZaxis().SetRangeUser(disc_min, disc_max)

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
    hist.SetTitleSize(18.0)
	
    # isomasses = [20,55,80,120,200]
    # isomasses = [20,65,80,125,200]
    # isomasses = range(40,200,20)
    isomasses = [40,80,120,200]
    str_isomass = "%.2f*TMath::Exp(-x/2)"
    tf1_isomasses = []
    for i in range(len(isomasses)):
        msd = isomasses[i]
        new_isomass = TF1('isomass_%i'%int(msd),str_isomass%msd,rho_min,rho_max)
        # new_isomass.SetLineColor(920+i)
        new_isomass.SetLineColorAlpha(1,0.4)
        tf1_isomasses.append(new_isomass)


    c1 = TCanvas("c1","c1",700,600)
    c1.cd()
    hist.Draw("colz")
    latex = TLatex()
    latex.SetNDC(1)
    latex.SetTextColor(1)
    latex.SetTextFont(43)
    latex.SetTextSize(15.5)
    latex.SetTextAngle(297)
    for i in range(len(isomasses)):
        x_pitch = (1-left_margin-right_margin)/(rho_max-rho_min)
        y_pitch = (1-bottom_margin-top_margin)/(pt_max-pt_min)
        msd = isomasses[i]
        pt=msd*np.exp(-rho_min/2)
        pt = 0.5*(pt_max-pt_min)+pt_min if pt > pt_max else pt
        angle = ((180/np.pi)*np.arctan(2*x_pitch/(pt*y_pitch)))+272
        latex.SetTextAngle(angle)
        tf1_isomass = tf1_isomasses[i]
        tf1_isomass.Draw('SAME')
        x_pos = left_margin+((2*np.log(msd/pt)-rho_min))*x_pitch+0.01
        y_pos = bottom_margin+(pt-pt_min)*y_pitch +0.01
        latex.DrawLatex(x_pos,y_pos,'m_{SD} = %.1f GeV'%msd)
    c1.SaveAs(out_file_path+".pdf")
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

    parser.add_argument('-s','--smooth', action='store_true', help='specify if you want to postprocess the maps with an gaussian filter' )
    parser.add_argument('-g','--gaus_sigma', type=float, default=1.0, help='specify what sigma should be used by the smoothing gaus-kernel.')
    parser.add_argument('--justplots',action='store_true', help='just make some ugly plots. -i acts here as inputfiles containing the finished 2D maps')
    parser.add_argument('-p','--makeplots',action='store_true', help='make some ugly plots')
    
    args = parser.parse_args()

    if(args.hist_dir is None):
        hist_dirs = [None for k in args.hist_names]
    elif(len(args.hist_dir) == len(args.hist_names)):
        hist_dirs = args.hist_dir
    elif(len(args.hist_dir) == 1):
        hist_dirs = [args.hist_dir[0] for k in args.hist_names]
    else:
        raise ValueError('number of provided hist_names (%i) and hist_dirs (%i) are not compatible! You must provide either no hist_dir, one hist_dir in total or one hist_dir per hist_name!'%(len(args.hist_names),len(args.hist_dir)))

    hist_names = zip(hist_dirs,args.hist_names)

    if(args.smooth):
        smooth_suffix = ("_smooth_gaus%.2fsigma"%args.gaus_sigma).replace('.','p')
        args.output = args.output.replace('.root','')+smooth_suffix +".root"
    else:
        smooth_suffix = ''

    if args.justplots:
        hist_names = zip(hist_dirs,args.hist_names)
        for hist_dir,hist_name in hist_names:
            for wp in args.working_point:
                plot_map_root(args.input,hist_name+'_'+str(wp).replace('.','p')+smooth_suffix+'_'+hist_dir)
        exit(0)
    th3_infile = uproot.open(args.input)
   

    if 'list' in args.hist_names:
        print('listing available hists:')
        if(hist_dirs[0] is None):
                print(th3_infile.keys())
        for hist_dir in hist_dirs:
            print(th3_infile[hist_dir].keys())
        exit(0)

        
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
            print("misstag rate:",wp)
            smooth_map = build_ddt_map(th3, float(wp), smooth = args.smooth, gaus_sigma=args.gaus_sigma)
            smooth_map_th2 = numpy_to_TH2(smooth_map,hist_name+'_'+str(wp).replace('.','p')+smooth_suffix+'_'+hist_dir,th3)

            th2_outfile.Write()
    th2_outfile.Close()
    if args.makeplots:
        hist_names = zip(hist_dirs,args.hist_names)
        for hist_dir,hist_name in hist_names:
            for wp in args.working_point:
                plot_map_root(args.output,hist_name+'_'+str(wp).replace('.','p')+smooth_suffix+'_'+hist_dir,wp)
