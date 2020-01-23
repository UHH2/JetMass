import sys
if sys.version_info < (3,6):
    raise BaseException("This script has been tested in python3.6 with special setup including coffea. If you want to try it with different python versions, feel free to disabel this exception at the top of this script!")

import uproot
import ROOT
import numpy as np
import matplotlib.pyplot as plt
# import coffea.hist as hist
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
    

def build_ddt_map(hqcd, qcd_eff, smooth = True):
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
    smooth_map = filters.gaussian_filter(quantile_map,1)
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
            smooth_map = build_ddt_map(th3, float(wp), smooth = True)
            smooth_map_th2 = numpy_to_TH2(smooth_map,hist_name+'ddt_map_smooth_'+str(wp).replace('.','p'),th3)
            
            plot_map(edges_x, edges_y, smooth_map.T,'ddt_map_'+hist_name+'_'+str(wp).replace('.','p'))
            th2_outfile.Write()
    th2_outfile.Close()
