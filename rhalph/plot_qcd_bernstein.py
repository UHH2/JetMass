import ROOT, argparse, sys, json, os, math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom

rho_min = -6.0
rho_max = -2.1
pt_min = 200
pt_max = 1200
disc_min = 0.12
disc_max = 0.3
left_margin = 0.14
right_margin = 0.16
top_margin = 0.08
bottom_margin = 0.12

def bernstein(x,lamb, nu):
    return binom(nu,lamb) * (x**lamb) * (1-x)**(nu-lamb)

def bernstein_poly(pt, rho, n_pt, n_rho, parameters):
    f = 0
    for i_pt in range(n_pt + 1):
        for i_rho in range(n_rho + 1):
            f += parameters[i_pt][i_rho] * bernstein(rho, i_rho, n_rho) * bernstein(pt, i_pt, n_pt)
    return f
            
def extract_fit_pars(file_path):
    fit_diagnostics = ROOT.TFile(file_path)
    fitargs = fit_diagnostics.Get('fit_s').floatParsFinal()
    bernstein_pars=[]
    for i in range(fitargs.getSize()):
        if('tf' in fitargs.at(i).GetTitle()):
            bernstein_pars.append(fitargs.at(i).getVal())
    return np.array(bernstein_pars).reshape(order[0]+1,order[1]+1)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("config", type=str, help="path to json with config")
#     parser.add_argument("--order", type=int ,nargs='+', default=[3,3], help="tuple containing (pt_degree,rho_degree) for bernstein-polynom")
#     args = parser.parse_args()
    # order = (args.order[0],args.order[1])
    # try:
    #     config = json.load(open(args.config))
    # except IndexError:
    #     print("You must specify a configuration JSON!")
    #     sys.exit(0)
def plot_fit_result(config={'ModelName':'WMassModel'}, order=(3,3)):
    
    bernstein_pars = extract_fit_pars(config['ModelName']+'/fitDiagnostics.root')
    out_dir = config['ModelName']+'/plots/'
    msd_min = 50.
    msd_max = 200.
    
    pt_min = 500.
    pt_max = 1200.

    rho_min = -6
    rho_max = -1

    pt_bins = np.array([500,550,600,675,800,1200])
    min_msd, max_msd = (50,170)
    binwidth = 4
    nbins = int(np.floor((max_msd-min_msd)/binwidth))
    msd_bins = np.linspace(min_msd, nbins*binwidth+min_msd, nbins+1)

    pt_pts, msd_pts = np.meshgrid(pt_bins[:-1]+0.3*np.diff(pt_bins),msd_bins[:-1]+0.5*np.diff(msd_bins),indexing = 'ij')

    rho_pts = 2*np.log(msd_pts/pt_pts)

    pt_scaled = (pt_pts-pt_min)/(pt_max - pt_min)
    rho_scaled = (rho_pts-rho_min)/(rho_max - rho_min)

    valid_rho_bins = (rho_scaled <=1) & (rho_scaled >=0)
    rho_scaled[~valid_rho_bins] = 0

    def eval_bernstein(pt,rho):
        return bernstein_poly(pt, rho, order[1], order[0], bernstein_pars)
    eval_bernstein_map = np.vectorize(eval_bernstein)
    bernstein_map = eval_bernstein_map(pt_scaled.flatten(),rho_scaled.flatten()).reshape(pt_scaled.shape)

    plt.switch_backend('agg')
    plt.rc('text', usetex=True)

    figure, ax = plt.subplots()
    ax.set_title('Bernstein map (%i,%i)'%(order[0],order[1]))
    ax.set_xlabel(r'$m_{SD}$ [GeV]')
    ax.set_ylabel(r'$p_{T}$ [GeV]')
    levels = np.linspace(np.min(bernstein_map)-0.01, np.max(bernstein_map)+0.01, 500)
    plt_cont = ax.contourf(msd_pts, pt_pts, bernstein_map, levels=levels, cmap="rainbow")
    z_bar = figure.colorbar(plt_cont,format='%.2f')
    z_bar.set_label(r'$R_{p/f}(\rho,p_T)$',rotation=270,labelpad=15)
    for f_type in ['.pdf','.png']:
        plt.savefig(out_dir+'bernstein%s'%f_type)
