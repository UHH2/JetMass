#!/nfs/dust/cms/user/albrechs/python/coffea/bin/python3
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from coffea.nanoevents.methods import vector
import hist
import os
from math import ceil
from scipy.optimize import curve_fit
from scipy.stats import moyal
import logging
from typing import Union, Callable
import mplhep as hep
import matplotlib as mpl
hep.style.use("CMS")
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[
    '#e66101',
    '#fdb863',
    '#5e3c99',
    '#b2abd2',
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
])
logger = logging.getLogger(__name__)


def cms_label(ax, fs=20):
    hep.cms.label(label=', Work in Progress', year=2017, ax=ax, fontsize=fs)


def fax(w=9, h=9):
    return plt.subplots(figsize=(w, h))


def calc_mean(x, y):
    return sum(x * y) / sum(y)


def calc_sigma(x, y):
    mean = calc_mean(x, y)
    return np.sqrt(sum(y * (x - mean) ** 2) / sum(y))


def gauss(x, p0, p1, x0, sigma):
    return p0 + p1 * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def landau(x, H, A, x0, sigma):
    return H+A*moyal.pdf(x, x0, sigma)


def landau_gauss(x, p0, p1, x0, sigma):
    # return moyal.pdf(x,x0,sigma)*gauss(x,H,A,x0,sigma)
    return p0+p1*np.exp(-(x-x0)**2 / (2 * sigma**2))*moyal.pdf(x, x0, sigma)


def find_bin(arr, x):
    return np.argmin(abs(arr-x))


def iterative_fit(h: hist.Hist,
                  iterations: int = 2,
                  fit_func: Callable = gauss,
                  plot: Union[plt.Axes, bool] = False
                  ):
    width_multiplier = 1
    edges = h.axes[0].edges
    x = h.axes[0].centers
    y = h.values()
    mean = calc_mean(x, y)
    sigma = calc_sigma(x, y)

    ax = plot  # plot is either True/False or mpl Axes object
    if ax:
        if isinstance(plot, plt.Axes):
            ax = plot
        else:
            f, ax = fax()
        hep.histplot((y, edges), label='init', ax=ax)

    fit_results = []
    for i in range(iterations):
        sigma = abs(sigma)
        xmin = x[0]
        xmax = x[-1]
        if (i > 0):
            xmin = mean-sigma*width_multiplier
            xmax = mean+sigma*width_multiplier
            min_bin = find_bin(x, xmin)
            max_bin = find_bin(x, xmax)
            x = x[min_bin:max_bin]
            y = y[min_bin:max_bin]
            edges = edges[min_bin:max_bin+1]
        # print(sigma)
        if len(y) < 4:
            logger.log(0, "Number of bins got reduced to less than 4." +
                       " skipping fit-iteration."
                       )
            continue
        try:
            popt, pcov = curve_fit(fit_func, x, y,
                                   p0=[min(y),
                                       max(y),
                                       calc_mean(x, y),
                                       calc_sigma(x, y)])
        except BaseException(f"Fit iteration #{i} failed.") as e:
            logger.exception(e)

            fit_results.append({
                'popt': [np.nan, np.nan, np.nan, np.nan],
                'pcov': np.nan,
                'xlim': (xmin, xmax),
            })
            continue

        mean = popt[2]
        sigma = popt[3]
        fit_results.append({
            'popt': popt,
            'pcov': pcov,
            'xlim': (xmin, xmax),
        })
        if plot:
            hep.histplot((y, edges), linestyle='--', label=f'crop {i}', ax=ax)
            x_fine = np.linspace(xmax, xmin, 1000)
            ax.plot(x_fine, fit_func(x_fine, *popt),
                    label=f'fit {i} mean={popt[2]:.2f}')
    if plot:
        ax.set_xlim(*fit_results[0]['xlim'])
    return fit_results, ax


class JMSExtractor(object):

    def __init__(self, tree_file):
        self.tree = tree_file
        self.fs = 18

    @property
    def tree(self):
        return self._tree

    @tree.setter
    def tree(self, tree_file):
        # getting raw tree without High-Level behaviour
        raw_tree = ak.from_parquet(tree_file)

        raw_tree.Jets = ak.zip(
            {
                "pt": raw_tree.Jets.pt,
                "eta": raw_tree.Jets.eta,
                "phi": raw_tree.Jets.phi,
                "mass": raw_tree.Jets.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )

        self._tree = raw_tree

    @tree.deleter
    def tree(self):
        del self._tree

    @property
    def pt_binning(self):
        return self._pt_binning

    @pt_binning.setter
    def pt_binning(self, edges):
        self._pt_binning = edges

    @pt_binning.deleter
    def pt_binning(self):
        del self._pt_binning

    def pt_reco_tex_str(self, i):
        return r'$%.0f < p_{T,\mathrm{reco}} \leq %.0f$ GeV' % (self.pt_binning[:-1][i], self.pt_binning[1:][i])

    def construct_hists(self, groomed=True, JEC=True):
        m_min = 30
        m_max = 200
        n_bins = 2*(m_max-m_min)

        h = hist.Hist(
            hist.axis.Variable(self.pt_binning, name='pt_reco'),
            hist.axis.Regular(n_bins, m_min, m_max, name='msd_reco'),
            hist.axis.Regular(n_bins, m_min, m_max, name='msd_gen'),
            storage=hist.storage.Weight()
        )

        mreco = self.tree.mjet if groomed else self.tree.Jets.mass[:, 0]
        if JEC:
            mreco = mreco*self.tree.jecfactor[:, 0]

        h.fill(
            pt_reco=self.tree.Jets.pt[:, 0]*self.tree.jecfactor[:, 0],
            msd_reco=mreco,
            msd_gen=self.tree.msd_gen_ak8,
            weight=self.tree.weight
        )

        h_response = hist.Hist(
            hist.axis.Variable(self.pt_binning, name='pt_reco'),
            hist.axis.Regular(300, 0, 2, name='msd_response'),
            storage=hist.storage.Weight()
        )

        m_over_30 = (mreco > m_min) & (self.tree.msd_gen_ak8 > m_min) \
            & (mreco < m_max) & (self.tree.msd_gen_ak8 < m_max)

        h_response.fill(
            pt_reco=(self.tree.Jets.pt[:, 0] *
                     self.tree.jecfactor[:, 0])[m_over_30],
            msd_response=(mreco/self.tree.msd_gen_ak8)[m_over_30],
            weight=self.tree.weight[m_over_30]
        )
        return h, h_response

    def create_hist_dict(self):
        self.hists = {f'{grooming}_{jec}':
                      self.construct_hists(groomed=(grooming == 'g'),
                                           JEC=(jec == 'jec'))
                      for grooming in ['g', 'u']
                      for jec in ['jec', 'nojec']}

    def save_control_plots(self, out_directory='.'):
        if not os.path.isdir(out_directory):
            os.makedirs(out_directory)

        groomed_jec_combinations = {
            'g_jec': 'groomed (JEC)',
            'g_nojec': 'groomed (no JEC)',
            'u_jec': 'ungroomed (JEC)',
            'u_nojec': 'ungroomed (no JEC)',
        }

        # plot msd responses
        f, ax = fax()
        pt_id = {'pt_reco': sum}
        for hid, hlabel in groomed_jec_combinations.items():
            hep.histplot(self.hists[hid][1][pt_id],
                         ax=ax, label=hlabel,
                         density=True)

        ax.legend(loc='upper left', fontsize=self.fs)
        cms_label(ax)
        ax.set_xlabel(r"$m_{SD,\mathrm{reco}}/m_{SD,\mathrm{gen}}$")
        ax.set_ylabel(r"$\Delta N / N$")
        f.savefig(f'{out_directory}/msd_responses.pdf', bbox_inches='tight')

        # plot msd distributions
        f, ax = fax()
        for hid, hlabel in groomed_jec_combinations.items():
            hep.histplot(self.hists[hid][0][{**pt_id, 'msd_gen': sum}],
                         ax=ax, label='reco'+hlabel)
        hep.histplot(self.hists['g_jec'][0][{**pt_id, 'msd_reco': sum}],
                     ax=ax, label='gen')
        ax.legend(loc='upper left', fontsize=self.fs)
        cms_label(ax)
        ax.set_xlabel(r"$m_{SD}$ [GeV]")
        ax.set_ylabel(r"Events/$\Delta m_{SD}$")
        f.savefig(f'{out_directory}/msd_distributions.pdf',
                  bbox_inches='tight')

    def extract_jms(self):

        groomings = ['u', 'g']
        jecs = ['jec', 'nojec']
        fit_variables = ['mean_reco', 'mean_gen', 'response']
        xlabels = {
            'mean_reco': r'$m_{SD,\mathrm{reco}}$',
            'mean_gen': r'$m_{SD,\mathrm{gen}}$',
            'response': r'$m_{SD,\mathrm{reco}}/m_{SD,\mathrm{gen}}$',
        }
        fitfunc_dict = {'means_reco': gauss,
                        'means_gen': gauss,
                        'response': gauss
                        }

        self.jms_from_mc = {}
        iterations = 2
        for grooming in groomings:
            for jec in jecs:
                h_ = self.hists[f"{grooming}_{jec}"]
                n_pt_bins = h_[0].axes['pt_reco'].size
                fit_results = {var: [] for var in fit_variables}
                for fit_variable in fit_variables:
                    f = plt.figure(figsize=(20, 12))
                    columns = 3
                    rows = ceil(n_pt_bins/columns)
                    grid = f.add_gridspec(rows, columns, hspace=0.3)
                    grid_overlap = f.add_gridspec(1, 1)
                    cms_label_ax = f.add_subplot(grid_overlap[0])
                    cms_label_ax.axis("off")
                    cms_label(cms_label_ax)

                    for ipt in range(n_pt_bins):
                        ax_ipt = f.add_subplot(grid[ipt])

                        if (fit_variable == 'mean_reco'):
                            fit_result = iterative_fit(
                                h_[0][{'pt_reco': ipt,
                                       'msd_gen': sum}],
                                iterations,
                                fit_func=fitfunc_dict['means_reco'],
                                plot=ax_ipt
                            )
                        elif (fit_variable == 'mean_gen'):
                            fit_result = iterative_fit(
                                h_[0][{'pt_reco': ipt, 'msd_reco': sum}],
                                iterations, fit_func=fitfunc_dict['means_gen'],
                                plot=ax_ipt
                            )
                        elif (fit_variable == 'response'):
                            fit_result = iterative_fit(
                                h_[1][{'pt_reco': ipt}],
                                iterations,
                                fit_func=fitfunc_dict['response'],
                                plot=ax_ipt)
                        ax_ipt.text(ax_ipt.get_xlim()[1]*0.6, ax_ipt.get_ylim()[1]*0.8,
                                    self.pt_reco_tex_str(ipt), fontsize=self.fs-4)
                        ax_ipt.set_xlabel(xlabels[fit_variable])
                        ax_ipt.legend(loc='upper left', fontsize=self.fs-3)

                        fit_results[fit_variable].append(
                            fit_result[0][-1]['popt'][-2] if len(fit_result[0]) > 0 else -999.
                        )

                    f.savefig(f"{grooming}_{jec}_{fit_variable}_control_plots.pdf", bbox_inches="tight")
                self.jms_from_mc[f"response_{grooming}_{jec}"] = \
                    np.array(fit_results['response'])
                self.jms_from_mc[f"means_{grooming}_{jec}"] = \
                    np.array(fit_results['response']) / np.array(fit_results['response'])
        return


if __name__ == '__main__':

    print('loading WJets Tree')
    wjets_jms_extractor = JMSExtractor("WJetsToQQ_tinyTree.parquet")

    wjets_jms_extractor.pt_binning = np.array([500, 550, 650, 725, 800, 1000, 1200])
    print('set pt binning to', wjets_jms_extractor.pt_binning)
    print('constructing and filling hists')
    wjets_jms_extractor.create_hist_dict()
    print('making control plots')
    # wjets_jms_extractor.save_control_plots('test')
    wjets_jms_extractor.extract_jms()
    q = __import__("functools").partial(__import__("os")._exit, 0)  # FIXME
    __import__("IPython").embed()  # FIXME
