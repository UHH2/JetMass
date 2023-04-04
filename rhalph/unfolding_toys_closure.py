from __future__ import print_function
import os
import json
import ROOT
import hist
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import moyal
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
hep.set_style(hep.style.CMS)


def calc_mean(x, y):
    return sum(x * y) / sum(y)


def calc_sigma(x, y):
    mean = calc_mean(x, y)
    return np.sqrt(sum(y * (x - mean) ** 2) / sum(y))


def gauss(x, p0, p1, x0, sigma):
    return p0 + p1 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def landau(x, H, A, x0, sigma):
    return H + A * moyal.pdf(x, x0, sigma)


def landau_gauss(x, p0, p1, x0, sigma):
    # return moyal.pdf(x,x0,sigma)*gauss(x,H,A,x0,sigma)
    return p0 + p1 * np.exp(-((x - x0) ** 2) / (2 * sigma**2)) * moyal.pdf(
        x, x0, sigma
    )


class UnfoldingToyClosure(object):
    def __init__(self, modeldir, workdir, name="ToyClosure"):
        self._modeldir = os.path.abspath(modeldir)
        self._workdir = os.path.abspath(workdir)
        self._cmssw_dir = (
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/"
            "CMSSW_11_3_4/src"
        )

        self._name = name

        self._config = json.load(open(self._modeldir + "/config.json", "r"))

        self._pois = ["r_" + genbin for genbin in self._config.get("genbins", [""])]

        if self._pois == ["r_"]:
            raise RuntimeError("could not find any genbins in config.")

        self._initial_seed = 123456

        if not os.path.exists(workdir):
            os.makedirs(workdir)

    def submit_toys(self, n_toys):
        os.chdir(self._workdir)
        for i in range(n_toys):
            cmd = (
                "combineTool.py -M FitDiagnostics -d {}/model_combined.root  --cminDefaultMinimizerStrategy 0 ".format(
                    self._modeldir
                )
            )
            cmd += "--redefineSignalPOIs {} ".format(",".join(self._pois))
            cmd += "--setParameters {}=1 ".format("=1,".join(self._pois))
            cmd += "--robustFit 1 "
            cmd += "-n .{}{} ".format(self._name, i)
            cmd += "--seed {} ".format(int(self._initial_seed*10**(np.floor(np.log10(n_toys))+1) + i))
            cmd += "-t 1 "
            cmd += "--job-mode condor "
            cmd += "--task-name {}_job_{} ".format(self._name, i)
            os.system(cmd)

    def get_parameters(self, itoy, fit_result_name="fit_s"):
        fname = "{}/fitDiagnostics.{}{}.root".format(self._workdir, self._name, itoy)
        file_ = ROOT.TFile.Open(fname, "READ")
        try:
            fit_result = file_.Get(fit_result_name)
            fit_result_pars = fit_result.floatParsFinal()
        except BaseException:
            print("could not find fit result ({}) in provided file (itoy={}).".format(fit_result_name, itoy))
            return None

        poi_result_tmp = {}
        for p in fit_result_pars:
            if p.GetName() in self._pois or len(self._pois) == 0:
                poi_result_tmp[p.GetName()] = [p.getVal(), p.getErrorHi(), p.getErrorLo()]
        return [np.array([poi_result_tmp[name][i] for name in self._pois]) for i in [0, 1, 2]]

    def plot_toys(self, n_toys_submitted):
        toy_postfit = np.array([arr for arr in [self.get_parameters(itoy) for itoy in range(n_toys_submitted)] if arr])
        n_toys = len(toy_postfit[:, 0, 0])
        print("extracted fitresults for {}/{} toys".format(n_toys, n_toys_submitted))

        for ipoi, poi in enumerate(self._pois):
            postfit_results = toy_postfit[:, 0, ipoi]

            # xmin = max(0, np.floor(postfit_results.min()))
            # xmax = min(10, np.ceil(postfit_results.max()))
            # nbins = int(np.ceil(np.sqrt(n_toys)))
            xmin = 0.0
            xmax = 2.0
            nbins = 20

            h = hist.Hist.new.Reg(nbins, xmin, xmax, name=poi).Int64()
            h.fill(postfit_results)

            # fit hist with some distribution
            x = h.axes[0].centers
            y = h.values()
            fit_func = gauss

            try:
                popt, pcov = curve_fit(
                    fit_func, x, y, p0=[min(y), max(y), calc_mean(x, y), calc_sigma(x, y)]
                )
            except RuntimeError:
                print("fit failed for poi: {}".format(poi))

            f, ax = plt.subplots()

            hep.label._exp_label(
                exp="",
                llabel="Private work (CMS data/simulation)",
                year=self._config.get("year", "UL17"),
                ax=ax,
                fontsize=20,
            )

            hep.histplot(h, label=r"$N_\mathrm{toys}=%i$ (mean=%.2f, width=%.2f)" % (n_toys, calc_mean(x,y), calc_sigma(x,y)), ax=ax, yerr=False)

            x_fine = np.linspace(xmin, xmax, 1000)
            ax.plot(x_fine, fit_func(x_fine, *popt), label=r"fit ($\mu=%.2f, \sigma=%.2f$)" % (popt[2], popt[3]))

            ax.set_xlabel(poi)
            ax.set_ylabel("toys")
            ax.legend()
            f.savefig("{}/{}_toys.pdf".format(self._modeldir, poi), bbox_inches="tight")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--model", help="model dir containing workspace 'model_combined.root'.")
    parser.add_argument("-w", "--workdir", help="workdir to perform fits in")
    parser.add_argument("-t", "--toys", help="number of toys to produce/fit", type=int)
    parser.add_argument("-s", "--submit", action="store_true", help="submit toys")
    parser.add_argument("-p", "--plot", action="store_true", help="plot finished toys")
    args = parser.parse_args()

    uct = UnfoldingToyClosure(modeldir=args.model, workdir=args.workdir, name="ToyTest")
    if args.submit:
        uct.submit_toys(args.toys)
    if args.plot:
        uct.plot_toys(args.toys)