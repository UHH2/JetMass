#!/usr/bin/env pythonJMS.sh
import os
import json
import uproot


class FitResults:
    def __init__(self, config):
        self._config = config
        self.name = self._config["ModelName"]
        self.fit_result = json.load(open(f"{self.name}/{self.name}fitResult.json", "r"))

        shapes_file_path = self.name + "/fit_shapes.root"
        if os.path.isfile(shapes_file_path):
            self.fit_shapes_uproot = uproot.open(shapes_file_path)
        else:
            self.fit_shapes = None
        self.results = {}

    def get_massScales(self):
        massScales = {}
        print("getting massscales")
        if self.fit_result is None:
            print("No fit result. skipping")
            return {}
        massScales_parameters = [
            [name, *self.fit_result[name]] for name in self.fit_result.keys() if "massScale" in name
        ]
        print(massScales_parameters)
        for param in massScales_parameters:
            pt_edges = None
            for channel_name, channel in self._config["channels"].items():
                pt_bin = channel_name.lower()
                pt_bin = pt_bin.split("pt")[-1]
                if pt_bin in param[0]:
                    pt_edges = channel["pt_bin"].split("to")
            factor = 0.5 if "Unscaled" in self.name else 1.0
            central = (100 + param[1] * factor) / 100
            error_up = (param[2] * factor) / 100
            error_down = (param[3] * factor) / 100
            massScales.update({param[0]: {"edges": pt_edges, "vals": [central, error_up, error_down]}})
        self.results.update({"jms": massScales})
        return massScales

    def get_peak_positions_new(self):
        import numpy as np
        from scipy.optimize import curve_fit

        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt
        import mplhep as hep

        plt.style.use(hep.style.CMS)

        font_size = 20
        mpl.rcParams["axes.labelsize"] = font_size
        mpl.rcParams["axes.labelsize"] = font_size

        def iterative_gaussian_fit(th1, name="", outDir="./", title=""):
            hist, bin_edges = th1.to_numpy()
            hist = hist / hist.sum()
            bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

            mean = np.average(bin_centres, weights=hist)

            def gauss(x, *p):
                A, mu, sigma = p
                return A * np.exp(-((x - mu) ** 2) / (2.0 * sigma**2))

            p0 = [1.0, mean, 10.0]

            coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
            # sometimes curve_fit returns sigma as negative
            coeff[2] = abs(coeff[2])

            xmin = coeff[1] - 2.5 * coeff[2]
            xmax = coeff[1] + 2.5 * coeff[2]
            bin_min = np.searchsorted(bin_centres, xmin, side="left")
            bin_max = np.searchsorted(bin_centres, xmax, side="right")
            # if number of bins is smaller than 3 increase by 1 in both directions
            # (curve_fit will run into problems otherwise)
            if bin_max - bin_min < 3:
                bin_min -= 1
                bin_max += 1

            coeff1, var_matrix1 = curve_fit(gauss, bin_centres[bin_min:bin_max], hist[bin_min:bin_max], p0=coeff)
            # sometimes curve_fit returns sigma as negative
            coeff1[2] = abs(coeff1[2])

            # getting drawable arrays from fits
            msd = np.linspace(bin_edges[0], bin_edges[-1], 1000)
            hist_fit = gauss(msd, *coeff)
            msd1 = np.linspace(bin_edges[bin_min], bin_edges[bin_max], 1000)
            hist_fit1 = gauss(msd1, *coeff1)

            # setting up canvas
            fig, ax = plt.subplots(figsize=(10, 8))
            hep.cms.text(",work in progress", ax=ax, fontsize=font_size)
            # ax.plot([],[], label=name,color='white')
            for subtitle in title.split("//"):
                ax.plot([], [], label=r"$\mathbf{%s}$" % subtitle, color="white")

            # plotting hist
            hep.histplot((hist, bin_edges), ax=ax)
            # plotting first fit
            ax.plot(msd, hist_fit, label=r"first fit: $\mu = %.2f, \sigma=%.2f$" % (coeff[1], coeff[2]))
            # plotting second fit
            ax.plot(msd1, hist_fit1, label=r"second fit: $\mu = %.2f, \sigma=%.2f$" % (coeff1[1], coeff1[2]))

            plt.xticks(fontsize=font_size - 2)
            plt.yticks(fontsize=font_size - 2)
            ax.set_xlabel(r"$m_{SD}$ [GeV]")
            ax.set_ylabel(r"$\Delta N / N$")
            # ax.legend(fontsize = font_size-2,loc='upper left',bbox_to_anchor=(1, 1))
            ax.legend(fontsize=font_size - 2, loc="upper right")  # ,bbox_to_anchor=(1, 1))

            # print('Fitted mean = ', coeff[1],coeff1[1])
            # print('Fitted standard deviation = ', coeff[2],coeff1[2])
            if not os.path.isdir(outDir):
                os.makedirs(outDir)

            outname = outDir + name + ".pdf"
            plt.savefig(outname)
            return coeff1[1], coeff1[2]

        peak_positions = {}
        for ch in self._config["channels"]:
            selection = self._config["channels"][ch]["selection"]
            regions = self._config["channels"][ch]["regions"]
            regions.remove("fail")
            peak_positions.update({ch: {}})
            for region in regions:
                peak_positions[ch].update({region: {}})
                for sample in self._config["channels"][ch]["signal"]:
                    peak_positions[ch][region].update({sample: {}})
                    for fit in ["prefit", "postfit"]:
                        hist_dir = ch + region + "_" + fit + "/" + sample
                        th1 = self.fit_shapes_uproot[hist_dir]
                        name = self.name + "_" + selection + "_" + hist_dir.replace("/", "_")
                        title = selection + "~" + region + "~" + sample.replace("_", "~") + "//" + ch + " " + fit
                        outDir = "iterative_gaussians/"
                        mu, sigma = iterative_gaussian_fit(th1, name, outDir, title)
                        pt_edges = [float(pt) for pt in self._config["channels"][ch]["pt_bin"].split("to")]
                        peak_positions[ch][region][sample].update({fit: {"values": [mu, sigma], "pt_edges": pt_edges}})
        self.results.update({"peak_positions": peak_positions})


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fits", nargs="+", default=[], help="List of fits (separated by spaces).", required=True)
    parser.add_argument("--fit-peaks", action="store_true")
    parser.add_argument("-o", "--output", default="fitResults.json")
    args = parser.parse_args()
    configs = [(fit_dir if "config.json" else f"{fit_dir}/config.json") for fit_dir in args.fits]

    all_results = {}
    for config_path in configs:
        config = json.load(open(config_path))
        name = config["ModelName"]
        print(name)
        FR = FitResults(config)
        sf = FR.get_massScales()

        if args.fit_peaks:
            FR.get_peak_positions_new()

        all_results.update({name: FR.results})

    open(args.output, "w").write(json.dumps(all_results, indent=2))
