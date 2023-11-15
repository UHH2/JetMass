#!/usr/bin/env python3
import ROOT
from ROOT import gSystem
import numpy as np
import rhalphalib as rl
import json
from CombineUtils import import_config
ROOT.PyConfig.IgnoreCommandLineOptions = True


def correlation_coefficient(fname, pois):
    file_ = ROOT.TFile.Open(fname, "READ")
    fit_s = None
    corr_coeff = np.array([-1])
    try:
        fit_s = file_.Get("fit_s")
        n_pars = fit_s.floatParsFinal().getSize()

        poi_ind = [fit_s.floatParsFinal().index(p) for p in pois]

        cov = np.array([[fit_s.covarianceMatrix()[row][col] for col in range(n_pars)] for row in range(n_pars)])
        cov_inv = np.linalg.inv(cov)
        corr_coeff = np.sqrt(1 - 1 / (cov_inv[poi_ind, poi_ind] * cov[poi_ind, poi_ind]))
        file_.Close()
        del fit_s

        del cov, cov_inv
    except BaseException as e:
        print(e)
    return corr_coeff


def scale_lumi(hist, lumiscale):
    if lumiscale == 1.0:
        return hist
    hist_old = hist.Clone()
    hist.Reset()
    for i in range(1, hist.GetNbinsX() + 1):
        hist.SetBinContent(i, hist_old.GetBinContent(i) * lumiscale)
        hist.SetBinError(i, hist_old.GetBinError(i) * np.sqrt(lumiscale))
    return hist


def build_pseudo(samples, hist_file, hist_dir, vary_pseudo_like, MINIMAL_MODEL=False, seed=123456, msd_bins=None):
    gSystem.ProcessEvents()
    lumiScale = 1.0
    if len(vary_pseudo_like) > 0 and "lumiScale" in vary_pseudo_like[0]:
        lumiScale = float(vary_pseudo_like[0].split(":")[-1])
    print("LumiScale:", lumiScale)
    if len(vary_pseudo_like) > 1:
        print("building pseudo data from ", samples, " (varied like", vary_pseudo_like, ")")
        first_hist_dir = (hist_dir + "__%s") % (samples[0], vary_pseudo_like[1], vary_pseudo_like[2])
        print(first_hist_dir)
        print(hist_file)
        pseudo_data_hist = hist_file.Get(first_hist_dir)
        for sample in samples[1:]:
            this_hist_dir = (hist_dir + "__%s") % (sample, vary_pseudo_like[1], vary_pseudo_like[2])
            pseudo_data_hist.Add(hist_file.Get(this_hist_dir))
        pseudo_data_hist = scale_lumi(pseudo_data_hist, lumiScale)
        return pseudo_data_hist
    else:
        print("building pseudo data from ", samples, " (unvaried)")
        toys = False
        if len(vary_pseudo_like) > 0 and "toys" in vary_pseudo_like[0]:
            toys = True
        first_hist_dir = hist_dir % (samples[0], "")
        pseudo_data_hist = hist_file.Get(str(first_hist_dir))
        # print(pseudo_data_hist)
        for sample in samples[1:]:
            if ("qcd" not in sample.lower() or "WJetsMatched" not in sample.lower()) and MINIMAL_MODEL:
                continue
            this_hist_dir = hist_dir % (sample, "")
            print(this_hist_dir)
            pseudo_data_hist.Add(hist_file.Get(str(this_hist_dir)))

        pseudo_data_hist.GetXaxis().SetTitle("msd")
        pseudo_data_hist.SetName("msd")
        print(msd_bins)
        if msd_bins is not None:
            pseudo_data_hist = pseudo_data_hist.Rebin(len(msd_bins) - 1, "msd", msd_bins)

        if toys:
            N_entries_toys = pseudo_data_hist.Integral() * lumiScale
            toy_pseudo_data = pseudo_data_hist.Clone()
            toy_pseudo_data.Reset()
            n_toys = 100
            for i in range(n_toys):
                ROOT.gRandom.SetSeed(seed + i)
                toy_pseudo_data.FillRandom(pseudo_data_hist, int(N_entries_toys / n_toys))
            toy_pseudo_data.Scale(N_entries_toys / toy_pseudo_data.Integral())
            return toy_pseudo_data
        else:
            # poissonify uncertainty
            for i in range(1, pseudo_data_hist.GetNbinsX() + 1):
                pseudo_data_hist.SetBinError(i, np.sqrt(pseudo_data_hist.GetBinContent(i)))

            pseudo_data_hist = scale_lumi(pseudo_data_hist, lumiScale)
            return pseudo_data_hist


def build_mass_scale_variations(configs, args):
    grid_hist_file_name = configs["gridHistFileName"]
    # print('reading grid for nuisance parameter:')
    grid_hist_file = ROOT.TFile(grid_hist_file_name, "READ")
    grid_hist = grid_hist_file.Get("grid")
    grid_axes = dict(item.strip().split("=") for item in grid_hist.GetTitle().split(","))
    x_bins = range(grid_hist.GetNbinsX())
    y_bins = range(grid_hist.GetNbinsY())

    categories_hist = grid_hist_file.Get("categories")
    particle_categories = []
    for i in range(1, categories_hist.GetNbinsX() + 1):
        particle_categories.append(categories_hist.GetXaxis().GetBinLabel(i))

    grid_hist_file.Close()

    # building dictionary to help setup the massScales in actual workspace creator.
    # dict hold the massScale substrings containing pT-bin info
    # and application info sub-dict, which tells the code on which samples/regions
    # the respective nuisance should act on.
    # and can also be separated into W/Z/top scales.
    #
    # for the simple cases of one scale per pT bin (or inclusive in pT) we want the scale(s) to act on either all
    # or just the signal templates:
    base_application_info = {
        "regions": ["pass", "passW", "fail"],
        "samples": "signal" if args.VaryOnlySignal else "all",
    }
    if args.pTdependetMassScale:
        mass_scale_setup_dict = {
            ("PT" + ch_name.split("Pt")[-1]): base_application_info for ch_name in configs["channels"].keys()
        }
    else:
        mass_scale_setup_dict = {"inclusive": base_application_info}

    # get set of selections to decide which massScale names to return
    selections = set([configs["channels"][c]["selection"] for c in configs["channels"].keys()])

    if args.separateMassScales:
        separate_mass_scale_setup_dict = {}
        for mass_scale_suffix in mass_scale_setup_dict.keys():
            separate_mass_scale_setup_dict.update(
                {
                    "W_"
                    + mass_scale_suffix: {
                        "regions": ["pass", "passW", "fail"],
                        "samples": ["WJetsMatched", "TTToSemiLeptonic_mergedW"],
                    }
                }
            )
            # if any(s in selections for s in ["W", "Zbb"]):
            if "Zbb" in selections:
                separate_mass_scale_setup_dict.update(
                    {"Z_" + mass_scale_suffix: {"regions": ["pass", "passW", "fail"], "samples": ["ZJetsMatched"]}}
                )
            if "top" in selections:
                separate_mass_scale_setup_dict.update(
                    {
                        "top_"
                        + mass_scale_suffix: {
                            "regions": ["pass", "passW", "fail"],
                            "samples": ["TTToSemiLeptonic_mergedTop"],
                        }
                    }
                )
        mass_scale_setup_dict = separate_mass_scale_setup_dict

    # setting up nuisances correspondig to consituent-variation according to categories from grid
    grid_nuisances = []
    mass_scale_names = []

    for category in particle_categories:
        for x_bin in x_bins:
            for y_bin in y_bins:
                mass_scale_names = [
                    "massScale_%s%i_%s%i_%s_%s"
                    % (grid_axes["x"], x_bin, grid_axes["y"], y_bin, category, mass_scale_suffix)
                    for mass_scale_suffix in mass_scale_setup_dict.keys()
                ]
                grid_nuisances.append(
                    [
                        {
                            mass_scale_setup_dict.keys()[i]: {
                                "application_info": list(mass_scale_setup_dict.values())[i],
                                "nuisance": rl.NuisanceParameter(mass_scale_names[i], "shape", 0, -0.5, 0.5),
                            }
                            for i in range(len(mass_scale_names))
                        },
                        x_bin,
                        y_bin,
                        category,
                    ]
                )
    return grid_nuisances, mass_scale_names


def extract_fit_results(configs, return_result=False):
    '''
    extracting postfit parameters from fitDiagnostics.root
    '''
    if isinstance(configs, str):
        configs = import_config(configs)
    model_dir = configs.get("ModelDir", configs.get("ModelName"))
    fit_succeded = True
    try:
        fit_diagnostics = ROOT.TFile(model_dir + "/fitDiagnostics.root", "READ")
        fit_result = fit_diagnostics.Get("fit_s")

        fit_result_parameters = {}
        for p in fit_result.floatParsFinal():
            value = p.getVal()
            err_hi = p.getErrorHi()
            err_low = p.getErrorLo()
            err = p.getError()
            if err_hi == 0 or err_low == 0:
                err_low = -1.0 * err
                err_hi = err
            fit_result_parameters[p.GetName()] = [value, err_hi, err_low]
        fit_succeded = fit_result.status() <= 3
        if return_result and fit_succeded:
            return fit_result_parameters
        open(model_dir + "/" + configs["ModelName"] + "fitResult.json", "w").write(
            json.dumps(fit_result_parameters, sort_keys=True, indent=2)
        )
        return fit_succeded
    except BaseException as e:
        print("fit failed. only plotting prefit distributions from fitDiangnostics (beware weird shape uncertainties)")
        print(e)
        fit_succeded = False
    return fit_succeded


def extract_correlation_matrix(configs):
    if isinstance(configs, str):
        configs = import_config(configs)
    model_dir = configs.get("ModelDir", configs.get("ModelName"))
    file_ = ROOT.TFile.Open(model_dir + "/fitDiagnostics.root", "READ")
    fit_s = None
    pois = ["r_{}".format(genbin) for genbin in configs.get("genbins", [])]
    if pois == []:
        genbininfo = configs.get("unfolding_bins", {})
        msdgen_edges = genbininfo.get("msdgen", [])
        ptgen_edges = genbininfo.get("ptgen", [])
        pois = [
            "r_ptgen{}_msdgen{}".format(ipt, imsd)
            for ipt in range(len(ptgen_edges)-1)
            for imsd in range(len(msdgen_edges)-1)
        ]
    print(pois)
    if pois == []:
        raise AttributeError("config does not have genbin info")
    try:
        fit_s = file_.Get("fit_s")
        poi_ind = [fit_s.floatParsFinal().index(p) for p in pois]
        corr_matrix = np.array([[fit_s.correlationMatrix()[row][col] for col in poi_ind] for row in poi_ind])
        cov_matrix = np.array([[fit_s.covarianceMatrix()[row][col] for col in poi_ind] for row in poi_ind])
        np.save("{}/poi_correlation_matrix.npy".format(model_dir),
                {
                    "pois": pois,
                    "correlationMatrix": corr_matrix,
                    "covarianceMatrix": cov_matrix,
                })
        del fit_s, corr_matrix
    except BaseException as e:
        print(e)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--corrmatrix", action="store_true")
    parser.add_argument("-i", "--input", default="")

    args = parser.parse_args()

    if args.extract:
        extract_fit_results(args.input)
    if args.corrmatrix:
        extract_correlation_matrix(args.input)
