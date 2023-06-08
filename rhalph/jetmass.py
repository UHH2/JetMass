#!/bin/env python
from __future__ import print_function
import sys
import os
import numpy as np
import common_configs
import ROOT  # type: ignore
ROOT.PyConfig.IgnoreCommandLineOptions = True
from jetmass_scale_fit_utils import scale_lumi, build_pseudo, build_mass_scale_variations, extract_fit_results  # noqa
import rhalphalib as rl         # noqa
rl.util.install_roofit_helpers()


def nuisance_name_(name, configs):
    result = (
        "%s_%s" % (name, configs["year"])
        if any(nuisance_name in name for nuisance_name in configs["nuisance_year_decorrelation"])
        else name
    )
    return result


def jet_mass_producer(args, configs):
    """
    configs: configuration dict including:
    ModelName,gridHistFileName,channels,histLocation
      -> channels: dict with dict for each channels:
        -> includes histDir,samples,NormUnc,signal,regions,QcdEstimation
    """
    rebin_msd = True
    binnings = {"W": np.linspace(50, 300, 26), "top": np.linspace(50, 300, 26)}
    binning_info = configs.get("binning", binnings["W"])
    min_msd, max_msd = (binning_info[0], binning_info[1])
    binwidth = binning_info[2]
    nbins = int(np.floor((max_msd - min_msd) / binwidth))
    msd_bins = np.linspace(min_msd, nbins * binwidth + min_msd, nbins + 1)

    year_str = configs["year"]

    def nuisance_name(name):
        return nuisance_name_(name, configs)

    xsec_priors = configs.get("xsec_priors", {})

    # channels for combined fit
    channels = configs["channels"]
    qcd_estimation_channels = {
        k: v for k, v in channels.items() if "QcdEstimation" in v and v["QcdEstimation"] == "True"
    }

    if args.verbose > 0:
        print("channels:", channels.keys())

    # getting path of dir with root file from config
    if args.verbose > 0:
        print("histfilename", configs["histLocation"])

    hist_file = ROOT.TFile(configs["histLocation"])
    aux_hist_files = {
        # "barrel": ROOT.TFile(configs["histLocation"].replace(".root", "_barrel.root")),
        # "endcap": ROOT.TFile(configs["histLocation"].replace(".root", "_endcap.root")),
    }

    for var in [
            "jec_up", "jec_down",
            "jer_up", "jer_down",
            "isr_up", "isr_down",
            "fsr_up", "fsr_down",
            "toppt_off",
    ]:
        fname = configs["histLocation"].replace(".root", "_{}.root".format(var))
        if os.path.isfile(fname):
            aux_hist_files[var] = ROOT.TFile(fname, "READ")

    do_qcd_estimation = len(qcd_estimation_channels) > 0
    do_initial_qcd_fit = configs.get("InitialQCDFit", "False") == "True"
    qcd_fail_region_constant = configs.get("QCDFailConstant", "False") == "True"
    qcd_fail_sigma_scale = configs.get("QCDSigmaScale", 1.0)
    TF_ranges = (-50, 50)
    qcdparam_lo, qcdparam_hi = (-50, 50)
    QCDFailUnbound = False

    lumi_scale = 1.0
    if "Pseudo" in configs and len(configs["Pseudo"]) > 0 and "lumiScale" in configs["Pseudo"][0]:
        lumi_scale = float(configs["Pseudo"][0].split(":")[-1])

    def get_hist(hist_dir, hist_file_="nominal"):
        # print(hist_dir)
        if hist_file_ == "nominal":
            hist = hist_file.Get(str(hist_dir))
        else:
            hist = aux_hist_files[hist_file_].Get(str(hist_dir))

        for sample_name_pattern, xsec_prior in xsec_priors.items():
            if sample_name_pattern in hist_dir:
                hist.Scale(xsec_prior)

        hist.SetName("msd")
        hist.GetXaxis().SetTitle("msd")
        if rebin_msd > 0:
            hist = hist.Rebin(len(msd_bins) - 1, "msd", msd_bins)
        if lumi_scale != 1.0:
            hist = scale_lumi(hist, lumi_scale)
        return hist

    model_name = configs.get("ModelName", "Jet_Mass_Model")  # get name from config, or fall back to default
    model_dir = "{}/{}".format(args.workdir, model_name)

    # specify if QCD estimation (using Bernstein-polynomial as TF) should be used
    ################
    # QCD Estimation#
    ################
    # derive pt bins from channel names for the pt,rho grid for the Bernstein-Polynomial

    def get_process_normalizations(process):
        norms = {"channels": [], "fail": [], "pass": []}
        for channel_name, config in qcd_estimation_channels.items():
            if args.verbose > 0:
                print(config["selection"] + "_%s__mjet_" % process + config["pt_bin"] + additional_bin + "_fail")
            fail_hist = get_hist(
                config["selection"] + "_%s__mjet_" % process + config["pt_bin"] + additional_bin + "_fail"
            )
            pass_hist = get_hist(
                config["selection"] + "_%s__mjet_" % process + config["pt_bin"] + additional_bin + "_pass"
            )
            # fail_hist = hist_file.Get('W_%s__mjet_'%process+config['pt_bin']+ additional_bin+'_fail')
            # pass_hist = hist_file.Get('W_%s__mjet_'%process+config['pt_bin']+ additional_bin+'_pass')
            # if(rebin_msd > 0):
            #     fail_hist = fail_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            #     pass_hist = pass_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            norms["channels"].append(channel_name)
            norms["fail"].append(fail_hist.Integral())
            norms["pass"].append(pass_hist.Integral())
        norms["fail"] = np.array(norms["fail"])
        norms["pass"] = np.array(norms["pass"])
        norms["eff"] = norms["pass"] / norms["fail"]
        channels_sorted = norms["channels"]
        channels_sorted.sort()
        # norms["eff_arr"] = np.array([[norms["eff"][i]] * (len(msd_bins) - 1) for i in range(len(norms["eff"]))])
        norms["eff_arr"] = np.array(
            [[norms["eff"][norms["channels"].index(ch)]] * (len(msd_bins) - 1) for ch in channels_sorted]
        )
        return norms

    if do_qcd_estimation:
        if args.verbose > 0:
            print("Doing some preparations for data driven QCD Estimate (Bernstein TF)")
        bernstein_orders = tuple(configs.get("BernsteinOrders", [2, 2]))
        qcd_model = rl.Model("qcdmodel")
        qcd_pass, qcd_fail = 0.0, 0.0
        qcd_estimation_relevant_selection = "W"
        for channel_name, config in qcd_estimation_channels.items():
            qcd_estimation_relevant_selection = config["selection"]
            fail_ch = rl.Channel(str(channel_name + "fail"))
            pass_ch = rl.Channel(str(channel_name + "pass"))
            fail_ch._renderExtArgs = not args.skipExtArgRender
            pass_ch._renderExtArgs = not args.skipExtArgRender
            qcd_model.addChannel(fail_ch)
            qcd_model.addChannel(pass_ch)
            additional_bin = config.get("additional_bin", "")

            fail_hist = get_hist(
                qcd_estimation_relevant_selection + "_QCD__mjet_" + config["pt_bin"] + additional_bin + "_fail"
            )
            pass_hist = get_hist(
                qcd_estimation_relevant_selection + "_QCD__mjet_" + config["pt_bin"] + additional_bin + "_pass"
            )

            # empty_hist = fail_hist.Clone()
            # empty_hist.Reset()
            # adding signal process just for the sake of being complete. won't be used in fit of qcd-model
            signal_fail_hist = get_hist(
                qcd_estimation_relevant_selection
                + "_"
                + config["signal"][0]
                + "__mjet_"
                + config["pt_bin"]
                + additional_bin
                + "_fail"
            )

            signal_pass_hist = get_hist(
                qcd_estimation_relevant_selection
                + "_"
                + config["signal"][0]
                + "__mjet_"
                + config["pt_bin"]
                + additional_bin
                + "_pass"
            )

            signal_fail = rl.TemplateSample(channel_name + "fail" + "_" + "Signal", rl.Sample.SIGNAL, signal_fail_hist)
            fail_ch.addSample(signal_fail)
            signal_pass = rl.TemplateSample(channel_name + "pass" + "_" + "Signal", rl.Sample.SIGNAL, signal_pass_hist)
            pass_ch.addSample(signal_pass)

            # creating pseudo data from qcd mc
            h_np_fail = rl.util._to_numpy(fail_hist)
            h_np_pseudo_fail = (np.random.poisson(h_np_fail[0]),) + h_np_fail[1:]
            h_np_pass = rl.util._to_numpy(pass_hist)
            h_np_pseudo_pass = (np.random.poisson(h_np_pass[0]), ) + h_np_pass[1:]

            fail_ch.setObservation(h_np_pseudo_fail)
            pass_ch.setObservation(h_np_pseudo_pass)

            # fail_
            qcd_fail += fail_ch.getObservation().sum()
            qcd_pass += pass_ch.getObservation().sum()
        qcd_eff = qcd_pass / qcd_fail
        qcd_norms = get_process_normalizations("QCD")
        data_norms = get_process_normalizations("Data")

        # get all lower edges from channel names
        # pt_edges = configs.get('pt_edges',[500,550,600,675,800,1200])
        pt_edges = configs.get("pt_edges", [500, 650, 800, 1200])
        pt_bins = np.array([(1400 if pt == "Inf" else pt) for pt in pt_edges])
        # pt_bins = np.array([500, 550, 600, 675, 800, 1200])
        msd = rl.Observable("msd", msd_bins)

        # here we derive these all at once with 2D array
        ptpts, msdpts = np.meshgrid(
            pt_bins[:-1] + 0.3 * np.diff(pt_bins), msd_bins[:-1] + 0.5 * np.diff(msd_bins), indexing="ij"
        )
        rhopts = 2 * np.log(msdpts / ptpts)
        ptscaled = (ptpts - pt_bins[0]) / (pt_bins[-1] - pt_bins[0])

        # remove lower bound on rho, since we removed it also from selection
        # for now we do this by simply setting rho_min to bounds specified by pt and mass bound
        rho_min = 2*np.log(msd_bins[0]/pt_bins[-1])
        rho_max = -2.1
        rhoscaled = (rhopts - (rho_min)) / ((rho_max) - (rho_min))
        validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
        rhoscaled[~validbins] = 1  # we will mask these out later

        TF_suffix = configs.get("TFSuffix", "")

        tf1_basis = "Bernstein"
        tf2_basis = "Bernstein"

        if do_initial_qcd_fit:
            initial_qcd_fit_orders = tuple(configs.get("InitialQCDFitOrders", [2, 2]))
            if not os.path.exists(model_dir):
                os.makedirs(model_dir)
            if args.verbose > 0:
                print("QCD eff:", qcd_eff)
            tf_MCtempl = rl.BasisPoly(
                "tf_MCtempl_%s" % year_str,
                initial_qcd_fit_orders,
                ["pt", "rho"],
                basis=tf1_basis,
                init_params=np.ones((initial_qcd_fit_orders[0] + 1, initial_qcd_fit_orders[1] + 1)),
                limits=TF_ranges,
            )
            tf_MCtempl_params = qcd_norms["eff_arr"] * tf_MCtempl(ptscaled, rhoscaled)
            for channel_name, config in qcd_estimation_channels.items():
                # ptbin = np.where(pt_bins==float(channel_name.split('Pt')[-1]))[0][0]
                ptbin = np.where(pt_bins == float(config["pt_bin"].split("to")[0]))[0][0]
                failCh = qcd_model[channel_name + "fail"]
                passCh = qcd_model[channel_name + "pass"]
                failObs = failCh.getObservation()
                if qcd_fail_region_constant and args.verbose > 0:
                    print("Setting QCD parameters in fail region constant")
                qcdparams = np.array(
                    [
                        rl.IndependentParameter(
                            "qcdparam_ptbin%d_msdbin%d_%s" % (ptbin, i, year_str),
                            0,
                            constant=qcd_fail_region_constant,
                            lo=qcdparam_lo,
                            hi=qcdparam_hi,
                        )
                        for i in range(msd.nbins)
                    ]
                )
                for param in qcdparams:
                    param.unbound = QCDFailUnbound

                # scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
                scaledparams = failObs * (1 + qcd_fail_sigma_scale / 100.0) ** qcdparams
                fail_qcd = rl.ParametericSample(
                    "%sfail_qcd_%s" % (channel_name, year_str), rl.Sample.BACKGROUND, msd, scaledparams
                )
                failCh.addSample(fail_qcd)
                pass_qcd = rl.TransferFactorSample(
                    "%spass_qcd_%s" % (channel_name, year_str),
                    rl.Sample.BACKGROUND,
                    tf_MCtempl_params[ptbin, :],
                    fail_qcd,
                )
                passCh.addSample(pass_qcd)

                failCh.mask = validbins[ptbin]
                passCh.mask = validbins[ptbin]

            qcd_model.renderCombine(model_dir + "/qcdmodel")

            qcdfit_ws = ROOT.RooWorkspace("w")
            simpdf, obs = qcd_model.renderRoofit(qcdfit_ws)
            ROOT.Math.MinimizerOptions.SetDefaultPrecision(1e-18)
            ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
            ROOT.Math.MinimizerOptions.SetDefaultTolerance(0.0001)
            # ROOT.Math.MinimizerOptions.SetDefaultPrecision(-1.0)

            qcdfit = simpdf.fitTo(
                obs,
                ROOT.RooFit.Extended(True),
                ROOT.RooFit.SumW2Error(False),
                ROOT.RooFit.Strategy(2),
                ROOT.RooFit.Save(),
                ROOT.RooFit.Minimizer("Minuit2", "migrad"),
                ROOT.RooFit.Offset(True),
                ROOT.RooFit.PrintLevel(1),
            )
            # qcdfit = simpdf.fitTo(
            #     obs,
            #     ROOT.RooFit.Extended(True),
            #     ROOT.RooFit.SumW2Error(True),
            #     ROOT.RooFit.Strategy(1),
            #     ROOT.RooFit.Save(),
            #     ROOT.RooFit.Minimizer("Minuit2", "migrad"),
            #     # ROOT.RooFit.PrintLevel(-1),
            #     ROOT.RooFit.PrintLevel(1),
            #     ROOT.RooFit.Minos(0),
            # )

            qcdfit_ws.add(qcdfit)
            qcdfit_ws.writeToFile(model_dir + ("/qcdfit_%s_" % year_str) + model_name + TF_suffix + ".root")
            if qcdfit.status() != 0:
                raise RuntimeError("Could not fit qcd")

            qcd_model.readRooFitResult(qcdfit)

            param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
            decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + "_deco", qcdfit, param_names)
            tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
            tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

            tf_dataResidual = rl.BasisPoly(
                "tf_dataResidual_%s" % year_str, bernstein_orders, ["pt", "rho"],
                basis=tf2_basis, limits=TF_ranges
            )
            tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
            tf_params = data_norms["eff_arr"] * tf_MCtempl_params_final * tf_dataResidual_params
        else:
            tf_params = None  # define later

    # Reading categories of consituent-variations for nuisance paramters from gridHist

    grid_nuisances, _ = build_mass_scale_variations(configs, args)

    # setting up rhalphalib roofit model
    model = rl.Model(model_name)

    if args.unfolding:
        for channel_name, config in channels.items():
            # Signal sample division in genbins
            from copy import deepcopy

            signal_samples = deepcopy(config["signal"])
            unfolding_bins = configs["unfolding_bins"]
            genbins = [
                    "ptgen%i_msdgen%i" % (iptgen, imsdgen)
                    for iptgen in range(len(unfolding_bins["ptgen"]) - 1)
                    for imsdgen in range(len(unfolding_bins["msdgen"]) - 1)
                ]
            configs["genbins"] = genbins
            for signal_sample in signal_samples:
                config["samples"].remove(signal_sample)
                config["signal"].remove(signal_sample)
                for genbin in genbins:
                    sample_genbin_name = "{SAMPLE}_{GENBIN}".format(SAMPLE=signal_sample, GENBIN=genbin)
                    config["samples"].append(sample_genbin_name)
                    config["signal"].append(sample_genbin_name)

    # setting up nuisances for systematic uncertainties
    lumi = rl.NuisanceParameter(nuisance_name("CMS_lumi"), "lnN")
    lumi_effect = 1.027

    pt_bins_for_jmr_nuisances = list(configs["channels"].keys())
    pt_bins_for_jmr_nuisances.append("inclusive")
    jmr_nuisances = {
        parname: rl.NuisanceParameter(parname, "shape", 0, -10, 10)
        for parname in [
            nuisance_name("jmr_variation_PT%s" % ch_name.split("Pt")[-1]) for ch_name in pt_bins_for_jmr_nuisances
        ]
    }

    norm_nuisances = {}
    initial_qcd_from_data = {}
    for channel_name in channels.keys():
        if args.minimalModel:
            break
        if "NormUnc" not in channels[channel_name]:
            continue
        norm_uncertainty_pars = {}
        norm_uncertainties = channels[channel_name]["NormUnc"]
        for name, norm_unc in norm_uncertainties.items():
            norm_uncertainty_regions = [""]
            norm_unc_val = 0.0
            if isinstance(norm_unc, dict):
                if norm_unc.get("decorrelateRegions", False):
                    print("decorrelating normUnc {} accross regions.".format(name))
                    norm_uncertainty_regions = channels[channel_name]["regions"]
                    norm_unc_val = norm_unc["value"]
                else:
                    raise NotImplementedError(
                        'If NormUnc is a dict it should be of the form {"value":XX.XX,"decorrelateRegions":True}.'
                        "Anything else is not implemented at the moment."
                    )
            elif isinstance(norm_unc, float):
                norm_unc_val = norm_unc
            else:
                raise NotImplementedError("provided form of NormUnc parameter is not implemented!")
            norm_uncertainty_pars.update(
                {
                    name: {
                        region: [
                            rl.NuisanceParameter(nuisance_name("%s%s_normUnc" % (name, region)), "lnN"),
                            norm_unc_val,
                        ]
                        for region in norm_uncertainty_regions
                    }
                }
            )

        for sample in channels[channel_name]["samples"]:
            if sample in norm_nuisances:
                continue
            for nuisance_par_name, nuisance_par_dict in norm_uncertainty_pars.items():

                if nuisance_par_name in sample:
                    norm_nuisances[sample] = nuisance_par_dict

    # jec variation nuisances
    jec_var_nuisance = rl.NuisanceParameter(nuisance_name("jec_variation"), "shape", 0, -10, 10)
    extra_nuisances = {
        "isr": rl.NuisanceParameter(nuisance_name("isr_variation"), "shape", 0, -10, 10),
        "fsr": rl.NuisanceParameter(nuisance_name("fsr_variation"), "shape", 0, -10, 10),
        "toppt": rl.NuisanceParameter(nuisance_name("toppt_reweight"), "shape", 0, -10, 10),
    }
    # tagging eff sf for ttbar semileptonic samples

    # top_tag_eff = rl.IndependentParameter("top_tag_eff_sf",1.,-10,10)
    # W_tag_eff = rl.IndependentParameter("W_tag_eff_sf",1.,-10,10)
    for channel_name, config in channels.items():
        top_tag_eff = rl.IndependentParameter(nuisance_name("%stop_tag_eff_sf" % channel_name), 1.0, -4, 4)
        W_tag_eff = rl.IndependentParameter(nuisance_name("%sW_tag_eff_sf" % channel_name), 1.0, -4, 4)

        # using hists with /variable/ in their name (default: Mass, if defined get from config)
        variable = "mjet" if "variable" not in config else config["variable"]
        # getting list of samples from config
        if args.minimalModel:
            config["samples"] = ["QCD", "WJetsMatched"]
        samples = config["samples"]
        # for WMass fit there are multiple regions per sample
        regions = [""] if "regions" not in config else config["regions"]
        if args.verbose > 0:
            print("setting up channel:", channel_name)
            print("getting template of variable:", variable)
            print("samples:", samples)
            print("regions:", regions)

        for region in regions:
            additional_bin = config.get("additional_bin", "")
            region_suffix = "_" + region if len(region) > 0 else ""
            hist_dir = (
                config["selection"] + "_%s__" + variable + "_%s" + config["pt_bin"] + additional_bin + region_suffix
            )
            # setting up channel for fit (name must be unique and can't include any '_')
            region_name = str(channel_name + region)
            ch = rl.Channel(region_name)
            ch._renderExtArgs = not args.skipExtArgRender
            model.addChannel(ch)
            if args.verbose > 0:
                print("hist_dir:", hist_dir)
                print("rl.Channel:", ch)
            # if not args.noNormUnc:
            #     ch.addParamGroup(
            #         "NormUnc", set([norm_nuisances[sample_name][0] for sample_name in norm_nuisances.keys()])
            #     )

            for sample_name in samples:
                # do not include QCD template here, but rather use qcd estimation below
                if ("QcdEstimation" in config and config["QcdEstimation"] == "True") and "qcd" in sample_name.lower():
                    continue

                # specify if sample is signal or background type
                sample_type = rl.Sample.SIGNAL if sample_name in config["signal"] else rl.Sample.BACKGROUND
                sample_hist = get_hist(hist_dir % (sample_name, ""))

                # setup actual rhalphalib sample
                sample = rl.TemplateSample(ch.name + "_" + sample_name, sample_type, sample_hist)

                # if we want to use massScales in fit, we add the actual ParamEffect to NuisanceParamters,
                # this will get them rendered into the workspace.
                if args.massScales or True:
                    # setting effects of constituent variation nuisances (up/down)
                    for grid_nuisance_dict, x, y, category in grid_nuisances:
                        for grid_nuisance_name in grid_nuisance_dict.keys():
                            # this should strictly not be necessary
                            if sample.sampletype != rl.Sample.SIGNAL and args.VaryOnlySignal:
                                continue

                            # at this point this just decides to skip massScale based on pT bin since
                            # mass_scale_substring hold just info about pT bin
                            mass_scale_substring = (
                                "PT" + channel_name.split("Pt")[-1] if args.pTdependetMassScale else "inclusive"
                            )
                            if mass_scale_substring not in grid_nuisance_name:
                                continue

                            # getting rhalphalib paramter from setup_dict
                            grid_nuisance = grid_nuisance_dict[grid_nuisance_name]["nuisance"]

                            # getting and "applying" info on which region and sample
                            # the respective massScale shoud act on
                            application_info = grid_nuisance_dict[grid_nuisance_name]["application_info"]
                            if sample_name not in application_info["samples"] and application_info["samples"] not in [
                                    "all",
                                    "signal",
                            ]:
                                continue
                            if region not in application_info["regions"]:
                                continue

                            # getting hists with up/down varaition for Parameter Effect
                            hist_up = get_hist(
                                hist_dir % (sample_name, str(x) + "_" + str(y) + "_" + category + "_") + "__up"
                            )
                            hist_down = get_hist(
                                hist_dir % (sample_name, str(x) + "_" + str(y) + "_" + category + "_") + "__down"
                            )

                            sample.setParamEffect(
                                grid_nuisance, hist_up, hist_down, scale=configs.get("massScaleFactor", 1.0)
                            )

                if not args.noNuisances:
                    if args.JECVar and (
                        sample.sampletype == rl.Sample.SIGNAL or sample.sampletype == rl.Sample.BACKGROUND
                    ):
                        hist_jec_up = get_hist(hist_dir % (sample_name, ""), "jec_up")
                        hist_jec_down = get_hist(hist_dir % (sample_name, ""), "jec_down")
                        sample.setParamEffect(jec_var_nuisance, hist_jec_up, hist_jec_down)

                    # ISR down
                    if sample.sampletype == rl.Sample.SIGNAL and (
                        "isr_up" in aux_hist_files and "isr_down" in aux_hist_files
                    ):
                        hist_isr_up = get_hist(hist_dir % (sample_name, ""), "isr_up")
                        hist_isr_down = get_hist(hist_dir % (sample_name, ""), "isr_down")
                        sample.setParamEffect(extra_nuisances["isr"], hist_isr_up, hist_isr_down)
                    # FSR down
                    if sample.sampletype == rl.Sample.SIGNAL and (
                        "fsr_up" in aux_hist_files and "fsr_down" in aux_hist_files
                    ):
                        hist_fsr_up = get_hist(hist_dir % (sample_name, ""), "fsr_up")
                        hist_fsr_down = get_hist(hist_dir % (sample_name, ""), "fsr_down")
                        sample.setParamEffect(extra_nuisances["fsr"], hist_fsr_up, hist_fsr_down)

                    if "TTTo" in sample.name and "toppt_off" in aux_hist_files:
                        hist_toppt_off = get_hist(hist_dir % (sample_name, ""), "toppt_off")
                        sample.setParamEffect(extra_nuisances["toppt"], effect_up=hist_toppt_off)  # , scale=0.5)

                    # setting effects of JMR variation nuisance(s)
                    if args.JMRparameter and sample.sampletype == rl.Sample.SIGNAL:
                        hist_jmr_up = get_hist(hist_dir % (sample_name, ""), "jer_up")
                        hist_jmr_down = get_hist(hist_dir % (sample_name, ""), "jer_down")
                        if args.pTdependetJMRParameter:
                            sample.setParamEffect(
                                jmr_nuisances.get(
                                    nuisance_name("jmr_variation_PT%s_%s" % channel_name.split("Pt")[-1])
                                ),
                                hist_jmr_up,
                                hist_jmr_down,
                            )
                        else:
                            sample.setParamEffect(
                                jmr_nuisances.get(nuisance_name("jmr_variation_PTinclusive_%s")),
                                hist_jmr_up,
                                hist_jmr_down,
                            )
                    # other nuisances (lumi, norm unc)
                    sample.setParamEffect(lumi, lumi_effect)
                    if not args.noNormUnc:
                        if sample_name in norm_nuisances.keys():
                            for norm_nuisance_region in norm_nuisances[sample_name].keys():
                                if norm_nuisance_region == region or norm_nuisance_region == "":
                                    sample.setParamEffect(
                                        norm_nuisances[sample_name][norm_nuisance_region][0],
                                        norm_nuisances[sample_name][norm_nuisance_region][1],
                                    )

                ch.addSample(sample)

            PseudoData = "Pseudo" in configs and len(configs["Pseudo"]) > 0
            if PseudoData:
                # try:
                data_hist = build_pseudo(
                    samples, hist_file, hist_dir, configs["Pseudo"], args.minimalModel, msd_bins=msd_bins
                )
                # except:
                # data_hist = build_pseudo_uproot(samples,hist_dir,msd_bins,configs)
            else:
                if args.verbose > 0:
                    print("using data!!!!!")
                data_hist = get_hist(hist_dir % ("Data", ""))

            # if(rebin_msd > 0):
            #     data_hist = data_hist.Rebin(len(msd_bins)-1, 'msd', msd_bins)
            # data_hist.SetName('msd')
            ch.setObservation(data_hist, read_sumw2=PseudoData)
            if "QcdEstimation" in config and config["QcdEstimation"] == "True":
                mask = validbins[np.where(pt_bins == float(config["pt_bin"].split("to")[0]))[0][0]]
                # dropped_events = np.sum(ch.getObservation().astype(float)[~mask])
                # percentage = dropped_events/np.sum(ch.getObservation().astype(float))
                # print(
                #     "dropping due to mask: %.2f events (out of %.2f -> %.2f%%)"
                #     % (dropped_events, np.sum(ch.getObservation().astype(float)), percentage * 100)
                # )
                ch.mask = mask

        # setting effect for tagging eff scale factors
        # top tagging
        if args.TTbarTaggingEff and config["selection"] == "top":
            top_pass_sample = model[channel_name + "pass"]["TTToSemiLeptonic_mergedTop"]
            top_passW_sample = model[channel_name + "passW"]["TTToSemiLeptonic_mergedTop"]
            top_fail_sample = model[channel_name + "fail"]["TTToSemiLeptonic_mergedTop"]
            rpf_top_Wfail = top_pass_sample.getExpectation(nominal=True).sum() / (
                top_passW_sample.getExpectation(nominal=True).sum() + top_fail_sample.getExpectation(nominal=True).sum()
            )
            top_pass_sample.setParamEffect(top_tag_eff, 1.0 * top_tag_eff)
            top_passW_sample.setParamEffect(top_tag_eff, (1 - top_tag_eff) * rpf_top_Wfail + 1.0)
            top_fail_sample.setParamEffect(top_tag_eff, (1 - top_tag_eff) * rpf_top_Wfail + 1.0)
            # W tagging
            W_pass_sample = model[channel_name + "pass"]["TTToSemiLeptonic_mergedW"]
            W_passW_sample = model[channel_name + "passW"]["TTToSemiLeptonic_mergedW"]
            W_fail_sample = model[channel_name + "fail"]["TTToSemiLeptonic_mergedW"]
            rpf_W_topfail = W_passW_sample.getExpectation(nominal=True).sum() / (
                W_pass_sample.getExpectation(nominal=True).sum() + W_fail_sample.getExpectation(nominal=True).sum()
            )
            W_passW_sample.setParamEffect(W_tag_eff, 1.0 * W_tag_eff)
            W_pass_sample.setParamEffect(W_tag_eff, (1 - W_tag_eff) * rpf_W_topfail + 1.0)
            W_fail_sample.setParamEffect(W_tag_eff, (1 - W_tag_eff) * rpf_W_topfail + 1.0)

    if do_qcd_estimation:
        # QCD TF
        if not do_initial_qcd_fit:
            tf_params = rl.BasisPoly(
                "tf_params_%s" % year_str, bernstein_orders, ["pt", "rho"],
                basis=tf1_basis, limits=TF_ranges
            )

            # for a in tf_params.parameters:
            #     for params in a:
            #         params.unbound = True
            if args.verbose > 0:
                print("Using QCD efficiency (N2-ddt) of %.2f%% to scale initial QCD in pass region" % (qcd_eff * 100))
            tf_params = data_norms["eff_arr"] * tf_params(ptscaled, rhoscaled)

        for channel_name, config in channels.items():
            if "QcdEstimation" not in config or config["QcdEstimation"] == "False":
                continue
            fail_ch = model[channel_name + "fail"]
            pass_ch = model[channel_name + "pass"]
            ptbin = np.where(pt_bins == float(config["pt_bin"].split("to")[0]))[0][0]
            if qcd_fail_region_constant and args.verbose > 0:
                print("Setting QCD parameters in fail region constant")
            qcd_params = np.array(
                [
                    rl.IndependentParameter(
                        "qcdparam_ptbin%i_msdbin%i_%s" % (ptbin, i, year_str),
                        0,
                        constant=qcd_fail_region_constant,
                        lo=qcdparam_lo,
                        hi=qcdparam_hi,
                    )
                    for i in range(msd.nbins)
                ]
            )

            for param in qcd_params:
                param.unbound = QCDFailUnbound

            initial_qcd = (
                fail_ch.getObservation()[0].astype(float)
                if isinstance(fail_ch.getObservation(), tuple)
                else fail_ch.getObservation().astype(float)
            )
            for sample in fail_ch:
                initial_qcd -= sample.getExpectation(nominal=True)
            if np.any(initial_qcd < 0.0):
                initial_qcd = np.where(initial_qcd <= 0.0, 0, initial_qcd)
                if args.verbose > 0:
                    print("negative bins in initial_qcd in ", channel_name)
                # continue
                minimum = np.amin(initial_qcd)
                initial_qcd = np.where(initial_qcd == 0, minimum, initial_qcd)
                initial_qcd += abs(minimum)
                raise ValueError("inital qcd (fail qcd from data - mc) negative at least one bin")

            # scaledparams = initial_qcd * (1 + sigmascale / np.maximum(1.0, np.sqrt(initial_qcd))) ** qcd_params
            # scaledparams = (
            #     initial_qcd * (1 + qcd_fail_sigma_scale / np.maximum(1.0, np.sqrt(initial_qcd))) ** qcd_params
            # )

            initial_qcd_from_data[channel_name] = initial_qcd
            scaledparams = initial_qcd * (1 + qcd_fail_sigma_scale / 100.0) ** qcd_params
            fail_qcd = rl.ParametericSample(
                "%sfail_qcd_%s" % (channel_name, year_str), rl.Sample.BACKGROUND, msd, scaledparams
            )
            fail_ch.addSample(fail_qcd)
            # fail_ch.addParamGroup("QCDFail",fail_qcd.parameters)
            pass_qcd = rl.TransferFactorSample(
                "%spass_qcd_%s" % (channel_name, year_str), rl.Sample.BACKGROUND, tf_params[ptbin, :], fail_qcd
            )
            pass_ch.addSample(pass_qcd)
            raw_tf_params = []
            for param in pass_qcd.parameters:
                if "tf" in param.name:
                    raw_tf_params.append(param)
            # pass_ch.addParamGroup("QCDPass",raw_tf_params)

    if args.prefitAsimov:
        for c in model:
            prefit_asimov = initial_qcd_from_data[c.name.replace("pass", "").replace("fail", "")].copy()
            if "pass" in c.name:
                prefit_asimov *= data_norms["eff"][list(configs["channels"].keys()).index(c.name.replace("pass", ""))]
            for sample in c.samples:
                if "qcd" in sample.name.lower():
                    continue
                prefit_asimov += sample.getExpectation(nominal=True)
            prefit_asimov_data = (prefit_asimov, c.observable.binning, c.observable.name)
            c.setObservation(prefit_asimov_data)

    model.renderCombine(model_dir)


if __name__ == "__main__":
    import json
    import argparse
    import fitplotter
    from CombineWorkflows import CombineWorkflows

    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with config")
    parser.add_argument("--workdir", type=str, help="path to workdir containing modeldir", default=".")
    parser.add_argument("--justplots", action="store_true", help="just make plots.")
    parser.add_argument("--build", action="store_true", help="just build workspace for combine")
    parser.add_argument("--job_index", default="", type=str)
    parser.add_argument("--minimalModel", action="store_true")
    parser.add_argument("--forceCommonConfig", action="store_true")
    parser.add_argument("--skipTemplatePlots", action="store_true")
    parser.add_argument("--customCombineWrapper", action="store_true")
    parser.add_argument("--noNuisances", action="store_true")
    parser.add_argument("--JMRparameter", action="store_true")
    parser.add_argument("--prefitAsimov", action="store_true")
    parser.add_argument("--pTdependetJMRParameter", action="store_true")
    parser.add_argument("--noNormUnc", action="store_true")
    parser.add_argument("--skipExtArgRender", action="store_true")
    parser.add_argument("--seed", type=str, default="42")
    parser.add_argument("--freezeParameters", nargs="+", default=[])
    parser.add_argument(
        "--initialQCDTF",
        action="store_true",
        help="setup dual TF with inital TF from fit to QCD MC and residual Data TF",
    )
    parser.add_argument(
        "--combineOptions", type=str, help="string with additional cli-options passed to combine", default=""
    )
    parser.add_argument("--verbose", type=int, default=-1)

    parser.add_argument("-M", "--mode", choices=["unfolding", "jms", "default"],
                        help="choose what fit should be run. 'unfolding' performs MaxLik. unfolding, jms measures "
                        "JetMassScale variations,'default' measures signal cross section.",
                        required=True)
    parser.add_argument(
        "--uncertainty-breakdown",
        default=[],
        nargs="+",
        help="specify nuisances for which uncertainty breakdown plots should be created",
    )

    parser.add_argument(
        "--tagger",
        type=str,
        help=(
            "tagger of main selection. The string will be insterted into histLocation in the configs dict,"
            "until i come up with a better solution."  # TODO: come up with a better solution
        ),
        default=""
    )

    args = parser.parse_args()

    args.unfolding = args.mode == "unfolding"
    args.massScales = args.mode == "jms"
    args.defaultPOI = args.mode == "default"

    try:
        if args.config.endswith(".json"):
            configs = json.load(open(args.config))
        else:
            execfile(args.config)  # noqa # type: ignore
        existing_config = args.workdir + "/" + configs["ModelName"] + "/config.json"
        if os.path.isfile(existing_config) and args.unfolding:
            use_existing_config = (
                raw_input(      # noqa # type: ignore
                    "There already is a directory corresponding to this config. "
                    "Do you want to load the existing config? [Y/N]"
                ).lower()
                == "y"
            )
            if use_existing_config:
                configs = json.load(open(existing_config))
    except IndexError:
        print("You must specify a configuration JSON!")
        sys.exit(0)

    configs["ModelName"] = configs["ModelName"] + str(args.job_index)

    model_dir = "{}/{}".format(args.workdir, configs["ModelName"])
    configs["ModelDir"] = model_dir

    particlenet = "particlenet" in args.tagger
    if "BernsteinOrders" not in configs:
        print("taking default BernsteinOrders from common configs..")
        configs["BernsteinOrders"] = common_configs.bernstein_orders(
            configs["year"], TF="Data", particlenet=particlenet
        )
        if args.initialQCDTF:
            configs["BernsteinOrders"] = common_configs.bernstein_orders(
                configs["year"], TF="Data2TF", particlenet=particlenet
            )

    if args.initialQCDTF:
        if "InitialQCDFitOrders" not in configs or args.forceCommonConfig:
            configs["InitialQCDFitOrders"] = common_configs.bernstein_orders(
                configs["year"], TF="QCD", particlenet=particlenet
            )
        configs["InitialQCDFit"] = "True"

    if not args.tagger.startswith("_") and args.tagger != "":
        args.tagger = "_" + args.tagger
    configs["histLocation"] = configs["histLocation"].replace(configs["year"], configs["year"] + args.tagger)

    configs["nuisance_year_decorrelation"] = [
        "CMS_lumi",
        "jec_variation", "isr_variation", "fsr_variation",
        "toppt_reweight", "tag_eff_sf", "jec_variation"
    ]

    args.freezeParameters = [nuisance_name_(par_name, configs) for par_name in args.freezeParameters]
    args.TTbarTaggingEff = configs.get("TTbarTaggingEff", "True") == "True"
    args.pTdependetMassScale = configs.get("pTdependentMassScale", "True") == "True"
    args.separateMassScales = configs.get("separateMassScales", "False") == "True"
    args.VaryOnlySignal = configs.get("VaryOnlySignal", "False") == "True"
    args.JECVar = configs.get("JECVar", "True") == "True"

    if not args.justplots:
        jet_mass_producer(args, configs)
        open(model_dir + "/config.json", "w").write(json.dumps(configs, sort_keys=False, indent=2))
        if not args.customCombineWrapper:
            cw = CombineWorkflows(build_only=args.build)
            if args.skipTemplatePlots:
                cw._skip_plotting = True
            cw.workspace = model_dir + "/model_combined.root"
            if args.unfolding:
                cw.method = "unfolding"
                cw.write_wrapper()
            else:
                cw.extraOptions = args.combineOptions
                if args.massScales:
                    cw.freezeParameters = ["r"]+args.freezeParameters
                    cw.extraOptions += " --preFitValue 0"
                    cw.POIRange = (-100, 100)
                if args.defaultPOI:
                    cw.POIRange = (0.01, 100.0)
                cw.POI = "r" if args.defaultPOI else build_mass_scale_variations(configs, args)[1]
                cw.method = "diagnostics"
                cw.write_wrapper()
                # if args.massScales:
                #     cw.method = "FastScanMassScales"
                #     cw.write_wrapper(append=True)
        else:
            if not os.path.isfile(model_dir + "/wrapper.sh"):
                import warnings
                warnings.warn(
                    (
                        "\033[93mYou used the option --CustomCombineWrapper,"
                        + " but no wrapper can be found in the modeldir!\033[0m"
                    ),
                    RuntimeWarning,
                )
            # write_wrapper(
            #     [configs["ModelName"]],
            #     "r" if use_r_poi else mass_scale_names,
            #     additional_options=args.combineOptions,
            #     combineWorkflow=args.combineWorkflow,
            # )
        if args.build:
            os.chdir(model_dir)
            os.system("bash build.sh")
            exit(0)
        # from runFit import runFits
        # runFits([configs['ModelName']])
        # exedir = os.getcwd()
        os.system("bash " + model_dir + "/wrapper.sh")
        # os.system("cd "+exedir)

    if args.customCombineWrapper:
        exit(0)

    do_postfit = extract_fit_results(configs)

    if not args.skipTemplatePlots:
        fitplotter.plot_fit_result(
            configs,
            plot_total_sig_bkg=False,
            do_postfit=do_postfit,
            use_config_samples=args.unfolding,
            pseudo_data=False,
        )
        fitplotter.plot_fit_result(
            configs,
            logY=True,
            plot_total_sig_bkg=False,
            do_postfit=do_postfit,
            use_config_samples=args.unfolding,
            pseudo_data=False,
        )
    if do_postfit and args.massScales:
        fitplotter.plot_mass_scale_nuisances(configs)

    qcd_estimation_channels = {
        k: v for k, v in configs["channels"].items() if "QcdEstimation" in v and v["QcdEstimation"] == "True"
    }
    if len(qcd_estimation_channels) > 0 and do_postfit:
        fitplotter.plot_qcd_bernstein(configs, do_3d_plot=False)
        if configs.get("QCDFailConstant", "False") == "False":
            fitplotter.plot_qcd_fail_parameters(configs)

    if args.JECVar and not args.unfolding and all("jec_variation" not in p for p in args.freezeParameters):
        args.uncertainty_breakdown.append(nuisance_name_("jec_variation", configs))
    if len(args.uncertainty_breakdown) > 0:
        os.system(
            "./scan_poi.py {} --poi all --plots --lasteffect rest -f {}".format(
                model_dir, " ".join(args.uncertainty_breakdown)
            )
        )
