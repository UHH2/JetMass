#!/usr/bin/env pythonJMS.sh
from GenericJobSubmitter import JobSubmitter
from JMSPlotter import create_plotter, finalize_ax, setup_ax
import argparse
import os
import json


def exec_bash(cmd, debug=False):
    print(cmd)
    if not debug:
        os.system(cmd)


def add_fits(taggers, samples, years, common_wrapper_lines):
    for tagger in taggers:
        # tagger_argument = " --tagger particlenet --initialQCDTF" if tagger == "ParticleNet" else ""
        tagger_argument = " --tagger particlenet " if tagger == "ParticleNet" else ""

        for sample in samples:

            # separate years
            for year in years:
                # JEC Nuisance
                workdir = "fits_years_{}_noJEC".format(tagger)
                wrapper_lines = common_wrapper_lines + [
                    "python jetmass.py -M jms configs/noJEC/{}_{}_noJEC.py --workdir {}/{}".format(
                        sample, year, base_dir, workdir
                    )
                    + tagger_argument,
                ]
                submitter.add_job("{}{}".format(sample, year), workdir, wrapper_lines)

                # JEC Nuisance
                workdir = "fits_years_{}_JECNuis".format(tagger)
                wrapper_lines = common_wrapper_lines + [
                    "python jetmass.py -M jms configs/nominal/{}_{}.py --workdir {}/{}".format(
                        sample, year, base_dir, workdir
                    )
                    + tagger_argument,
                ]
                submitter.add_job("{}{}".format(sample, year), workdir, wrapper_lines)

                # JEC Unc agnostic
                workdir = "fits_years_{}_JECVar".format(tagger)
                wrapper_lines = common_wrapper_lines + [
                    "python jetmass.py -M jms configs/nominal/{}_{}.py --workdir {}/{}".format(
                        sample, year, base_dir, workdir
                    )
                    + " --freezeParameters jec_variation"
                    + tagger_argument,
                ]
                submitter.add_job("{}{}".format(sample, year), workdir, wrapper_lines)

                # JEC Unc Var as input
                for direction in ["up", "down"]:
                    wrapper_lines = common_wrapper_lines + [
                        "python jetmass.py -M jms configs/JECVar/{}_{}_JEC_{}.py --workdir {}/{}".format(
                            sample, year, direction, base_dir, workdir
                        )
                        + tagger_argument,
                    ]
                    submitter.add_job("{}{}JEC{}".format(sample, year, direction.upper()), workdir, wrapper_lines)

                # separate MassScales
                workdir = "fits_years_{}_separateMassScale".format(tagger)
                wrapper_lines = common_wrapper_lines + [
                    (
                        "python jetmass.py -M jms configs/nominal/{}_{}.py --workdir {}/{}".format(
                            sample, year, base_dir, workdir
                        )
                        + " --separateMassScales "
                        + " --freezeParameters jec_variation"
                        + tagger_argument
                    )
                ]
                submitter.add_job("{}{}".format(sample, year), workdir, wrapper_lines)

            # # combined Run II
            # # JEC Nuisance
            # workdir = "fits_RunII_{}_JECNuis".format(tagger)
            # wrapper_lines = common_wrapper_lines + [
            #     (
            #         "python jetmass_combination.py -M jms --configs configs/nominal/{}*.py"
            #         " --workdir {}/{} --name {}RunII"
            #     ).format(sample, base_dir, workdir, sample)
            #     + ' --extra-options "{} --uncertainty_breakdown"'.format(tagger_argument)
            # ]
            # submitter.add_job("{}RunII".format(sample), workdir, wrapper_lines)

            # # JEC Unc agnostic
            # workdir = "fits_RunII_{}_JECVar".format(tagger)
            # wrapper_lines = common_wrapper_lines + [
            #     (
            #         "python jetmass_combination.py -M jms --configs configs/nominal/{}*.py"
            #         " --workdir {}/{} --name {}RunII"
            #     ).format(sample, base_dir, workdir, sample)
            #     + ' --extra-options "--freezeParameters jec_variation {} "'.format(tagger_argument),
            # ]
            # submitter.add_job("{}RunII".format(sample), workdir, wrapper_lines)

            # # JEC Unc Var as input
            # for direction in ["up", "down"]:
            #     wrapper_lines = common_wrapper_lines + [
            #         (
            #             "python jetmass_combination.py -M jms --configs configs/JECVar/{}*JEC_{}.py"
            #             " --workdir {}/{} --name {}RunIIJEC{}"
            #         ).format(sample, direction, base_dir, workdir, sample, direction.upper())
            #         + ("" if tagger_argument == "" else " --extra-options \"{}\"".format(tagger_argument))
            #     ]
            #     submitter.add_job("{}RunIIJEC{}".format(sample, direction.upper()), workdir, wrapper_lines)

            # # separate MassScales
            # workdir = "fits_RunII_{}_separateMassScale".format(tagger)
            # wrapper_lines = common_wrapper_lines + [
            #     (
            #         "python jetmass_combination.py -M jms --configs configs/nominal/{}*.py --workdir {}/{}"
            #         " --name {}RunII"
            #     ).format(
            #         sample, base_dir, workdir, sample
            #     )
            #     + ' --extra-options "--separateMassScales --freezeParameters jec_variation {} "'.format(tagger_argument),
            # ]
            # submitter.add_job("{}RunII".format(sample), workdir, wrapper_lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", "-s", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--extract", "-e", action="store_true")
    parser.add_argument("--plots", "-p", action="store_true")
    args = parser.parse_args()

    rhalph_path = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/"

    # separate year fits
    samples = ["Combined", "WJets", "TTBar"]
    years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
    taggers = ["Substructure", "ParticleNet"]

    # date_str = "04-07-23"
    # date_str = "05-07-23_noRobustHesse_range10_Strat1"
    # date_str = "05-07-23_noRobustHesse_range10_Strat1"
    date_str = "05-07-23_noRobustHesse_range15_Strat0_rrange"
    # date_str = "05-07-23_noRobustHesse_range15_Strat0"
    base_dir = "/nfs/dust/cms/user/albrechs/JetMassFits/JMSFits/{}/".format(date_str)
    submitter = JobSubmitter(base_dir)
    common_wrapper_lines = [
        "#!/bin/bash",
        "source /cvmfs/cms.cern.ch/cmsset_default.sh",
        "cd {}/CMSSW_11_3_4/src".format(rhalph_path),
        "eval `scramv1 runtime -sh`",
        "cd {}".format(rhalph_path),
    ]
    add_fits(taggers, samples, years, common_wrapper_lines)
    if args.submit:
        submitter.submit_jobs(debug=args.debug)

    if args.extract:
        for tagger in taggers:
            for fit_collection in ["noJEC", "JECVar", "JECNuis", "separateMassScale"]:
                fits = [
                    "{}/{}".format(fit._workdir, fit._name)
                    for fit in submitter._jobs
                    if all(attribute in fit._workdir for attribute in [tagger, fit_collection])
                ]
                exec_bash(
                    "./extractMassScales.py --fits {} -o fitResults_{}_{}_{}.json".format(
                        " ".join(fits), date_str, tagger, fit_collection
                    ),
                    debug=args.debug,
                )
    if args.plots:
        plotter = {}
        full_years = years  # + ["RunII"]
        for tagger in taggers:
            print(tagger)

            # loading separateMassScale fits
            print("creating plotter for  separate MassScale fits")
            for jms_scale in ["W", "top"]:
                plotter["{}_separateMassScale_{}".format(tagger, jms_scale)] = create_plotter(
                    "fitResults_{}_{}_separateMassScale.json".format(date_str, tagger),
                    full_years,
                    samples,
                    jms_pattern="_{}_".format(jms_scale),
                )

            fit_collections = [
                "{}_JECNuis".format(tagger),
                "{}_JECVar".format(tagger),
                "{}_noJEC".format(tagger),
            ]

            # loading fits with JECNuis
            print("creating plotter for noJEC fits")
            plotter[fit_collections[2]] = create_plotter(
                "fitResults_{}_{}_noJEC.json".format(date_str, tagger), full_years, samples
            )

            # loading fits with JECNuis
            print("creating plotter for JECNuis fits")
            plotter[fit_collections[0]] = create_plotter(
                "fitResults_{}_{}_JECNuis.json".format(date_str, tagger), full_years, samples
            )

            # loading fits without JECNuis
            print("creating plotter for JECVar fits")
            plotter[fit_collections[1]] = create_plotter(
                "fitResults_{}_{}_JECVar.json".format(date_str, tagger), full_years, samples
            )

            # print("creating plotter for JECNVarAsInput fits")
            # plotter[fit_collections[1]].load_fit_results(
            #     json.load(open("fitResults_{}_{}_JECVar.json".format(date_str, tagger))),
            #     [
            #         f"{sample}{year}JEC{jecdir}"
            #         for year in full_years
            #         for jecdir in ["UP", "DOWN"]
            #         for sample in samples
            #     ],
            # )
            plotter[fit_collections[1]].load_JEC_fits(
                json.load(open("fitResults_{}_{}_JECVar.json".format(date_str, tagger)))
            )

            plotter[fit_collections[1]].construct_hists()

        colors = {
            "TTBar": {
                "W": "#f26f6d",
                "top": "#33a02c",
            },
            "WJets": {
                "W": "#8ac4e3",
                "top": "#b2df8a"
            },
            "Combined": {
                "W": "#52a2bf",
                "top": "#4db37b",
            },
            # "TTBar": {
            #     "W": "#1f78b4",
            #     "top": "#33a02c",
            # },
            # "WJets": {
            #     "W": "#a6cee3",
            #     "top": "#b2df8a",
            # },
            # "Combined": {
            #     "W": "#1f78b4",
            #     "top": "#33a02c",
            # },
        }

        color_palette = [
            "tab:blue",
            "tab:orange",
            "tab:green",
            "tab:red",
        ]
        markers = ["o", "^", "v", "x"]
        markersize = 6
        sep_years = {
            "Substructure": ["UL16preVFP", "UL16postVFP", "UL17", "UL18"],
            "ParticleNet": ["UL16preVFP", "UL16postVFP", "UL17", "UL18"],
        }
        # years = {tagger: sep_years[tagger] + ["RunII"] for tagger in taggers}
        years = {tagger: sep_years[tagger] for tagger in taggers}
        for tagger in taggers:
            print(tagger)
            for year in years[tagger]:
                f, ax = setup_ax(10, 7)
                iplot = 0
                legend_alias = {"WJets": r"$W(q\bar{q})$+jets regions", "TTBar": r"$t\bar{t}$ regions"}
                for sample in ["WJets", "TTBar"]:
                    jet_scales = ["W"]
                    # jet_scales = ["W", "top"] if sample == "TTBar" else ["W"]
                    for jms in jet_scales:
                        print(year, sample, jms)
                        jms_plotter = plotter[f"{tagger}_separateMassScale_{jms}"]
                        if f"{sample}{year}" in jms_plotter.dHists.keys():
                            jms_plotter.dHists[f"{sample}{year}"].plot_errorbar(
                                ax=ax,
                                split_uncertainty=True,
                                alpha=0.99,
                                label=f"{legend_alias[sample]} ({jms} scale)",
                                fmt=markers[iplot],
                                linewidth=0.9,
                                color=colors[sample][jms],
                                markersize=markersize,
                            )
                            iplot += 1
                finalize_ax(ax, fname=f"JMSSF_{date_str}_{tagger}/{year}_comparison_sep_scales.pdf", year=year)
                f, ax = setup_ax(10, 7)
                iplot = 0
                for jms in ["W", "top"]:
                    print(year, sample, jms)
                    jms_plotter = plotter[f"{tagger}_separateMassScale_{jms}"]
                    if f"Combined{year}" in jms_plotter.dHists.keys():
                        jms_plotter.dHists[f"Combined{year}"].plot_errorbar(
                            ax=ax,
                            split_uncertainty=True,
                            alpha=0.99,
                            label=f"All regions ({jms} scale)",
                            fmt=markers[iplot],
                            markersize=markersize,
                            linewidth=0.9,
                            color=colors["Combined"][jms],
                        )
                        iplot += 1
                # plotting combined mass scale on same plot
                jms_plotter = plotter[f"{tagger}_JECVar"]
                if f"Combined{year}" in jms_plotter.dHists.keys():
                    jms_plotter.dHists[f"Combined{year}"].plot_errorbar(
                        ax=ax,
                        split_uncertainty=True,
                        alpha=0.99,
                        label="All regions and W and top scale 100% correlated",
                        fmt="v",
                        markersize=markersize,
                        linewidth=0.9,
                        color="#bf77c8",
                        skip_syst_band=True,
                    )

                finalize_ax(ax, fname=f"JMSSF_{date_str}_{tagger}/{year}_Combined_sep_scales.pdf", year=year)

        for tagger in taggers:
            fit_collections = ["{}_noJEC".format(tagger), "{}_JECVar".format(tagger), "{}_JECNuis".format(tagger)]
            fit_info = {
                fit_collections[0]: {
                    "label": r"no JEC applied to $m_{SD}$",
                    "color": "#a6cee3",
                    "marker": "s"
                },
                fit_collections[1]: {
                    "label": "nominal JEC applied",
                    "color": "#bf77c8",
                    "marker": "v",
                },
                fit_collections[2]: {
                    "label": "nominal JEC applied and JEC nuisance parameter",
                    "color": "#33a02c",
                    "marker": "o",
                },
            }
            # plot ratio of JMS over UL17 ones
            for fit_collection in fit_collections:
                jms_plotter = plotter[fit_collection]
                for sample in samples:
                    for denom_year in ["UL16postVFP", "UL17", "UL18"]:
                        fr, axr = setup_ax(10, 7)
                        for iy, year in enumerate(years[tagger]):
                            label_ = f"{fit_info[fit_collection]['label']} {year}"
                            if (
                                f"{sample}{year}" in jms_plotter.dHists.keys()
                                and f"{sample}{denom_year}" in jms_plotter.dHists.keys()
                            ):
                                jms_plotter.dHists[f"{sample}{year}"].plot_ratio(
                                    jms_plotter.dHists[f"{sample}{denom_year}"],
                                    ax=axr,
                                    alpha=0.8,
                                    label=label_,
                                    fmt=fit_info[fit_collection]["marker"],
                                    markersize=markersize,
                                    linewidth=0.9,
                                    color=color_palette[iy],
                                )
                        finalize_ax(
                            axr,
                            fname=f"JMSSF_{date_str}_{tagger}/{sample}_{fit_collection}_ratio_{denom_year}_comparison.pdf",
                            ylabel=f"jms correction ratio (Year/{denom_year})",
                        )

            for sample in samples:
                for year in years[tagger]:
                    f, ax = setup_ax(10, 7)
                    print(sample, year)
                    legend_handlers = []
                    legend_labels = []
                    for fit_collection in fit_collections:
                        jms_plotter = plotter[fit_collection]
                        if f"{sample}{year}" in jms_plotter.dHists.keys():
                            label_ = fit_info[fit_collection]["label"]
                            # if "JECVar" in fit_collection:
                            #     label_ += " (unc. incl. JEC var syst. as shaded bar)"
                            legend_handlers.append(
                                jms_plotter.dHists[f"{sample}{year}"].plot_errorbar(
                                    ax=ax,
                                    split_uncertainty=True,
                                    alpha=0.8,
                                    label=label_,
                                    fmt=fit_info[fit_collection]["marker"],
                                    markersize=markersize,
                                    linewidth=0.9,
                                    color=fit_info[fit_collection]["color"],
                                    # linestyle="--" if "JECVar" in fit_collection else "-",
                                )
                            )
                            legend_labels.append(label_)

                        # if "JECVar" in fit_collection:
                        #     # nominal_color = jms_plotter.dHists[f"{sample}{year}"].color
                        #     extend = True
                        #     if f"{sample}{year}JECDOWN" in jms_plotter.dHists.keys():
                        #         jms_plotter.dHists[f"{sample}{year}JECDOWN"].plot_uncertainty_band(
                        #             ax,
                        #             extend=extend,
                        #             hatch="//",
                        #             alpha=0.4,
                        #             linewidth=0.6,
                        #             label="JEC down variation applied",
                        #             color="#a6cee3",
                        #             edgecolor="k",
                        #         )
                        #         # var_color = jms_plotter.dHists[f"{sample}{year}JECDOWN"].color
                        #     if f"{sample}{year}JECUP" in jms_plotter.dHists.keys():
                        #         jms_plotter.dHists[f"{sample}{year}JECUP"].plot_uncertainty_band(
                        #             ax,
                        #             extend=extend,
                        #             hatch="\\\\",
                        #             linewidth=0.6,
                        #             alpha=0.4,
                        #             label="JEC up variation applied",
                        #             color="#fb9a99",
                        #             edgecolor="k",
                        #         )

                    finalize_ax(
                        ax,
                        fname=f"JMSSF_{date_str}_{tagger}/{sample}_{year}_comparison.pdf",
                        year=year,
                        sort_legend_alphabet=True,
                    )
            for sample in samples:
                f, ax = setup_ax(10, 7)
                for iy, year in enumerate(reversed(sep_years[tagger])):
                    jms_plotter = plotter[f"{tagger}_JECVar"]
                    if f"{sample}{year}" in jms_plotter.dHists.keys():
                        jms_plotter.dHists[f"{sample}{year}"].plot_errorbar(
                            ax=ax,
                            split_uncertainty=False,
                            alpha=0.7,
                            label=f"{year}",
                            fmt=markers[iy],
                            linewidth=0.9,
                            color=color_palette[iy],
                            # marker=markers[iy],
                            markersize=markersize,
                            skip_syst_band=True
                        )
                finalize_ax(ax, fname=f"JMSSF_{date_str}_{tagger}/{sample}_year_comparison.pdf")
