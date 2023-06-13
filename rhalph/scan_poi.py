#!/usr/bin/env python3

import os
import argparse
import json
from copy import deepcopy
import uproot
import numpy as np

pois = {
    "ttbar": [
        "massScale_pt0_eta0_all_PT200",
        "massScale_pt0_eta0_all_PT300",
        "massScale_pt0_eta0_all_PT400",
        "massScale_pt0_eta0_all_PT500",
    ],
    "vjets": [
        "massScale_pt0_eta0_all_PT500",
        "massScale_pt0_eta0_all_PT650",
        "massScale_pt0_eta0_all_PT800"
    ],
    "combined": [
        "massScale_pt0_eta0_all_PT200",
        "massScale_pt0_eta0_all_PT300",
        "massScale_pt0_eta0_all_PT400",
        "massScale_pt0_eta0_all_PT500",
        "massScale_pt0_eta0_all_PT650",
        "massScale_pt0_eta0_all_PT800",
    ],
}


class MassScalePOIScanner(object):
    def __init__(self, args):
        if "Combined" in args.config:
            self.selection = "combined"
        elif "VJets" in args.config or "WJets" in args.config:
            self.selection = "vjets"
        elif "TTBar" in args.config:
            self.selection = "ttbar"

        self.fit_result_path = os.path.abspath(f"{args.config}/{os.path.basename(args.config)}fitResult.json")

        with open(self.fit_result_path, "r") as fit_result_file:
            self.fit_results_copy = deepcopy(
                json.load(fit_result_file)
            )
            fit_result_file.close()

        # self.pois = pois[self.selection] if args.poi == "all" else [args.poi]
        if args.poi == "all":
            self.pois = pois[self.selection]
        elif "*" in args.poi:
            self.pois = []
            for p in self.fit_results_copy.keys():
                if args.poi.replace("*", "") in p:
                    self.pois.append(p)
        else:
            self.pois = [args.poi]
        print(self.pois)        
        self.split_uncs = {
            poi: {"up": [], "down": [], "source": args.freeze_nuisances + [args.lasteffect]}
            for poi in self.pois
        }
        self.poi_minimas = {poi: [] for poi in self.pois}
        self.args = args
        self.n_grid_points = 50

        os.chdir(args.config)

    def exec_bash(self, cmd):
        print(cmd)
        if not self.args.debug:
            os.system(cmd)

    def run(self):
        self.nominal_fit()

        for poi in self.pois:
            self.scan_poi(poi)

        self.dump_split_unc()

    def nominal_fit(self):
        if not os.path.isfile("model_combined.root"):
            os.system("source ./build.sh")

        snapshot_fit_cmd = (
            f"combine -M MultiDimFit -d model_combined.root "
            f"--redefineSignalPOIs  {','.join(self.pois)} "
            f"--setParameterRanges {':'.join('%s=-15,15'%poi for poi in self.pois)} "
            f"--cminDefaultMinimizerStrategy 0 --robustFit 1 "
            f"-n .nominal_fit --saveWorkspace "
            "--saveFitResult "
            f"--preFitValue {1 if self.args.poi == 'r' else 0} "
            # f"--cminDefaultMinimizerStrategy 0 "
        )
        if self.args.force or not os.path.isfile("higgsCombine.nominal_fit.MultiDimFit.mH120.root"):
            self.exec_bash(snapshot_fit_cmd)

    def scan_poi(self, poi):
        # poi_range = self.args.poi_range
        poi_range = None
        if poi_range is None:
            try:
                poi_fitDiag_result = self.fit_results_copy[poi]
                poi_fitDiag_unc = abs(poi_fitDiag_result[1]), abs(poi_fitDiag_result[2])
                # take minimum if all unc greater 0, otherwise take maximum (hoping this will never be not a number)
                poi_max_range = min(*poi_fitDiag_unc) if all(u > 0 for u in poi_fitDiag_unc) else max(*poi_fitDiag_unc)
                poi_range = [poi_fitDiag_result[0] - poi_max_range * 2, poi_fitDiag_result[0] + poi_max_range * 2]
            except BaseException as e:
                print(e)
                print("Could not extract POI range from FitDiagnostics fit results. Sticking with default (-10,10)")
                poi_range = [-10.0, 10.0]

        def perform_multidimfit(step, algo, poi_range=[-10, 10]):
            base_cmd = (
                f"combine -M MultiDimFit -d higgsCombine.nominal_fit.MultiDimFit.mH120.root --snapshotName MultiDimFit "
                # "--skipInitialFit "
                # f"combine -M MultiDimFit -d model_combined.root "
                f"-P {poi}  --setParameterRanges {poi}={','.join(map(str,poi_range))} "
                f"--redefineSignalPOIs {','.join(pois[self.selection])} "
                # f"--algo=grid --points 50 "
                # f"--algo=singles "
                f"--cminDefaultMinimizerStrategy 0 --robustFit 1 "
                f" --floatOtherPOIs 1 "
                f"--preFitValue {1 if poi == 'r' else 0} "
                # f"--cminDefaultMinimizerStrategy 0 "
            )

            if algo == "singles":
                base_cmd += "--algo=singles "
            elif algo == "grid":
                # base_cmd += f"--algo=grid --points {self.n_grid_points} --autoRange 2"
                base_cmd += f"--algo=grid --points {self.n_grid_points} "

            name = f".{step}_{algo}_{poi}"

            freeze = ""
            if step != "nominal":
                freeze += "--freezeParameters "
                freeze_pars_cumsum = []
                for par in self.args.freeze_nuisances:
                    freeze_pars_cumsum.append(par)
                    if par == step:
                        break
                freeze += ",".join(freeze_pars_cumsum)

            if self.args.force or not os.path.isfile(f"higgsCombine{name}.MultiDimFit.mH120.root"):
                self.exec_bash(f"{base_cmd} -n {name} {freeze}")

        # perform singles fits
        perform_multidimfit(step="nominal", algo="singles")
        for nuisance in self.args.freeze_nuisances:
            perform_multidimfit(step=nuisance, algo="singles")

        def get_scan_results(scan_name):
            file_ = uproot.open(f"higgsCombine{scan_name}.MultiDimFit.mH120.root")
            arrays = file_["limit"].arrays([poi, "deltaNLL"])
            try:
                arrays = {k.decode(): v for k, v in arrays.items()}
            except BaseException as e:
                print(e)
            arrays = {k: v.astype(float) for k, v in arrays.items()}
            return arrays

        # todo implement singles call and make scan optional if one wants a plot
        nominal_singles = get_scan_results(f".nominal_singles_{poi}")  # central down up
        freeze_nuisances_singles = {
            nuisance: get_scan_results(f".{nuisance}_singles_{poi}")
            for nuisance in self.args.freeze_nuisances
        }

        def construct_error_composition(errors_list):
            errors_sq = np.array(errors_list) ** 2
            error_comp = np.concatenate([errors_sq[:-1] - errors_sq[1:], errors_sq[-1:]])
            return np.sqrt(error_comp)

        self.poi_minimas[poi].append(nominal_singles[poi][0])
        indep_error_down = [abs(nominal_singles[poi][1]-nominal_singles[poi][0])]
        indep_error_up = [abs(nominal_singles[poi][2]-nominal_singles[poi][0])]
        for freeze_nuisances, results in freeze_nuisances_singles.items():
            self.poi_minimas[poi].append(results[poi][0])
            # self.split_uncs[poi] += list(results[poi][1:])
            # self.split_uncs[poi]["down"].append(abs(results[poi][1]-results[poi][0]))
            # self.split_uncs[poi]["up"].append(abs(results[poi][2]-results[poi][0]))
            indep_error_down.append(abs(results[poi][1]-results[poi][0]))
            indep_error_up.append(abs(results[poi][2]-results[poi][0]))

        self.split_uncs[poi]["down"] = list(
            construct_error_composition(indep_error_down)
        )
        self.split_uncs[poi]["up"] = list(
            construct_error_composition(indep_error_up)
        )

        poi_min = max(
            nominal_singles[poi][0]-1.5*abs(nominal_singles[poi][0]-nominal_singles[poi][1]),
            -10.
        )
        poi_max = min(
            nominal_singles[poi][0]+1.5*abs(nominal_singles[poi][0]-nominal_singles[poi][2]),
            10.
        )

        if args.plots:
            # perform grid scan fits
            perform_multidimfit(step="nominal", algo="grid", poi_range=[poi_min, poi_max])
            for nuisance in self.args.freeze_nuisances:
                perform_multidimfit(step=nuisance, algo="grid", poi_range=[poi_min, poi_max])

            nominal_scan = get_scan_results(f".nominal_grid_{poi}")
            nominal_scan = {k: v[1:] for k, v in nominal_scan.items()}
            freeze_nuisances_scans = {
                nuisance: get_scan_results(
                    f".{nuisance}_grid_{poi}"
                )
                for nuisance in self.args.freeze_nuisances
            }
            freeze_nuisances_scans = {
                name: {k: v[1:] for k, v in scan.items()} for name, scan in freeze_nuisances_scans.items()
            }
            import matplotlib as m
            m.use('Agg')
            import matplotlib.pyplot as plt
            import mplhep as hep
            hep.set_style("CMS")
            import re

            config_pattern = re.compile("(?P<fit>(Combined|TTBar|VJets))(?P<year>UL1[678]{1}((pre|post){1}VFP)*)")
            pattern_match = config_pattern.match(args.config)
            fit_year_dict = pattern_match.groupdict() if pattern_match else {"fit": "", "year": ""}
            
            def legend_par_str(singles):
                return ("POI = "
                        + "${%.3g}" % singles[poi][0]
                        + "^{+%.3g}" % (singles[poi][2]-singles[poi][0])
                        + "_{-%.3g}$" % (singles[poi][0]-singles[poi][1])
                        )

            f, ax = plt.subplots()
            ax.plot(
                nominal_scan[poi],
                2 * nominal_scan["deltaNLL"],
                "k-.",
                label=f"nominal best-fit {legend_par_str(nominal_singles)}",
            )

            alpha_aux_lines = 0.6
            ax.plot(ax.get_xlim(), [1, 1], "k--", alpha=alpha_aux_lines)

            ax.plot(2 * [nominal_singles[poi][1]], [0, 1], "k--", alpha=alpha_aux_lines)
            ax.plot(2 * [nominal_singles[poi][2]], [0, 1], "k--", alpha=alpha_aux_lines)

            errors_up = [nominal_singles[poi][2] - nominal_singles[poi][0]]
            errors_down = [nominal_singles[poi][0] - nominal_singles[poi][1]]

            for name, scan in freeze_nuisances_scans.items():
                line = ax.plot(
                    scan[poi],
                    2 * scan["deltaNLL"],
                    linestyle="-.",
                    label=f"freeze {name.replace('_',' ')} {legend_par_str(freeze_nuisances_singles[name])}",
                )
                errors_up.append(freeze_nuisances_singles[name][poi][2] - freeze_nuisances_singles[name][poi][0])
                errors_down.append(freeze_nuisances_singles[name][poi][0] - freeze_nuisances_singles[name][poi][1])

                line_color = line[0].get_color()
                ax.plot(
                    2 * [freeze_nuisances_singles[name][poi][1]],
                    [0, 1],
                    linestyle="--",
                    alpha=alpha_aux_lines,
                    color=line_color,
                )
                ax.plot(
                    2 * [freeze_nuisances_singles[name][poi][2]],
                    [0, 1],
                    linestyle="--",
                    alpha=alpha_aux_lines,
                    color=line_color,
                )

            error_comp_up = construct_error_composition(errors_up)
            error_comp_down = construct_error_composition(errors_down)
            error_comp_names = [nuisance_name.replace("_", "~") for nuisance_name in self.args.freeze_nuisances]
            error_comp_names.append(self.args.lasteffect)
            error_comp_str = (
                r"$\sigma = "
                + "".join([
                    (
                        r"^{+%.2g}" % (error_comp_up[ierror])
                        + r"_{-%.2g}" % (error_comp_down[ierror])
                        + r"(\mathrm{%s})" % (error_comp_names[ierror])
                    )
                    for ierror in range(len(error_comp_names))
                ])
                + "$"
            )

            ax.text(
                ax.get_xlim()[0]+0.2*np.diff(ax.get_xlim()),
                ax.get_ylim()[1]*0.7,
                error_comp_str, fontsize=15)
            hep.label._exp_label(
                exp="",
                llabel="Private work (CMS data/simulation)",
                year=fit_year_dict["year"],
                ax=ax,
                fontsize=20,
            )
            leg = ax.legend(title=fit_year_dict["fit"]+" fit")
            leg._legend_box.align = "left"
            ax.set_xlabel(poi)
            ax.set_ylabel(r"$-2\Delta \mathrm{log} \mathcal{L}$")
            ax.set_ylim(0, None)
            f.savefig(
                f"likelihood_scan_{poi}_errorcomp_{'_'.join(self.args.freeze_nuisances)}.pdf", bbox_inches="tight"
            )

    def dump_split_unc(self):
        # fit_result = None
        # with open(self.fit_result_path, "r") as json_file:
        #     fit_result = json.load(json_file)
        #     json_file.close()

        # for par, uncs in self.split_uncs.items():
        #     for unc in uncs:
        #         if unc not in fit_result[par]:
        #             fit_result[par].append(unc)
        print("AAA", self.split_uncs)
        json.dump(
            self.split_uncs, open(self.fit_result_path.replace(".json", "SplitUnc.json"), "w"), sort_keys=True, indent=2
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="fit directory")
    parser.add_argument("-P", "--poi", default="massScale_pt0_eta0_all_PT500")
    # parser.add_argument("-r", "--poi-range", nargs="+", default=None)
    parser.add_argument(
        "-F", "--force", action="store_true", help="force to recreate perform fits, even if output already exists!"
    )
    parser.add_argument("-f", "--freeze-nuisances", nargs="+", default=["jec_variation"])
    parser.add_argument(
        "--lasteffect",
        default="Stat.",
        help=(
            "label in plot for the remaining uncertainty. Mostly this will be Stat. since the ultimate superset "
            "of frozen nuisances is mostly all nuisances."
        ),
    )
    parser.add_argument("--plots", action="store_true")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        raise RuntimeError(f"fit directory {args.config} does not exist!")

    scanner = MassScalePOIScanner(args)
    scanner.run()
