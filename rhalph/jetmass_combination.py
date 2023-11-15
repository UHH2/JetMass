#!/usr/bin/env python
from __future__ import print_function
import os
import subprocess
import time
from CombineUtils import import_config
from jetmass_scale_fit_utils import extract_fit_results, correlation_coefficient, extract_correlation_matrix
import numpy as np
import fitplotter
from copy import deepcopy
import re
import json
from ROOT import TFile, TH1
TH1.AddDirectory(False)


def get_hists(fname):
    """Getting hists from a ROOT file.
    Expects a list of TDirectory objects at the basedir, which holds a list of hists each.
    """
    file_ = TFile(fname)
    dirs = [key.GetName() for key in file_.GetListOfKeys()]
    hists = {
        hist_dir: {
            key.GetName(): file_.Get("{}/{}".format(hist_dir, key.GetName()))
            for key in file_.Get(hist_dir).GetListOfKeys()
        }
        for hist_dir in dirs
    }
    return hists


def dump_config(config, path):
    open(path, "w").write(
        json.dumps(config, sort_keys=False, indent=2)
    )


class JetMassCombination(object):
    def __init__(self, workdir, configs, name="FullRunII", year="RunII", mode="jms", config_update="{}", job_index=""):
        self.workdir = workdir
        self._job_index = job_index
        self._config_dicts = {}
        self._models = []
        self._config_update = config_update
        self._configs = {}
        self._split_postfitshapes = True
        self.configs = configs
        self._mode = mode
        self.extra_options = ""
        self._combination_dir = "{}/{}/".format(self.workdir, name)
        if not os.path.exists(self._combination_dir):
            os.makedirs(self._combination_dir)

        # create bigger workdir on dust for combineTool workflows
        basename = os.path.basename(self.workdir)
        self._dust_dir = "/nfs/dust/cms/user/albrechs/JetMassFits/CombineToolWorkdir/{}".format(basename)
        print(self._dust_dir)
        if not os.path.exists(self._dust_dir):
            os.makedirs(self._dust_dir)
        if not os.path.exists("{}/CombineToolWorkdir".format(self._combination_dir)):
            os.system("ln -sf {} {}/CombineToolWorkdir".format(self._dust_dir, self._combination_dir))
        for model in self._models:
            model_config = self._config_dicts[model]
            model_config["ModelDir"] = "{}/{}".format(self.workdir, model)
            model_config_json_path = "{}/{}config.json".format(self._combination_dir, model)
            if os.path.isfile(model_config_json_path):
                self._config_dicts[model] = import_config(model_config_json_path)
            else:
                dump_config(model_config, model_config_json_path)

        self.combined_config = deepcopy(self._config_dicts[self._models[-1]])
        self.combined_config.update(
            {
                "ModelName": name,
                "ModelDir": self._combination_dir,
                "year": year,
            }
        )
        dump_config(self.combined_config, "{}/config.json".format(self._combination_dir))

    @property
    def configs(self):
        return self._configs

    @configs.setter
    def configs(self, config_list):
        for config in config_list:
            config_dict = import_config(config)
            config_dict.update(eval(self._config_update))
            model_name = config_dict["ModelName"]
            model_name += self._job_index
            config_dict["ModelName"] = model_name
            self._models.append(model_name)
            self._config_dicts[model_name] = config_dict
            self._configs[model_name] = config

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, w):
        self._workdir = os.path.abspath(w)
        if not os.path.exists(self._workdir):
            os.makedirs(self._workdir)

    def build_models(self):
        model_processes = []
        for model_name in self._models:
            cmd = "python jetmass.py -M {} ".format(self._mode)
            cmd += "--workdir {} ".format(self.workdir)
            cmd += self.extra_options
            if self._job_index != "" and "--job_index" not in self.extra_options:
                cmd += " --job_index {} ".format(self._job_index)
            cmd += " --build "
            cmd += " --skipTemplatePlots "
            if self._config_update != "{}":
                cmd += " --config-update \"{}\" ".format(self._config_update.replace('"', '\\"'))
            cmd += "{} ".format(self.configs[model_name])
            print(cmd)
            with open("{}/{}.log".format(self.workdir, model_name), "wb") as log:
                model_processes.append(subprocess.Popen(cmd, stdout=log, stderr=log, shell=True))
        print("building individual models...")
        while len(model_processes) > 0:
            model_processes[-1].wait()
            if model_processes[-1].returncode == 0:
                model_processes.pop()
            else:
                [p.terminate() for p in model_processes]
                raise RuntimeError("One of the processes had an issue. Check the logs!")
            time.sleep(2)

    def combine_models(self):
        for model in self._models:
            os.system("cp {}/{}/{}.root {}".format(self.workdir, model, model, self._combination_dir))
            # overwrite config with updated per year config
            # important for unfolding, since the list of signal processes is updated depending on the choice of gen-level bins
            os.system(
                "cp {}/{}/config.json {}/{}config.json".format(self.workdir, model, self._combination_dir, model)
            )

        with open("{}/build.sh".format(self._combination_dir), "wb") as build_script:
            # combine cards
            cmd = ""
            combination_cmd = "combineCards.py "
            for model in self._models:
                cmd += "combineCards.py {}={}/{}/model_combined.txt > {}_renamed.txt\n".format(
                    model, self.workdir, model, model
                )
                combination_cmd += "{}={}/{}/model_combined.txt ".format(model, self.workdir, model)
            cmd += combination_cmd + " > {}/model_combined.txt\n".format(self._combination_dir)
            build_script.write(cmd)

            # text2workspace
            if self._mode != "unfolding":
                cmd = "text2workspace.py {}/model_combined.txt\n".format(self._combination_dir)
                build_script.write(cmd)

        templates_model_name = self._models[-1]
        # take wrapper from last added config/model and edit that
        os.system(
            "cp {}/{}/wrapper.sh {}/wrapper_template.sh".format(
                self.workdir, templates_model_name, self._combination_dir
            )
        )

        # tf_parameter_postfit = {}
        # for model in self._models:
        #     fit_result = extract_fit_results(self._config_dicts[model], return_result=True)
        #     if not fit_result:
        #         print("Fit of model {} failed. Keeping initial TF params at 1.0.".format(model))
        #         continue
        #     for param, values in fit_result.items():
        #         if "tf_" in param:
        #             tf_parameter_postfit[param] = values[0]

        def edit_line(line_):
            sections = line_.split(" ")
            resulting_sections = []
            for section in sections:
                if templates_model_name in section:
                    section = "".join([self._combination_dir, section.split(templates_model_name)[-1]])
                resulting_sections.append(section)
            return " ".join(resulting_sections)

        with open("{}/wrapper.sh".format(self._combination_dir), "wb") as wrapper:
            with open("{}/wrapper_template.sh".format(self._combination_dir), "r") as template:
                for line in template:
                    if "PostFitShapesFromWorkspace" in line and self._split_postfitshapes:
                        continue
                    if "combineCards.py" in line and self._mode == "unfolding":
                        continue
                    else:
                        if templates_model_name in line:
                            line = edit_line(line)

                        # if "combine" in line:
                            # set_parameters_addition = ",".join(
                            #     "{}={}".format(param, value) for param, value in tf_parameter_postfit.items()
                            # )
                            # if "--setParameters" in line:
                            #     line = line.replace(
                            #         "--setParameters ", "--setParameters {},".format(set_parameters_addition)
                            #     )
                            # else:
                            #     line = line.strip() + " --setParameters {} \n ".format(set_parameters_addition)

                        if "combine" in line and "--freezeParameters" in line:
                            # we have to freeze Parameters, first lets get the list from the template (i.e. last model)
                            pattern = r"--freezeParameters\s+([\w,]+)"
                            match = re.search(pattern, line.strip())
                            freezeParams = set()
                            if match:
                                freezeParams = set(match.group(1).split(","))
                            else:
                                raise RuntimeError("Could not extract list of parameters to freeze")
                            # does any of them have year in their name?
                            # If so we need to get the ones from the other wrappers as well
                            template_year = self._config_dicts[self._models[-1]]["year"]
                            correlated_params = [param for param in freezeParams if template_year in param]
                            if len(correlated_params) > 0:
                                for model_name in self._models[:-1]:
                                    with open("{}/{}/wrapper.sh".format(self.workdir, model_name), "r") as sub_wrapper:
                                        for line_sub_wrapper in sub_wrapper.readlines():
                                            match = re.search(pattern, line_sub_wrapper)
                                            if match:
                                                for param in match.group(1).strip().split(","):
                                                    freezeParams.add(param)
                            print("full list of params to freeze", freezeParams)
                            line = line.replace(
                                re.search(pattern, line.strip()).group(),
                                "--freezeParameters {}".format(",".join(freezeParams))
                            )
                    if "model_combined.txt" in line:
                        for model in self._models:
                            wrapper.write(line.replace("model_combined", "{}_renamed".format(model)))
                    wrapper.write(line)
                    if self._mode == "unfolding" and "cd " in line:
                        wrapper.write("source build.sh\n")

        os.system("rm {}/wrapper_template.sh".format(self._combination_dir))

    def run_combined_fit(self):
        os.system("chmod +x {}/wrapper.sh".format(self._combination_dir))
        os.system("bash {}/wrapper.sh".format(self._combination_dir))

    def combine_fit_shapes(self):
        # if self._split_postfitshapes:
        self.combine_fit_shapes_from_split()
        # else:
        self.combine_fit_shapes_from_precompiled()

    def update_config(self):
        template_config = json.load(open("{}/{}config.json".format(self._combination_dir, self._models[-1]), "r"))
        signal_templates = template_config["channels"]["WJetsPt575"]["signal"]
        self.combined_config["genbins"] = template_config["genbins"]
        self.combined_config["histLocation"] = template_config["histLocation"]
        for channel in self.combined_config["channels"]:
            if "WJetsMatched" in self.combined_config["channels"][channel]["samples"]:
                self.combined_config["channels"][channel]["samples"].remove("WJetsMatched")
                self.combined_config["channels"][channel]["signal"].remove("WJetsMatched")
                for signal_template in signal_templates:
                    self.combined_config["channels"][channel]["samples"].append(signal_template)
                    self.combined_config["channels"][channel]["signal"].append(signal_template)
        dump_config(self.combined_config, "{}/config.json".format(self._combination_dir))

    def combine_fit_shapes_from_split(self):
        model_processes = []
        fit_shapes_files = []
        for model_name in self._models:
            fit_shape_filename = "{}/{}_fit_shapes.root".format(self._combination_dir, model_name)
            fit_shapes_files.append(fit_shape_filename)
            cmd = "PostFitShapesFromWorkspace  -w {}/{}_renamed.root ".format(self._combination_dir, model_name)
            cmd += "--postfit --sampling "
            cmd += "--output {} ".format(fit_shape_filename)
            cmd += "-f {}/fitDiagnostics.root:fit_s &> /dev/null".format(self._combination_dir)
            print(cmd)
            with open("{}/{}.log".format(self.workdir, model_name), "wb") as log:
                model_processes.append(subprocess.Popen(cmd, stdout=log, stderr=log, shell=True))
        print("getting postfit shapes...")
        while len(model_processes) > 0:
            model_processes[-1].wait()
            if model_processes[-1].returncode == 0:
                model_processes.pop()
            else:
                [p.terminate() for p in model_processes if p.returncode]
                raise RuntimeError("One of the processes had an issue. Check the logs!")
            time.sleep(2)

        for fname in fit_shapes_files:
            tfile = TFile(fname, "UPDATE")
            for dir_key in tfile.GetListOfKeys():
                hist_dir_name = dir_key.GetName()
                hist_dir = tfile.Get(hist_dir_name)
                for hist_key in hist_dir.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    if "qcd_" not in hist_name:
                        continue
                    hist_path = "{}/{}".format(hist_dir_name, hist_name)
                    hist = tfile.Get(hist_path)
                    tfile.cd(hist_dir_name)
                    hist.SetTitle("qcd")
                    hist.SetName("qcd")
                    hist.Write()
                    tfile.cd("../")
                    tfile.Delete("{};*".format(hist_path))
            tfile.Close()

        os.system("hadd {}/fit_shapes_split_years.root {}".format(self._combination_dir, " ".join(fit_shapes_files)))

    def combine_fit_shapes_from_precompiled(self):
        split_fit_shapes_filename = "{}/fit_shapes_split_years.root".format(self._combination_dir)
        combined_fit_shapes_filename = "{}/fit_shapes.root".format(self._combination_dir)

        hists = get_hists(split_fit_shapes_filename)

        common_hist_dirs = [
            hist_dir.replace(self._models[-1] + "_", "") for hist_dir in hists.keys() if self._models[-1] in hist_dir
        ]
        common_hist_keys = {hist_dir: set() for hist_dir in common_hist_dirs}
        for model_name in self._models:
            for hist_dir in common_hist_dirs:
                for hist_id in hists["{}_{}".format(model_name, hist_dir)].keys():
                    common_hist_keys[hist_dir].add(hist_id)

        combined_hists = {}
        for hist_dir, hist_keys in common_hist_keys.items():
            this_hists = {(lambda f: "qcd" if "qcd_" in f else f)(key): None for key in hist_keys}
            # beware hists.keys() != hist_keys e.g. [qcd_UL17, qcd_UL18] != [qcd]
            for hist_key in this_hists.keys():
                for model_name in self._models:
                    model_hist_dir = hists["{}_{}".format(model_name, hist_dir)]
                    lookup_key = hist_key
                    if hist_key == "qcd":
                        qcd_pattern = re.compile("qcd_*")
                        qcd_match_key = None
                        for key in model_hist_dir.keys():
                            qcd_match = qcd_pattern.search(key)
                            if qcd_match is not None:
                                qcd_match_key = key
                                break
                        if qcd_match_key is None:
                            raise RuntimeError(
                                "could not find qcd hist, but was expecting one. {}/{}".format(model_name, hist_dir)
                            )
                        qcd_key = qcd_match_key
                        if qcd_key not in hist_keys:
                            raise RuntimeError(
                                "Found qcd key, that is not in list of expected keys! {} [{}]".format(
                                    qcd_key, ",".join(hist_keys)
                                )
                            )
                        lookup_key = qcd_key
                    if lookup_key not in model_hist_dir:
                        print("Could not find {} in histdir".format(lookup_key))
                        continue
                    if this_hists[hist_key] is None:
                        this_hists[hist_key] = model_hist_dir[lookup_key]
                    else:
                        this_hists[hist_key].Add(model_hist_dir[lookup_key])
            combined_hists[hist_dir] = this_hists

        combined_fit_shapes_file = TFile(combined_fit_shapes_filename, "RECREATE")
        for hist_dir, hists in combined_hists.items():
            combined_fit_shapes_file.mkdir(hist_dir)
            combined_fit_shapes_file.cd(hist_dir)
            for hist_key, hist in hists.items():
                hist.SetName(hist_key)
                hist.Write()
            combined_fit_shapes_file.cd("..")
        combined_fit_shapes_file.Close()
        # os.system("cp {}/fit_shapes.root {}/fit_shapes_split.root".format(self._combination_dir, self._combination_dir))
        # os.system(
        #     "mv {}/combined_fit_shapes.root {}/fit_shapes.root".format(self._combination_dir, self._combination_dir)
        # )

    def dump_correlation_coefficient(self):
        pois = ["r_{}".format(genbin) for genbin in self.combined_config.get("genbins", [])]
        corr_coeff = correlation_coefficient(
            "{}/fitDiagnostics.root".format(self._combination_dir), pois
        )
        np.save("{}/corr_coeff.npy".format(self._combination_dir), corr_coeff)

    def impacts(self, step):
        cmd = "ulimit -s unlimited;"
        if step == "fits":
            cmd += "cd {} && combineTool.py -M Impacts -d {}/model_combined.root ".format(
                self._dust_dir, self._combination_dir
            )
            cmd += " --doInitialFit --robustFit 1 -m 0 --exclude rgx{qcdparam*} "
            cmd += " --cminDefaultMinimizerStrategy 0 "
            cmd += " && "
            cmd += "combineTool.py -M Impacts -d {}/model_combined.root ".format(self._combination_dir)
            cmd += "-m 0 --robustFit 1 --doFits --exclude rgx{qcdparam*} --parallel 100 --job-mode condor "
            cmd += "--cminDefaultMinimizerStrategy 0"
        elif step == "plots":
            cmd += "cd {} && combineTool.py -M Impacts -d {}/model_combined.root ".format(
                self._dust_dir, self._combination_dir
            )
            cmd += " --exclude rgx{qcdparam*} "
            cmd += " -m 0 -o impacts.json"
            cmd += "&& "
            cmd += r'for ipt in 0 1 2 3; do echo $ipt ;printf "%s\n" 0 1 2 3 |'
            cmd += "xargs -I{} -P 5 -n1 plotImpacts.py -i impacts.json -o impacts_r_ptgen${ipt}_msdgen{}"
            cmd += " --POI r_ptgen${ipt}_msdgen{} --per-page 35 "
            cmd += " --translate ${CMSSW_BASE}/../POI_rename.json --height 600;done"
        print(cmd)
        # os.system(cmd)


if __name__ == "__main__":
    import argparse
    # import logging
    # logging.basicConfig()

    # logger = logging.getLogger()
    # logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default="jetmass_combination", help="workdir for combined jetmass fits")
    parser.add_argument("-M", "--mode", choices=["unfolding", "jms", "default"],
                        help="choose what fit should be run. 'unfolding' performs MaxLik. unfolding, jms measures "
                        "JetMassScale variations,'default' measures signal cross section.",
                        required=True)
    parser.add_argument("--configs", nargs="+", required=True)
    parser.add_argument("--config-update", type=str, default="{}")
    parser.add_argument("--name", default="FullRunII")
    parser.add_argument("--correlationCoeffcients", action="store_true")
    parser.add_argument("--year", default="RunII")
    parser.add_argument("--extra-options", default="")
    parser.add_argument("--job_index", type=str, default="")
    parser.add_argument("--sumgenbins", action="store_true")
    parser.add_argument("--impacts", choices=["none", "fits", "plots"], default="none")
    parser.add_argument("--justplots", action="store_true", help="just redo the plots.")
    args = parser.parse_args()

    asimov = "prefitasimov" in args.extra_options.lower()

    JMS_Combination = JetMassCombination(
        args.workdir, args.configs, name=args.name, year=args.year,
        mode=args.mode, config_update=args.config_update,
        job_index=args.job_index
    )
    JMS_Combination.extra_options = args.extra_options
    print(args.impacts)
    if not args.justplots and args.impacts == "none":
        JMS_Combination.build_models()

        JMS_Combination.combine_models()
        JMS_Combination.run_combined_fit()

        JMS_Combination.combine_fit_shapes()
    if args.impacts == "none":
        do_postfit = extract_fit_results(JMS_Combination.combined_config)
        if do_postfit:
            extract_correlation_matrix(JMS_Combination.combined_config)
            if args.correlationCoeffcients:
                JMS_Combination.dump_correlation_coefficient()

        unfolding = args.mode == "unfolding"
        if unfolding:
            JMS_Combination.update_config()

        fitplotter.plot_fit_result(
            JMS_Combination.combined_config,
            plot_total_sig_bkg=False,
            do_postfit=do_postfit,
            use_config_samples=unfolding,
            pseudo_data=False,
            prefit_asimov=asimov,
            unfolding=unfolding,
            fit_shapes_root="fit_shapes.root",
            sum_genbins=args.sumgenbins,
        )
        # do_postfit=True
        for c_name_, config in JMS_Combination._config_dicts.items():
            this_config = deepcopy(config)
            this_config["ModelDir"] = JMS_Combination._combination_dir
            this_config["TFSuffix"] = this_config["year"]
            fitplotter.plot_qcd_bernstein(this_config, do_3d_plot=False)
            fitplotter.plot_qcd_fail_parameters(this_config)

        # print(
        #     "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/pretty_postfit.py {} --mctruth --year {}".format(JMS_Combination._combination_dir, args.year)
        # )
        # os.system(
        #     "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python/pretty_postfit.py {} --mctruth --year {}".format(JMS_Combination._combination_dir, args.year)
        # )
    else:
        JMS_Combination.impacts(args.impacts)
