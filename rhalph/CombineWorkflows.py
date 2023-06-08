#!/usr/bin/env python
from __future__ import print_function
import os
import argparse
import ROOT
import rhalphalib as rl
import numpy as np
import sys
# import cms_style

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/rhalphalib')
sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/python')

rl.util.install_roofit_helpers()

# cms_style.extra_text = "Preliminary Simulation"
# cms_style.cms_style()


def _RooFitResult_massScales(self):
    result = []
    for p in self.floatParsFinal():
        if "massScale" in p.GetName():
            result.append([p.GetName(), p.getVal(), p.getErrorHi(), p.getErrorLo()])
    return np.array(result, dtype=object)


ROOT.RooFitResult.massScales = _RooFitResult_massScales


cms_logo = False
global silence
silence = False


def exec_bash(command='echo "hello world"', debug=False):
    global silence
    if not silence:
        print(command)
    if not debug:
        os.system(command)
    return """%s\n""" % command


def get_parameters(args, query, workspace):
    w = ROOT.TFile(args.workspace, "READ").Get("w")
    allVars = w.allVars().contentsString().split(",")
    parameters = []
    for var in allVars:
        if "_In" in var:
            continue
        if query in var:
            parameters.append(var)
    return parameters


rMIN = -1.5
rMAX = 1.5


class CombineWorkflows(object):
    def __init__(self, build_only=False):
        # print('Nothing to do here. This class is just a wrapper for some worflow methods using combine')
        self.methods = [
            func for func in dir(CombineWorkflows) if (callable(getattr(CombineWorkflows, func)) and "__" not in func)
        ]
        self.methods.remove('write_wrapper')
        # self.methods.remove('method')
        self.externToys = False
        self.justplots = False
        self.skipplots = True
        self._poi = ""
        self.POIRange = (-10, 10)
        self.altmodel = None
        self.workers = 1
        self._freezeParameters = []
        self.lumi = 41.8
        self.name = ""
        self.seed = 123456
        self.toys = 50
        self.algo = "saturated"
        self.extraOptions = ""
        self.job_index = 0
        self.rhalphdir = os.getcwd()
        self.toysOptions = "--toysFrequentist"
        self.combineCMSSW = self.rhalphdir + '/CMSSW_10_2_13'
        self.modeldir = ""
        self._workspace = 'model_combined.root'
        self._method = ""
        self._build_prefix = "#" if build_only else ""
        self._skip_plotting = False

        def dummyMethod(debug=True):
            raise BaseException("You have not selected a CombineWorkflow method! Choose from: "+", ".join(self.methods))
        self.combineString = dummyMethod

    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, w):
        if sys.version_info.major < 3:
            if isinstance(w, unicode): # noqa
                w = str(w)

        if isinstance(w, str):
            self._workspace = os.path.abspath(w)
            self.modeldir = self.workspace.replace(self.workspace.split("/")[-1], "")
        elif isinstance(w, list):
            self._workspace = np.array([os.path.abspath(iw) for iw in w], dtype=object)
            self.modeldir = np.array([os.path.dirname(iw) for iw in self.workspace], dtype=object)

    @property
    def altmodel(self):
        return self._altmodel

    @altmodel.setter
    def altmodel(self, w):
        if isinstance(w, str):
            self._altmodel = os.path.abspath(w)
        elif w is None:
            self._altmodel = None

    @property
    def freezeParameters(self):
        if len(self._freezeParameters) == 0:
            return ""
        return "--freezeParameters " + ",".join(self._freezeParameters)

    @freezeParameters.setter
    def freezeParameters(self, pars):
        print("adding pars to freeze list:", pars)
        if isinstance(pars, str):
            self._freezeParameters += pars.split(",")
        elif isinstance(pars, list):
            self._freezeParameters += pars
        else:
            raise TypeError(
                "freezeParameters must be either string (',' as delimiter) or list!\nYou provided a ", type(pars)
            )

    @property
    def POI(self):
        return self._poi

    @POI.setter
    def POI(self, pois):
        if isinstance(pois, str):
            if pois.lower() == "r":
                self._poi = ""
            elif "--redefineSignalPOIs " in pois:
                self._poi = pois
            else:
                self._poi = "--redefineSignalPOIs " + pois
        elif isinstance(pois, list):
            self._poi = "--redefineSignalPOIs " + ",".join(pois)
            range_str = "=%i,%i:" % (self.POIRange[0], self.POIRange[1])
            self._poi += " --setParameterRanges " + range_str.join(pois) + range_str[:-1]
            return
        else:
            raise TypeError("POI must be either string (',' as delimiter) or list!\nYou provided a ", type(pois))

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method_name):
        self._method = method_name
        self.combineString = getattr(self, method_name)

    def write_wrapper(self, append=False):
        global silence
        silence = True
        # pathCMSSW = os.path.realpath(self.combineCMSSW)

        if not os.path.isdir(self.modeldir):
            os.makedirs(self.modeldir)
        wrapper_name = self.modeldir + "/wrapper.sh"

        append_wrapper = os.path.isfile(wrapper_name) and append

        with open(wrapper_name, "a" if append_wrapper else "w") as wrapper:
            if not append_wrapper:
                wrapper.write("#!/bin/bash\n")
                wrapper.write("ulimit -s unlimited\n")
            # wrapper.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            # wrapper.write("cd "+pathCMSSW+"/src/\n")
            # wrapper.write("eval `scramv1 runtime -sh`\n")
            wrapper.write(self.combineString(debug=True))
            wrapper.close()
        os.system("chmod u+x " + wrapper_name)
        silence = False

    def unfolding(self, debug=True):
        command_string = """#Unfolding Workflow\n"""
        command_string += exec_bash("cd " + self.modeldir + "\n", debug)
        import json

        configs = json.load(open(self.modeldir + "/config.json", "r"))

        bin_signal_strenght_constructor = "[1,0.1,10.0]"

        genbins = configs["genbins"]

        # adding lines to setup env variables helpful for later workflows
        command_string += exec_bash(
            "export POIS=({})".format(" ".join(["r_{GENBIN}".format(GENBIN=genbin) for genbin in genbins])),
            debug,
        )
        command_string += exec_bash('if [[ $1 == "env" ]];then return 0 ; fi', debug)

        regions_mapping = " ".join(
            [
                "{NAME}={NAME}.txt".format(NAME=channel_name + region_name)
                for channel_name, channel_config in configs["channels"].items()
                for region_name in channel_config["regions"]
            ]
        )
        command_string += exec_bash(
            "combineCards.py {REGIONS} > {DATACARD} \n".format(
                REGIONS=regions_mapping, DATACARD=self.workspace.replace(".root", ".txt")
            ),
            debug,
        )

        # unfolding_bins = configs["unfolding_bins"]
        # genbins = [
        #     "ptgen%i_msdgen%i" % (iptgen, imsdgen)
        #     for iptgen in range(len(unfolding_bins["ptgen"]) - 1)
        #     for imsdgen in range(len(unfolding_bins["msdgen"]) - 1)
        # ]

        # don't use first and last msd genbin as signal but treat as background
        # pt_uflow, msd_uflow = genbins[0].split("_")
        # pt_oflow, msd_oflow = genbins[-1].split("_")

        # def is_flow(bin_name):
        #     pt_bin, msd_bin = bin_name.split("_")
        #     # if pt_bin == pt_uflow or pt_bin == pt_oflow:
        #     #     return True
        #     # if msd_bin == msd_uflow or msd_bin == msd_oflow:
        #     if msd_bin == msd_uflow:
        #         return True

        # genbins = [bin_name for bin_name in genbins if not is_flow(bin_name)]

        regularization_strength = configs.get("regularizationStrength", 0.1)
        datacard = self.workspace.replace(".root", ".txt")
        command_string += exec_bash('sed -i "s/kmax/kmax * # /g" {DATACARD} \n'.format(DATACARD=datacard), debug)
        command_string += exec_bash('echo "" >> {DATACARD}'.format(DATACARD=datacard), debug)
        if configs.get("regularization", [""]) != [""]:
            command_string += exec_bash(
                'echo "# SVD regularization penalty terms" >> {DATACARD}'.format(DATACARD=datacard), debug
            )

        rs = ["r_"+b for b in genbins]

        def regul_bin(b, inf_repl):
            if isinstance(b, int):
                return float(b)
            elif isinstance(b, float):
                return b
            if "inf" in b.lower():
                b = b.lower().replace("inf", inf_repl)
            return float(b)
        pt_cutoff = str(configs.get("pt_cutoff", 1400))
        msd_cutoff = str(configs.get("msd_cutoff", 100))
        pt_edges = np.array([regul_bin(bin_, pt_cutoff) for bin_ in configs["unfolding_bins"]["ptgen"]])
        msd_edges = np.array([regul_bin(bin_, msd_cutoff) for bin_ in configs["unfolding_bins"]["msdgen"]])
        n_pt_bins = len(pt_edges) - 1
        n_msd_bins = len(msd_edges) - 1

        pt_centers = (pt_edges[1:] + pt_edges[:-1]) / 2
        msd_centers = (msd_edges[1:] + msd_edges[:-1]) / 2

        pt_widths = pt_edges[1:] - pt_edges[:-1]
        msd_widths = msd_edges[1:] - msd_edges[:-1]

        pt_width_average = pt_widths.sum() / n_pt_bins
        msd_width_average = msd_widths.sum() / n_msd_bins

        n_bins = {1: n_msd_bins, 2: n_pt_bins}
        centers = {1: msd_centers, 2: pt_centers}
        widths = {1: msd_widths, 2: pt_widths}
        Delta = {1: msd_width_average, 2: pt_width_average}
        # Delta = {1:msd_widths.max(), 2:pt_widths.max()}
        global_bins = np.arange(n_bins[1] * n_bins[2]).reshape((n_bins[2], n_bins[1])).T

        def delta(d, i, j):
            if i > j:
                i_ = j
                j = i
                i = i_
            return centers[d][j] - centers[d][i]

        def w(d, i):
            return widths[d][i]

        def get_sub_c(d, uniform=False):
            sub_c = np.zeros((n_bins[d], n_bins[d]))
            for i in range(n_bins[d]):
                if i == 0:
                    if uniform:
                        sub_c[i, 0] = -1
                        sub_c[i, 1] = 1
                    else:
                        sub_c[i, 0] = -Delta[d] / delta(d=d, i=0, j=1)
                        sub_c[i, 1] = +Delta[d] / delta(d=d, i=0, j=1)
                elif i == (n_bins[d]-1):
                    if uniform:
                        sub_c[i, i] = -1
                        sub_c[i, i - 1] = 1
                    else:
                        sub_c[i, i - 1] = +Delta[d] / delta(d=d, i=0, j=1)
                        sub_c[i, i] = -Delta[d] / delta(d=d, i=0, j=1)
                else:
                    if uniform:
                        sub_c[i, i - 1] = 1
                        sub_c[i, i] = -2
                        sub_c[i, i + 1] = 1
                    else:
                        sub_c[i, i - 1] = Delta[d] ** 2 / (
                            (delta(d=d, i=i, j=i - 1) + delta(d=d, i=i, j=i + 1)) * delta(d=d, i=i, j=i - 1)
                        )
                        sub_c[i, i] = (
                            -Delta[d] ** 2
                            / ((delta(d=d, i=i, j=i - 1) + delta(d=d, i=i, j=i + 1)))
                            * (1 / delta(d=d, i=i, j=i - 1) + 1 / delta(d=d, i=i, j=i + 1))
                        )
                        sub_c[i, i + 1] = Delta[d] ** 2 / (
                            (delta(d=d, i=i, j=i - 1) + delta(d=d, i=i, j=i + 1)) * delta(d=d, i=i, j=i + 1)
                        )
            return sub_c

        def get_C(dimensions=[1, 2], uniform=False, density=True):
            C = []

            for d in dimensions:
                d_other = {1: 2, 2: 1}[d]
                sub_c = get_sub_c(d, uniform)
                for i_other_dim in range(n_bins[d_other]):
                    if d == 2:
                        ws = np.ones_like(sub_c)/(w(1, i_other_dim)*w(2, slice(None,)))
                    elif d == 1:
                        ws = np.ones_like(sub_c)/(w(1, slice(None,))*w(2, i_other_dim))
                    if uniform:
                        sub_c_final = sub_c
                    else:
                        sub_c_final = np.multiply(sub_c, ws) if density else sub_c

                    for row in sub_c_final:
                        row_final = np.zeros((n_bins[1]*n_bins[2]))
                        for ic, c in enumerate(row):
                            if d == 1:
                                row_final[global_bins[ic, i_other_dim]] = c
                            elif d == 2:
                                row_final[global_bins[i_other_dim, ic]] = c
                        C.append(row_final)
            return np.array(C)

        constr_index = 0
        uniform_binning = configs.get("uniformGenbins", "False") == "True"
        scale_to_binwidth = configs.get("scaleGenBinWidth", "False") == "True"

        def constrain_snippet(iconst, r):
            snippet = "constr{} constr ".format(iconst)
            for ic, c in enumerate(r):
                if c == 0:
                    continue
                sign = "+" if c > 0 else "-"

                def number(n):
                    return ("%.0f*" if uniform_binning else "%.5f*") % abs(c)

                snippet += "{}{}{}".format(sign, number(c), rs[ic])
                # snippet += "{}{:.0f}*{}".format(sign, abs(c), rs[ic])
            snippet += " delta[{}]".format(regularization_strength)
            return snippet

        dims = [{"msd": 1, "pt": 2, "": -1}[d] for d in configs.get("regularization", [""])]
        if -1 not in dims:
            C = get_C(dims, uniform_binning, scale_to_binwidth)
            for row in C:
                command_string += exec_bash('echo "{CONSTR}" >> {DATACARD}\n'.format(
                    CONSTR=constrain_snippet(constr_index, row),
                    DATACARD=datacard
                ))
                constr_index += 1

        POMAPS = " ".join(
            [
                "--PO map='.*{GENBIN}.*:r_{GENBIN}{PARCONSTRUCT}'".format(
                    GENBIN=genbin, PARCONSTRUCT=bin_signal_strenght_constructor
                )
                for genbin in genbins
            ]
            # +
            # [
            #     "--PO map='.*ptgen.*:r[1.0,0.1,10.0]'"  # TODO use this as signal-strength modifier for all genbins?
            #     # Note 02-02-23 This will not work since this PO map wwill overwrite the unfolding rs
            # ]
        )

        command_string += exec_bash(
            "text2workspace.py -m 0 --X-allow-no-background -o {WORKSPACE} {DATACARD} "
            "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel {POMAPS}\n".format(
                WORKSPACE=self.workspace, DATACARD=self.workspace.replace(".root", ".txt"), POMAPS=POMAPS
            ),
            debug,
        )

        poi_defaults = ",".join(["r_{GENBIN}=1".format(GENBIN=genbin) for genbin in genbins])
        command_string += exec_bash(
            (
                "{BUILDPREFIX}combine -M FitDiagnostics -d {WORKSPACE} --saveShapes -n '' "
                # "--cminDefaultMinimizerStrategy 0 --robustFit 1 "
                "--robustHesse 1 "
                "{FREEZEPARAMS}"
            ).format(
                BUILDPREFIX=self._build_prefix,
                FREEZEPARAMS=self.freeze_Parameters,
                WORKSPACE=self.workspace,
            ),
            debug,
        )
        if not self._skip_plotting:
            command_string += exec_bash(
                "{BUILDPREFIX}PostFitShapesFromWorkspace -w {WORKSPACE} --output {MODELDIR}fit_shapes.root --postfit "
                "--sampling  -f {MODELDIR}fitDiagnostics.root:fit_s".format(
                    BUILDPREFIX=self._build_prefix,
                    WORKSPACE=self.workspace, MODELDIR=self.modeldir
                ),
                debug,
            )

        # command_string += exec_bash(
        #     "combine -M MultiDimFit --setParameters={PARAMETERS} -t -1 -m 0 {WORKSPACE} "
        #     "-n .Asimov --saveWorkspace\n".format(
        #         PARAMETERS=poi_defaults, WORKSPACE=self.workspace
        #     ),
        #     debug,
        # )

        # command_string += exec_bash(
        #     (
        #         "combine -M FitDiagnostics -d {WORKSPACE} --saveShapes -n '.PseudoAsimov_10kEvents' "
        #         "--cminDefaultMinimizerStrategy 0 -t -1 --expectSignal 1 --X-rtd TMCSO_PseudoAsimov=10000 "
        #         "{FREEZEPARAMS}"
        #     ).format(
        #         FREEZEPARAMS=""
        #         if len(freeze_Parameters) == 0
        #         else ("--freezeParameter " + ",".join(freeze_Parameters)),
        #         WORKSPACE=self.workspace,
        #     ),
        #     debug,
        # )

        # command_string += exec_bash(
        #     "PostFitShapesFromWorkspace -w {WORKSPACE} --output {MODELDIR}fit_shapes_pseudoasimov_10k.root --postfit "
        #     "--sampling  -f {MODELDIR}fitDiagnostics.PseudoAsimov_10kEvents.root:fit_s".format(
        #         WORKSPACE=self.workspace, MODELDIR=self.modeldir
        #     ),
        #     debug,
        # )

        # command_string += exec_bash(
        #     (
        #         "combine -M FitDiagnostics -d {WORKSPACE} --saveShapes -n '.PseudoAsimov_100kEvents' "
        #         "--cminDefaultMinimizerStrategy 0 -t -1 --expectSignal 1 --X-rtd TMCSO_PseudoAsimov=100000 "
        #         "{FREEZEPARAMS}"
        #     ).format(
        #         FREEZEPARAMS=""
        #         if len(freeze_Parameters) == 0
        #         else ("--freezeParameter " + ",".join(freeze_Parameters)),
        #         WORKSPACE=self.workspace,
        #     ),
        #     debug,
        # )

        # command_string += exec_bash(
        #     "PostFitShapesFromWorkspace -w {WORKSPACE} --output {MODELDIR}fit_shapes.root --postfit --sampling "
        #     "-f {MODELDIR}fitDiagnostics.PseudoAsimov_100kEvents.root:fit_s".format(
        #         WORKSPACE=self.workspace, MODELDIR=self.modeldir
        #     ),
        #     debug,
        # )

        command_string += exec_bash(
            "#combine -M MultiDimFit --setParameters={PARAMETERS}  -m 0 {WORKSPACE} -n .Data \n".format(
                PARAMETERS=poi_defaults, WORKSPACE=self.workspace
            ),
            debug,
        )

        poi_string = "--redefineSignalPOIs " + ",".join(["r_{GENBIN}".format(GENBIN=genbin) for genbin in genbins])
        command_string += exec_bash(
            "#combine -M FitDiagnostics -d {WORKSPACE} "
            "{POI} --setParameters={POIDEFAULTS} --saveShapes -n '' "
            "--cminDefaultMinimizerStrategy 0".format(
                WORKSPACE=self.workspace,
                POI=poi_string, POIDEFAULTS=poi_defaults
            ),
            debug,
        )
        command_string += exec_bash(
            "#PostFitShapesFromWorkspace -w {WORKSPACE} --output {MODELDIR}fit_shapes.root --postfit --sampling "
            "-f {MODELDIR}fitDiagnostics.root:fit_s".format(
                WORKSPACE=self.workspace, MODELDIR=self.modeldir
            ),
            debug,
        )

        return command_string

    def diagnostics(self, debug=True):
        command_string = """#FitDiagnostics Workflow\n"""
        command_string += exec_bash("cd " + self.modeldir + "\n", debug)
        command_string += exec_bash("source build.sh\n", debug)
        command_string += exec_bash(
            "combine -M FitDiagnostics {WORKSPACE} {POI} --saveShapes {EXTRA} -n '' "
            "--cminDefaultMinimizerStrategy 0 {FREEZE}".format(
                WORKSPACE=self.workspace, POI=self.POI, EXTRA=self.extraOptions,
                FREEZE=self.freezeParameters
            ),
            debug,
        )
        command_string += exec_bash(
            "PostFitShapesFromWorkspace -w {WORKSPACE} --output {MODELDIR}fit_shapes.root --postfit "
            "--sampling -f {MODELDIR}fitDiagnostics.root:fit_s".format(
                WORKSPACE=self.workspace, MODELDIR=self.modeldir
            ),
            debug,
        )
        return command_string

    def GOF(self, debug=True, merge=False):
        command_string = """#GOF test\n"""
        command_string += exec_bash("cd " + self.modeldir, debug)
        command_string += exec_bash("source build.sh", debug)

        if not self.justplots and not merge:
            command_string += exec_bash(
                'combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} '
                '{EXTRA} -n "{NAME}Baseline"'.format(
                    WORKSPACE=self.workspace,
                    FREEZEPARAMS=self.freezeParameters,
                    EXTRA=self.extraOptions,
                    NAME=self.name,
                    POI=self.POI,
                    ALGO=self.algo,
                ),
                debug,
            )
            if self.externToys:
                command_string += exec_bash(
                    "combine -M GenerateOnly -d {WORKSPACE} {POI} -m 0 -t {NTOYS} --toysFrequentist --saveToys -n "
                    "{NAME} --seed {SEED}".format(
                        WORKSPACE=self.workspace, NAME=self.name, POI=self.POI, SEED=self.seed, NTOYS=self.toys
                    ),
                    debug,
                )
                command_string += exec_bash(
                    'combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA}'
                    ' -n "{NAME}" -t {NTOYS} --toysFrequentist  '
                    '--toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}'.format(
                        WORKSPACE=self.workspace,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        SEED=self.seed,
                        NTOYS=self.toys,
                        ALGO=self.algo,
                    ),
                    debug,
                )
            else:
                command_string += exec_bash(
                    'combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} '
                    '-n "{NAME}" -t {NTOYS} {TOYSOPT} --seed {SEED}'.format(
                        WORKSPACE=self.workspace,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        SEED=self.seed,
                        NTOYS=self.toys,
                        TOYSOPT=self.toysOptions,
                        ALGO=self.algo,
                    ),
                    debug,
                )
        command_string += exec_bash(
            'python {RHALPHDIR}/CombinePlotter.py --method plot_gof_result '
            '--parameter "higgsCombine{NAME}Baseline.GoodnessOfFit.mH0.root;'
            'higgsCombine.{NAME}GoodnessOfFit.mH0.{SEED}.root;{ALGO};{LUMI}"'.format(
                RHALPHDIR=self.rhalphdir, NAME=self.name, SEED=self.seed, ALGO=self.algo, LUMI=self.lumi
            ),
            debug,
        )
        return command_string

    def FTestBatch(self, debug=True):
        command_string = "#FTest\n"
        qcd_fit = "qcdmodel" in self.modeldir
        import glob
        import json

        configs = [
            os.path.abspath(json_path)
            for json_path in glob.glob(self.modeldir + ("../" if qcd_fit else "") + "../*.json")
        ]
        alt_model_dirs = [
            os.path.dirname(config) + "/" + json.load(open(config))["ModelName"] + ("/qcdmodel/" if qcd_fit else "/")
            for config in configs
        ]
        alt_model_dirs.remove(self.modeldir)

        command_string += exec_bash("cd " + self.modeldir, debug)
        command_string += exec_bash("source build.sh", debug)

        GOF_extra = self.extraOptions  # + (" --fixedSignalStrength 1 " if "r" in self._freezeParameters else "")

        # # using snapshot
        # combine -M MultiDimFit -d workspace.root --saveWorkspace -m 0 -n '.snapshot' --setParameters r=1
        # --toysFrequentist --cminDefaultMinimizerStrategy 2 --cminDefaultMinim# izerTolerance 0.01
        # --seed 42 --freezeParameters r

        command_string += exec_bash(
            'combine -M MultiDimFit   -d {WORKSPACE} --saveWorkspace -m 0 {POI} {EXTRA} '
            '-n "{NAME}Snapshot" --seed {SEED} --trackParameters r'.format(
                WORKSPACE=self.workspace,
                EXTRA=self.extraOptions,
                NAME=self.name,
                POI=self.POI,
                SEED=self.seed,
            ),
            debug,
        )

        # extract r_bestfit to be able to throw toys from exact signal-bestfit:
        # this is necessary since combine apparently loads a snapshot but then resets r to the default
        # (or given) value of --expectSignal when throwing toys -> setting parameters to postfit but
        # setting r=0 by default.
        # to avoid throwing toys on s+b postfit nuisances but with s=0 we extract r here and
        # set it manually when throwing toys
        command_string += exec_bash(
            ('r_bestfit=$(python -c "from __future__ import print_function;import uproot;f=uproot.open(\\"'
             'higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root\\");print(f[\\"limit\\"]'
             '[\\"trackedParam_r\\"].array()[0])")').format(
                NAME=self.name,
                SEED=self.seed,
            ),
            debug,
        )
        # # snapshot + observed
        # combine -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit
        # --algo saturated --setParameters r=1 --freezeParameters r --cminDefaultMinimizerStrategy 2
        # --cminDefaultMinimizerTolerance 0.01 --seed 42 --fixedSignalStrength 1
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit --trackParameters r -m 0 {POI} --seed {SEED} '
            '{FREEZEPARAMS} {EXTRA} -n "{NAME}Baseline" --algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                ALGO=self.algo,
                EXTRA=GOF_extra,
            ),
            debug,
        )
        # # snapshot + generate-only
        # combine -M GenerateOnly -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -t 1
        # --toysFrequentist --seed 42 --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01
        # --setParameters r=1 --freezeParameters r --saveToys  -n '.snapshot2gen'
        command_string += exec_bash(
            'combine -M GenerateOnly  -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} '
            '--expectSignal $r_bestfit --trackParameters r '
            '{EXTRA} -n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  --saveToys  '.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                NTOYS=self.toys,
                TOYSOPTIONS=self.toysOptions,
                EXTRA=self.extraOptions,
            ),
            debug,
        )
        # combine -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit
        # --algo saturated -t 1 --toysFrequentist --setParameters r=1 --cminDefaultMinimizerStrategy 2
        # --cminDefaultMinimizerTolerance 0.01 --toysFile toys.root --seed 42 -n '.snapshot2gen' --freezeParameters r
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} '
            '{EXTRA} -n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root '
            '--expectSignal $r_bestfit --trackParameters r '
            '--algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                NTOYS=self.toys,
                TOYSOPTIONS=self.toysOptions,
                ALGO=self.algo,
                EXTRA=GOF_extra,
            ),
            debug,
        )

        # standalone
        # command_string += exec_bash(
        #     'combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
        #     '-n "{NAME}Baseline" --algo={ALGO}'.format(
        #         WORKSPACE=self.workspace,
        #         NAME=self.name,
        #         FREEZEPARAMS=self.freezeParameters,
        #         POI=self.POI,
        #         SEED=self.seed,
        #         ALGO=self.algo,
        #         EXTRA=GOF_extra,
        #     ),
        #     debug,
        # )
        # command_string += exec_bash(
        #     'combine -M GenerateOnly  -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n "{NAME}" '
        #     '-t {NTOYS} {TOYSOPTIONS}  --saveToys  '.format(
        #         WORKSPACE=self.workspace,
        #         NAME=self.name,
        #         FREEZEPARAMS=self.freezeParameters,
        #         POI=self.POI,
        #         SEED=self.seed,
        #         NTOYS=self.toys,
        #         TOYSOPTIONS=self.toysOptions,
        #         EXTRA=self.extraOptions,
        #     ),
        #     debug,
        # )
        # command_string += exec_bash(
        #     'combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
        #     '-n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root '
        #     '--algo={ALGO}'.format(
        #         WORKSPACE=self.workspace,
        #         NAME=self.name,
        #         FREEZEPARAMS=self.freezeParameters,
        #         POI=self.POI,
        #         SEED=self.seed,
        #         NTOYS=self.toys,
        #         TOYSOPTIONS=self.toysOptions,
        #         ALGO=self.algo,
        #         EXTRA=GOF_extra,
        #     ),
        #     debug,
        # )

        command_string += exec_bash("", debug)

        for model_dir in alt_model_dirs:
            command_string += exec_bash("", debug)
            command_string += exec_bash("cd " + model_dir, debug)
            command_string += exec_bash("source build.sh", debug)

            # using snapshot
            command_string += exec_bash(
                'combine -M MultiDimFit   -d *_combined.root --saveWorkspace -m 0 {POI} {EXTRA} '
                '-n "{NAME}Snapshot" --seed {SEED}'.format(
                    EXTRA=self.extraOptions,
                    NAME=self.name,
                    POI=self.POI,
                    SEED=self.seed,
                ),
                debug,
            )
            command_string += exec_bash(
                'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
                '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
                '-n "{NAME}Baseline" --algo={ALGO}'.format(
                    NAME=self.name,
                    FREEZEPARAMS=self.freezeParameters,
                    POI=self.POI,
                    SEED=self.seed,
                    ALGO=self.algo,
                    EXTRA=GOF_extra,
                ),
                debug,
            )
            command_string += exec_bash(
                'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
                '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
                '-n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  '
                '--expectSignal $r_bestfit --trackParameters r '
                '--toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}'.format(
                    NAME=self.name,
                    FREEZEPARAMS=self.freezeParameters,
                    POI=self.POI,
                    SEED=self.seed,
                    NTOYS=self.toys,
                    TOYSOPTIONS=self.toysOptions,
                    ALGO=self.algo,
                    EXTRA=GOF_extra,
                    BASEMODELDIR=self.modeldir,
                ),
                debug,
            )

            # standalone
            # command_string += exec_bash(
            #     'combine -M GoodnessOfFit -d *_combined.root -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA}'
            #     ' -n "{NAME}Baseline" --algo={ALGO}'.format(
            #         NAME=self.name,
            #         FREEZEPARAMS=self.freezeParameters,
            #         POI=self.POI,
            #         SEED=self.seed,
            #         ALGO=self.algo,
            #         EXTRA=GOF_extra,
            #     ),
            #     debug,
            # )
            # command_string += exec_bash(
            #     'combine -M GoodnessOfFit -d *_combined.root -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA}'
            #     ' -n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  '
            #     '--toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}'.format(
            #         NAME=self.name,
            #         FREEZEPARAMS=self.freezeParameters,
            #         POI=self.POI,
            #         SEED=self.seed,
            #         NTOYS=self.toys,
            #         TOYSOPTIONS=self.toysOptions,
            #         ALGO=self.algo,
            #         EXTRA=GOF_extra,
            #         BASEMODELDIR=self.modeldir,
            #     ),
            #     debug,
            # )

        command_string += exec_bash("cd " + self.modeldir + "\n", debug)

        return command_string

    def FTest(self, debug=True):
        command_string = "#FTest\n"
        print("MODELDIR", self.modeldir)
        if not debug:
            os.chdir(self.modeldir)
        command_string += exec_bash("cd " + self.modeldir, debug)
        # command_string += exec_bash("source build.sh",debug)

        GOF_extra = self.extraOptions + (" --fixedSignalStrength 1 " if "r" in self._freezeParameters else "")
        # # using snapshot
        command_string += exec_bash(
            'combine -M MultiDimFit   -d {WORKSPACE} --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}Snapshot" --seed {SEED}'.format(
                WORKSPACE=self.workspace,
                FREEZEPARAMS=self.freezeParameters,
                EXTRA=self.extraOptions,
                NAME=self.name,
                POI=self.POI,
                SEED=self.seed,
            ),
            debug,
        )
        # # snapshot + observed
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}Baseline" --algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                ALGO=self.algo,
                EXTRA=GOF_extra,
            ),
            debug,
        )
        # # snapshot + generate-only
        command_string += exec_bash(
            'combine -M GenerateOnly  -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  --saveToys'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                NTOYS=self.toys,
                TOYSOPTIONS=self.toysOptions,
                EXTRA=self.extraOptions,
            ),
            debug,
        )
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}" -t {NTOYS} {TOYSOPTIONS}  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root '
            '--algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                NTOYS=self.toys,
                TOYSOPTIONS=self.toysOptions,
                ALGO=self.algo,
                EXTRA=GOF_extra,
            ),
            debug,
        )

        command_string += exec_bash("", debug)

        # using snapshot
        command_string += exec_bash(
            'combine -M MultiDimFit   -d {WORKSPACE} --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}SnapshotAltModel" --seed {SEED}'.format(
                WORKSPACE=self.altmodel,
                FREEZEPARAMS=self.freezeParameters,
                EXTRA=self.extraOptions,
                NAME=self.name,
                POI=self.POI,
                SEED=self.seed,
            ),
            debug,
        )
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}SnapshotAltModel.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}BaselineAltModel" --algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                ALGO=self.algo,
                EXTRA=GOF_extra,
            ),
            debug,
        )
        command_string += exec_bash(
            'combine -M GoodnessOfFit -d higgsCombine{NAME}SnapshotAltModel.MultiDimFit.mH0.{SEED}.root '
            '--snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} '
            '-n "{NAME}AltModel" -t {NTOYS} {TOYSOPTIONS}  '
            '--toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}'.format(
                NAME=self.name,
                FREEZEPARAMS=self.freezeParameters,
                POI=self.POI,
                SEED=self.seed,
                NTOYS=self.toys,
                TOYSOPTIONS=self.toysOptions,
                ALGO=self.algo,
                EXTRA=GOF_extra,
                BASEMODELDIR=self.modeldir,
            ),
            debug,
        )

        command_string += exec_bash("cd " + self.modeldir + "\n", debug)

        return command_string

    def FastScan(self, debug=True):
        command_string = "#FastScan\n"
        command_string += "combineTool.py -M FastScan -w {WORKSPACE}:w".format(WORKSPACE=self.workspace)
        return command_string

    def FastScanMassScales(self, debug=True):
        command_string = "#FastScan\n"
        command_string += "combineTool.py -M FastScan -w {WORKSPACE}:w --match massScale".format(
            WORKSPACE=self.workspace
        )
        return command_string


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--justplots", action="store_true")
    parser.add_argument("--skipplots", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--method", type=str, choices=globals()["CombineWorkflows"]().methods, required=True)
    parser.add_argument("--POI", default="r")
    parser.add_argument("--workspace", "-w", default="WJetsOneScale_combined.root")
    parser.add_argument("--altmodel", "-a", default=None)
    parser.add_argument("--workers", default=5)
    parser.add_argument("--freezeParameters", default=None)
    parser.add_argument("--lumi", default=41.8)
    parser.add_argument("-n", "--name", default="")
    parser.add_argument("--seed", default="1234567")
    parser.add_argument("-t", "--toys", default=50)
    parser.add_argument("--algo", default="saturated")
    parser.add_argument("--extra", default="", help="pass extra arguments/options to combine commands")
    parser.add_argument("--job_index", default=0, type=int)
    parser.add_argument("--externToys", action="store_true")
    parser.add_argument(
        "--rhalphdir",
        type=str,
        default="/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.workspace):
        raise IOError("Could not find workspace file")

    args.modeldir = os.path.abspath("/".join(args.workspace.split("/")[:-1]) + "/") if "/" in args.workspace else ""
    if args.job_index > 0:
        args.seed = str(int(args.seed) + args.job_index)
        print("jobindex!=0. resetting seed to initial_seed+job_index:", args.seed)

    # print('workspace',args.workspace)
    # print('model_dir',args.modeldir)
    # print('using method',args.method)
    cw = CombineWorkflows()
    # setting up CombineWorkflow
    # (this is written with python2 in mind. So the property decorators defined above need to be "recreated" here)
    cw.method = args.method
    cw.POI = "" if args.POI == "r" else ("--redefineSignalPOIs" + args.POI)
    cw.workspace = os.path.abspath(args.workspace)
    cw.altmodel = os.path.abspath(args.altmodel)
    cw.freezeParameters = "" if args.freezeParameters == "" else ("--freezeParameters " + args.freezeParameters)
    cw.name = args.name
    cw.seed = args.seed
    cw.toys = args.toys
    cw.algo = args.algo
    cw.extra = args.extra
    cw.job_index = args.job_index
    cw.externToys = args.externToys
    cw.rhalphdir = args.rhalphdir
    cw.modeldir = args.modeldir

    method = getattr(cw, args.method)
    print(cw)
    command_string = method(args.debug)
    if args.debug:
        print()
        print()
        print(command_string)
