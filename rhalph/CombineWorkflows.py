#!/usr/bin/env python
from __future__ import print_function
import os,argparse,ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


import sys
sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/python')
import fitplotter
import cms_style
cms_style.extra_text="Preliminary Simulation"
cms_style.cms_style()
cms_logo = False
global silence
silence = False

def exec_bash(command='echo "hello world"',debug=False):
    global silence
    if(not silence):
        print(command)
    if(not debug):
        os.system(command)
    return """%s\n"""%command

def get_parameters(args,query, workspace):
    w = ROOT.TFile(args.workspace,"READ").Get('w')
    allVars = w.allVars().contentsString().split(',')
    parameters = []
    for var in allVars:
        if '_In' in var:
            continue
        if query in var:
            parameters.append(var)
    return parameters


rMIN=-1.5
rMAX=1.5

class CombineWorkflows:
    def __init__(self):
        # print('Nothing to do here. This class is just a wrapper for some worflow methods using combine')
        self.methods = [func for func in dir(CombineWorkflows) if (callable(getattr(CombineWorkflows, func)) and '__' not in func )]        
        self.methods.remove('write_wrapper')
        # self.methods.remove('method')
        self.externToys = False
        self.justplots = False
        self.skipplots = True
        self._poi = ""
        self.workspace = ""
        self.workers = 1
        self._freezeParameters = ""
        self.lumi = 41.8
        self.name = ""
        self.seed = 123456
        self.toys = 50
        self.algo = "saturated"
        self.extraOptions = ""
        self.job_index = 0
        self.bernsteinOrders = "2:2"
        self.rhalphdir = os.getcwd()#"/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph"
        self.toysOptions = "--toysFrequentist"

        self.workspace = ""

        self._method = ""

        def dummyMethod(debug=True):
            raise BaseException("You have not selected a CombineWorkflow method! Choose from: "+", ".join(self.methods))
        self.combineString = dummyMethod

    
    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, w):
        if(isinstance(w,str)):
            self._workspace = os.path.abspath(w)
            self.modeldir = self.workspace.replace(self.workspace.split('/')[-1],'')
        elif(isinstance(w,list)):
            self._workspace = np.array([os.path.abspath(iw) for iw in w],dtype=object)
            self.modeldir = np.array([os.path.dirname(iw) for iw in self.workspace],dtype=object)
            
    @property
    def freezeParameters(self):
        return self._freezeParameters

    @freezeParameters.setter
    def freezeParameters(self,pars):
        if(isinstance(pars,str)):
            if(pars == ""):
                self._pars = ""
            elif('--freezeParameters ' in pars):
                self._pars = pars
            else:
                self._pars = '--freezeParameters ' + pars
        elif(isinstance(pars,list)):
            self._poi = '--freezeParameters ' + ','.join(pars)
            return
        else:
            raise TypeError("freezeParameters must be either string (',' as delimiter) or list!\nYou provided a ",type(pars))

        
    @property
    def POI(self):
        return self._poi

    @POI.setter
    def POI(self,pois):
        if(isinstance(pois,str)):
            if(pois.lower() == "r"):
                self._poi = ""
            elif('--redefineSignalPOIs ' in pois):
                self._poi = pois
            else:
                self._poi = '--redefineSignalPOIs ' + pois
        elif(isinstance(pois,list)):
            self._poi = '--redefineSignalPOIs ' + ','.join(pois)
            return
        else:
            raise TypeError("POI must be either string (',' as delimiter) or list!\nYou provided a ",type(pois))

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self,method_name):
        self._method = method_name
        self.combineString = getattr(self,method_name)
        
    def write_wrapper(self):
        global silence
        silence = True
        pathCMSSW = os.path.realpath(self.rhalphdir + '/CMSSW_10_2_13')
        if(not os.path.isdir(self.modeldir)):
            os.makedirs(self.modeldir)
        wrapper_name = self.modeldir+'/wrapper.sh'
        with open(wrapper_name,'w') as wrapper:
            wrapper.write("#!/bin/bash\n")
            wrapper.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            wrapper.write("cd "+pathCMSSW+"/src/\n")
            wrapper.write("eval `scramv1 runtime -sh`\n")
            wrapper.write(self.combineString(debug=True))
            wrapper.close()
        os.system('chmod u+x '+ wrapper_name)
        silence = False


    def diagnostics(self, debug=True):
        command_string = """#FitDiagnostics Workflow\n"""
        command_string += exec_bash("cd " + self.modeldir + "\n",debug)
        command_string += exec_bash("source build.sh\n",debug)
        command_string += exec_bash("combine -M FitDiagnostics {WORKSPACE} {POI} --saveShapes {EXTRA}".format(WORKSPACE=self.workspace,POI=self.POI,EXTRA=self.extraOptions),debug)
        command_string += exec_bash("PostFitShapesFromWorkspace -w {WORKSPACE} -o {MODELDIR}fit_shapes.root --postfit --sampling -f {MODELDIR}fitDiagnostics.root:fit_s".format(WORKSPACE=self.workspace,MODELDIR=self.modeldir),debug)
        return command_string

    def GOF(self,debug=True,merge=False):
        command_string = """#GOF test\n"""
        command_string += exec_bash("cd " + self.modeldir,debug)
        command_string += exec_bash("source build.sh",debug)
        

        if(not self.justplots and not merge):
            command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,ALGO=self.algo),debug)
            if(self.externToys):
                command_string += exec_bash("combine -M GenerateOnly -d {WORKSPACE} {POI} -m 0 -t {NTOYS} --toysFrequentist --saveToys -n {NAME} --seed {SEED}".format(WORKSPACE=self.workspace,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys),debug)
                command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} --toysFrequentist  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,ALGO=self.algo),debug)
            else:
                command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPT} --seed {SEED}".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPT=self.toysOptions,ALGO=self.algo),debug)
        command_string += exec_bash('python {RHALPHDIR}/CombinePlotter.py --method plot_gof_result --parameter "higgsCombine{NAME}Baseline.GoodnessOfFit.mH0.root;higgsCombine.{NAME}GoodnessOfFit.mH0.{SEED}.root;{ALGO};{LUMI}"'.format(RHALPHDIR=self.rhalphdir,NAME=self.name,SEED=self.seed,ALGO=self.algo,LUMI=self.lumi),debug)
        return command_string

    def FTest(self, debug=True):
        command_string = "#FTest\n"
        qcd_fit = "qcdmodel" in self.modeldir
        import glob,json
        configs =[os.path.abspath(json_path) for json_path in  glob.glob(self.modeldir + ("../" if qcd_fit else "") + "../*.json")]
        alt_model_dirs = [os.path.dirname(config)+'/'+json.load(open(config))["ModelName"]+("/qcdmodel/" if qcd_fit else"/") for config in configs]
        alt_model_dirs.remove(self.modeldir)

        command_string += exec_bash("cd "+ self.modeldir,debug)
        command_string += exec_bash("source build.sh",debug)        

        command_string += exec_bash("combine -M MultiDimFit -d {WORKSPACE} --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI),debug)
        command_string += exec_bash("combine  -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,ALGO=self.algo),debug)
        command_string += exec_bash("combine  -M GenerateOnly -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} -t {NTOYS} {TOYSOPTIONS} {EXTRA} --saveToys -n \"{NAME}\" --seed {SEED}".format(WORKSPACE=self.workspace,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,EXTRA=self.extraOptions),debug)
        command_string += exec_bash("combine  -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS} --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,ALGO=self.algo),debug)

        
        # using the snapshot ensures that one can load the externally produced toys and still have non-zero extArgs in the ModelConfig (by default only the nuisances are saved into the toy snapshot -> specified by ModelConfig)
        # command_string += exec_bash("combine -M MultiDimFit -d {WORKSPACE} --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI),debug)
        # command_string += exec_bash("combine higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,ALGO=self.algo),debug)
        # command_string += exec_bash("combine higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -M GenerateOnly  -m 0 {POI} -t {NTOYS} {TOYSOPTIONS} --saveToys -n \"{NAME}\" --seed {SEED}".format(WORKSPACE=self.workspace,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions),debug)
        # command_string += exec_bash("combine higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} --toysFrequentist  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,ALGO=self.algo),debug)

        
        # command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,ALGO=self.algo),debug)
        # command_string += exec_bash("combine -M GenerateOnly -d {WORKSPACE} {POI} -m 0 -t {NTOYS} {TOYSOPTIONS} --saveToys -n \"{NAME}\" --seed {SEED}".format(WORKSPACE=self.workspace,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions),debug)
        # command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} --toysFrequentist  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,ALGO=self.algo),debug)

        command_string += exec_bash("",debug)
        
        for model_dir in alt_model_dirs:
            command_string += exec_bash("",debug)
            command_string += exec_bash("cd "+ model_dir,debug)
            command_string += exec_bash("source build.sh",debug)

            # not sure what to do here. One can not load the snapshot created above since then we'd just redo the command above. However just using the toysfile generated above seems to yield reasonable GOF values.
            # just create another snapshot with the alternative Model dummdumm 
            command_string += exec_bash("combine -M MultiDimFit -d *_combined.root --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\"".format(WORKSPACE=self.workspace,FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI),debug)
            command_string += exec_bash("combine higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\"".format(FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,ALGO=self.algo),debug)
            command_string += exec_bash("combine higgsCombine{NAME}Snapshot.MultiDimFit.mH0.root --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS} --toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}".format(FREEZEPARAMS=self.freezeParameters,EXTRA=self.extraOptions,NAME=self.name,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,ALGO=self.algo,BASEMODELDIR=self.modeldir),debug)


            
        # print(self.modeldir)
        # print(alt_model_dirs)
        
        command_string += exec_bash("cd " + self.modeldir + "\n",debug)

        return command_string
if(__name__=='__main__'):
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--justplots',action="store_true")
    parser.add_argument('--skipplots',action="store_true")
    parser.add_argument('--debug',action="store_true")
    parser.add_argument('--method',type=str,choices=globals()['CombineWorkflows']().methods,required=True)
    parser.add_argument('--POI',default="r")
    parser.add_argument('--workspace','-w',default = 'WJetsOneScale_combined.root')
    parser.add_argument('--workers',default=5)
    parser.add_argument('--freezeParameters',default="None")
    parser.add_argument('--lumi',default=41.8)
    parser.add_argument('-n','--name',default="")
    parser.add_argument('--seed',default="123456")
    parser.add_argument('-t','--toys',default=50)
    parser.add_argument('--algo',default="saturated")
    parser.add_argument('--extra',default="",help='pass extra arguments/options to combine commands')
    parser.add_argument('--job_index',default=0,type=int)
    parser.add_argument('--bernsteinOrders',default="2:2")
    parser.add_argument('--externToys',action="store_true")
    parser.add_argument('--rhalphdir',type=str,default="/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph")
    
    args = parser.parse_args()

    if(not os.path.isfile(args.workspace)):
        raise IOError('Could not find workspace file')
    
    args.model_dir=os.path.abspath('/'.join(args.workspace.split('/')[:-1])+'/') if '/' in args.workspace else ''

    if(args.job_index>0):
        args.seed = str(int(args.seed)+args.job_index)
        print('jobindex!=0. resetting seed to initial_seed+job_index:',args.seed)

    print('workspace',args.workspace)
    print('model_dir',args.model_dir)
    print('using method',args.method)
    method = getattr(globals()['CombineWorkflows'](),args.method)
    command_string = method(args,args.debug)
    
    if(args.debug):
        print()
        print()
        print(command_string)

    
    
    

    
