#!/usr/bin/env python
from __future__ import print_function
import os,subprocess,sys

def write_wrapper(model_dir,pathCMSSW):
    abs_dir=os.path.abspath(model_dir)
    #create wrapper script to setup environment, combine datacards and run the fit in combine
    model_name = model_dir.replace('/','')
    with open(model_dir+'/wrapper.sh','w') as wrapper:
        wrapper.write("""
        #!/bin/bash
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd """+pathCMSSW+"""/src/
        eval `scramv1 runtime -sh`
        cd """+abs_dir+"""
        source build.sh
        combine -M FitDiagnostics """+abs_dir+"""/"""+model_name+"""_combined.txt --plots --saveShapes

        PostFitShapesFromWorkspace -w """+abs_dir+"""/"""+model_name+"""_combined.root -o """+abs_dir+"""/fit_shapes.root --postfit --sampling -f """+abs_dir+"""/fitDiagnostics.root:fit_s
        """)
        wrapper.close()
    os.system('chmod u+x '+model_dir+'/wrapper.sh')

def runFits(models=['UHH_Model']):
    pathCMSSW = os.path.abspath('CMSSW_10_2_13')
    if(not os.path.isdir(pathCMSSW)):
        raise NotADirectoryError('CMSSW with combine is not installed!')
    fitProcesses=[]
    finishedProcesses=[]
    for model in models:
        if(not os.path.exists(model)):
            print('Model directory %s does not exists. skipping..'%model)
            continue
        write_wrapper(model,pathCMSSW)
        #execute wrapper in new environment
        command=['env','-i','bash',model+'/wrapper.sh']
        fitProcesses.append(subprocess.Popen(command, shell=None,stdout=open(model+'/log.o','w'),stderr=open(model+'/log.e','w')))
        print('Fitting', model,'...')

    finishedProcesses=[p.wait() for p in fitProcesses]
    for model in models:
        os.system('tail %s/log.o'%model)

if __name__ == "__main__":
    if(len(sys.argv)<3):
        runFits([sys.argv[1]])
    else:
        runFits([sys.argv[1]],sys.argv[2])
