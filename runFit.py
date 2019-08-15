#!/usr/bin/env python
from __future__ import print_function
import os,subprocess,sys

def write_wrapper(dir,pathCMSSW):
    absdir=os.path.abspath(dir)
    #create wrapper script to setup environment, combine datacards and run the fit in combine
    with open(dir+'/wrapper.sh','w') as wrapper:
        wrapper.write("""
        #!/bin/bash
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd """+pathCMSSW+"""/src/
        eval `scramv1 runtime -sh`
        cd """+absdir+"""
        source build.sh
        combine -M FitDiagnostics """+absdir+"""/"""+dir.replace('/','')+"""_combined.txt --plots --saveShapes
        """)
        wrapper.close()
    os.system('chmod u+x '+dir+'/wrapper.sh')

def runFits(models=['UHH_Model'],pathCMSSW='../CMSSW_8_1_0'):
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
