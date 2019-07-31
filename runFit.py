#!/usr/bin/env python
import os,subprocess

def write_wrapper(dir):
    with open(dir+'/wrapper.sh','w') as wrapper:
        wrapper.write("""
        #!/bin/bash
        #module use -a /afs/desy.de/group/cms/modulefiles/; module load cmssw
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd ../CMSSW_8_1_0/src/
        eval `scramv1 runtime -sh`
        cd ../../rhalphalib/"""+dir+"""
        source build.sh
        combine -M FitDiagnostics """+dir+"""_combined.txt --plots
        echo "Fit Done!!!"
        """)
        wrapper.close()
    os.system('chmod u+x '+dir+'/wrapper.sh')

def runFits(models=['UHH_Model']):
    fitProcesses=[]
    finishedProcesses=[]
    for model in models:
        write_wrapper(model)
        command=['env','-i','bash',model+'/wrapper.sh']
        fitProcesses.append(subprocess.Popen(command, shell=None,stdout=open(model+'/log.o','w'),stderr=open(model+'/log.e','w')))

if __name__ == "__main__":
    runFits()
