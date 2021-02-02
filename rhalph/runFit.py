#!/usr/bin/env python
from __future__ import print_function
import os,subprocess,sys

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
