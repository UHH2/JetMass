#!/usr/bin/env python

import os,subprocess,sys,threading

class ProgressBar(threading.Timer):
    def run(self):
        while True:
            self.finished.wait(self.interval)
            if self.finished.is_set():
                return
            else:
                self.function(*self.args, **self.kwargs)

def update_progress(iteration,complete):
    barLength = 30
    status = ""
    progress=float(iteration)/float(complete)
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rFits processed: %i/%i [%s] %s"%(int(iteration),int(complete),"#"*block+"-"*(barLength-block),status)
    sys.stdout.write(text)
    sys.stdout.flush()


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
        combine -M FitDiagnostics """+dir+"""_combined.txt --plots --saveShapes
        echo "Fit Done!!!"
        """)
        wrapper.close()
    os.system('chmod u+x '+dir+'/wrapper.sh')

def runFits(models=['UHH_Model']):
    fitProcesses=[]
    finishedProcesses=[]
    for model in models:
        if(not os.path.exists(model)):
            print('Model directory %s does not exists. skipping..')
            continue
        write_wrapper(model)
        command=['env','-i','bash',model+'/wrapper.sh']
        fitProcesses.append(subprocess.Popen(command, shell=None,stdout=open(model+'/log.o','w'),stderr=open(model+'/log.e','w')))

    bar=ProgressBar(0.1,update_progress,args=(len(finishedProcesses),len(fitProcesses)))
    bar.daemon=True
    bar.start()
    finishedProcesses=[p.wait() for p in fitProcesses]
    bar.cancel()


if __name__ == "__main__":
    runFits([sys.argv[1]])
