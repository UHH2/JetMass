#!/usr/bin/env python
from __future__ import print_function
import argparse, os, sys, subprocess, time

from multiprocessing import Pool, Value

def submit_job_unpack(args):
    return submit_job(*args)

def submit_job(instance,index):
    instance.submit_job(index)

def init(args):
    global n_finished
    n_finished = args

class JobManager(object):
    def __init__(self, options, workdir):
        self.workdir = workdir
        if(not os.path.isdir(self.workdir)):
            os.makedirs(self.workdir)

        self.worker = options.workers
        self.selection = options.selection
        self.macro = options.macro

        self.finished = [False]*self.worker
        

    def submit_job(self,job_index):
        nice = True
        if(nice):
            subprocess.call(['nice','-n','3',"./"+self.macro,self.selection,str(job_index),str(self.worker)], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
        else:
            subprocess.call(["./"+self.macro,self.selection,str(job_index),str(self.worker)], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
        global n_finished
        with n_finished.get_lock():
            n_finished.value += 1
        return 1

    def start(self):
        global n_finished
        n_finished = Value('i',0)

        print('starting pool with %i workers'%self.worker)

        pool = Pool(processes = self.worker, initializer = init, initargs = (n_finished,))

        result = pool.map_async(submit_job_unpack,[(self,i) for i in range(self.worker)], chunksize = 1)
        
        while not result.ready():
            with n_finished.get_lock():
                done = n_finished.value
            sys.stdout.write("\r(" + str(done)+"/"+str(self.worker)+") done.")
            sys.stdout.flush()
            time.sleep(5)

        print(result.get())
        pool.close()
        pool.join()
        
        
    def check_files(self):
        for i in range(self.worker):
            self.finished[i] = os.path.isfile(self.workdir + '/Histograms_%i.root'%i)

        return all(self.finished)

    def merge(self):
        hist_files = [self.workdir + "/Histograms_%i.root"%i for i in range(self.worker)]
        if(manager.finished.count(True) != manager.worker):
            raise BaseException("Not all expected root files were produced! Only (%i/%i) exist!"%(manager.finished.count(True),manager.worker))
        if(os.path.isfile("Histograms.root")):
            if("y" in raw_input('Histograms.root already exists! Do you want to replace it? [Y/N]').lower()):
                subprocess.call(["hadd","-f", "Histograms.root"]+hist_files)
            else:
                return 0
        else:
            subprocess.call(["hadd","Histograms.root"]+hist_files)

if(__name__ == "__main__"):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--selection',choices=["top","W","Zbb","all"],default = 'W',help='W/Zbb/top/all selection(s)')
    parser.add_argument('-n', '--workers', default = 4, type = int, help='number of workers')
    parser.add_argument('-m', '--merge', action='store_true', help='Merge hist root files' )
    parser.add_argument('--submit',action='store_true', help='submit workers')
    parser.add_argument('--macro',default='fillhists',help='specify which (compiled) macro to execute. (this also affects the name of workdir)')
    args = parser.parse_args()

    greeting = "== Fill Hists MT =="
    print("="*len(greeting))
    print(greeting)
    print("="*len(greeting))

    print("Selection:",args.selection)
    print("N Workers:",args.workers)
    
    manager = JobManager(args,"workdir_atlas_comparison" if 'atlas_comparison' in args.macro else "workdir")
    if(args.submit):
        manager.start()
    if(manager.check_files()):
        print("Success!")
    else:
        print("Only (%i/%i) root files were produced!"%(manager.finished.count(True),manager.worker))
    if(args.merge):
        manager.merge()
        

    
    
