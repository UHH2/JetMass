#!/usr/bin/python
from datetime import datetime
import sys, time, subprocess
from ROOT import *
from multiprocessing import Pool, Value
import random
import shutil
import glob
import os

Ndone = None # keep track of finished jobs by each worker

def main():
    if len(sys.argv) < 3:
        print 'Usage: ./SFrameMainMultiCore.py <Name Workdir> <Number of Jobs>'
        exit(0)

    workdir = sys.argv[1]
    nworkers = int(sys.argv[2])

    print '==============================='
    print '====== SFrame Multi Core ======'
    print '==============================='

    print '  -- Workdir: '+workdir
    print '  -- Number of Jobs:',nworkers

    if not os.path.isdir(workdir):
        print '[ERROR] workdir ('+workdir+') does not exist!'
        exit(0)

    if not os.path.exists(workdir+"/missing_files.txt"):
        print '[ERROR] not a proper sframe_batch workdir, missing_files.txt not found!'
        exit(0)

    # read missing_files.txt
    lines = read_missing_files(workdir)

    # remove everything but the name of the xml file
    xmls = []
    for line in lines:
        if "  sframe_main " in line:
            root, xml = line.split("  sframe_main ",1)
            xmls.append(workdir+'/'+xml)

    # get number of xmls and calculate xmls per job
    nxmls = len(xmls)
    print '  -- Number of xml files:',nxmls

    if nworkers > nxmls :
        print '[ERROR] you selected more workers than number of xmls'
        exit(0)

    # setup global counter
    Ndone = Value('i', 0)

    # start workers
    pool = Pool(processes = nworkers, initializer = init, initargs = (Ndone, ))
    # create lists for every worker including the xmls that should be processed
    # also pair it with index of worker (for keeping track of done jobs)
    job_lists = []
    for i in range(nworkers):
        job_lists.append( create_job_list(i, nworkers, xmls) )

    result = pool.map_async(submit_job,job_lists,chunksize=1)

    while not result.ready():
        with Ndone.get_lock():
            percent = Ndone.value*100/nxmls
            missing = nxmls - Ndone.value
        sys.stdout.write("\r" + str(missing) + " jobs left - " +loading_bar(percent))
        sys.stdout.flush()
        time.sleep(1)

    pool.close()
    pool.join()

    print ' '
    print '  -- sframe_main finished'
    print '  -- now: deleting temporary directories'
    delete_temp_dirs()
    print '  -- finished!'


# ------------------------------------------------------------------------------
# needed to setup global counter
def init(args):
    global Ndone
    Ndone = args
# ------------------------------------------------------------------------------
# this function submits the job and increments the global counter
# (global counter is not that easy with multiprocessing)
def submit_job(job_list):
    #print len(job_list)
    for xml in job_list:
        # subprocess.call(['sframe_main',xml],stdout=subprocess.PIPE)
        subprocess.call(['sframe_main',xml], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
        global Ndone
        with Ndone.get_lock():
            Ndone.value += 1

    return 1

# ------------------------------------------------------------------------------
# reads missing_files.txt and returns an array where every entry contains one line
def read_missing_files(workdir):
    with open(workdir+"/missing_files.txt", "r") as ins:
        array = []
        for line in ins:
            array.append(line.rstrip('\n'))
    return array

# ------------------------------------------------------------------------------
# for every worker an array containing different xmls is created
# it is taken care of distributing the job as equally as possible
def create_job_list(i, nworkers, xmls):
    integer_div = len(xmls) // nworkers
    rest = len(xmls) - (integer_div * nworkers)
    new_xml_list = []
    # the rest should not just be added to the last worker but distributed equally
    # therefore every worker with i+1 gets one additional job
    if i+1 <= rest:
        min_index = i*integer_div + i
        max_index = (i+1)*integer_div + 1 + i
    else:
        min_index = i*integer_div + rest
        max_index = (i+1)*integer_div + rest

    # loop over total list and only keep element if index is between min and max
    for index, val in enumerate(xmls):
        if index >= min_index and index < max_index:
            new_xml_list.append( xmls[index] )

    return new_xml_list
# ------------------------------------------------------------------------------
# this function deletes the temporary output of sframe_main
def delete_temp_dirs():
    path = os.getcwd()
    del_paths = glob.glob(os.path.join(path, 'jobTempOutput_*'))
    for del_path in del_paths:
        shutil.rmtree(del_path)
# ------------------------------------------------------------------------------
# this converts percent into a loading bar and returns a string containing the bar
def loading_bar(percent):
    p = float("{0:.2f}".format(percent)) # round to two digits
    bar = 'Progress: '+str(p)+'% ['
    levels = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95]
    for level in levels:
        if percent > level:
            bar += '=='
        else:
            bar += '  '
    bar += ']'
    return bar

# ------------------------------------------------------------------------------
# main function
if __name__ == "__main__":
    main()
