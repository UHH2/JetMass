from __future__ import print_function
import os, subprocess, copy, json


json_name_template = "WJets_scaleStudy_%s.json"


def submit_variations(grids,config_template):
    workdir = 'workdir_scaleStudy'
    rhalphdir = os.getcwd()
    
    if(not os.path.isdir(workdir)):
        os.makedirs(workdir)
        os.chdir(workdir)
        os.system('ln -sf '+rhalphdir+'/CMSSW_10_2_13')
        os.system('ln -sf '+rhalphdir+'/setup_ralph.sh')
        os.system('ln -sf '+rhalphdir+'/runFit.py')
        os.system('ln -sf '+rhalphdir+'/rhalphalib')
        os.system('ln -sf '+rhalphdir+'/jetmass.py')
        os.system('ln -sf '+rhalphdir+'/plotter')
        os.chdir(rhalphdir)

    variation_names = []
    for name, categories in grids.items():
        print(name)

        for ptbin in range(len(categories['ptbins'])-1):
            for etabin in range(len(categories['etabins'])-1):
                for pfflavour in categories['pfflavours']:
                    for direction in ['up','down']:
                        category = '%i_%i_%s_'%(ptbin,etabin,pfflavour)                        
                        variation_name = name+'_'+category+direction
                        config = copy.deepcopy(config_template)
                        config['ModelName'] = variation_name.replace('_','')
                        config['gridHistFileName'] = rhalphdir+'/../Histograms/grid_'+name+'.root'
                        config['histLocation'] = rhalphdir+'/../macros/scaleStudy/Histograms_'+name+'.root'
                        print(name,'-> vary pseudo like:','%i_%i_%s_'%(ptbin,etabin,pfflavour),direction)
                        config['Pseudo'].append(category)
                        config['Pseudo'].append(direction)
                        json.dump(config,open(workdir+'/'+json_name_template%(variation_name),'w'))
                        variation_names.append(variation_name)

    print(variation_names)

    with open(workdir+'/batch_wrapper.sh','w') as wrapper:
        wrapper.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source setup_ralph.sh
python jetmass.py $1
""")
        wrapper.close()
    os.system('chmod +x '+workdir+'/batch_wrapper.sh')
    
    with open(workdir+'/CondorSubmission.jdl','w') as condor_submit_yml:
        condor_submit_yml.write("""Requirements = ( OpSysAndVer == "CentOS7" )
universe          = vanilla
# #Running in local mode with 8 cpu slots
# universe          =  local
# request_cpus      =  8
notification      = Error
notify_user       = steffen.albrecht@desy.de
initialdir        = """+workdir+"""
output            = $(VaryPseudoLike).o$(ClusterId).$(Process)
error             = $(VaryPseudoLike).e$(ClusterId).$(Process)
log               = $(VaryPseudoLike).$(Cluster).log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
RequestMemory     = 4G
RequestDisk       = 4G
#getenv            = True
JobBatchName      = 
executable        = """+workdir+"""/batch_wrapper.sh
arguments         = """+json_name_template%'$(VaryPseudoLike)'+"""
""")
        condor_submit_yml.close()
    print(' '.join(['condor_submit','CondorSubmission.jdl','-a "queue VaryPseudoLike in ('+' '.join(variation_names)+')"']))
    proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(variation_names)+')"'],shell=True, stdout=subprocess.PIPE)
    print(proc_condor_sub.communicate()[0])


def submit_categories(grids,config_template):
    workdir = 'workdir_scaleStudy_vqq_normnuisance'
    rhalphdir = os.getcwd()
    
    if(not os.path.isdir(workdir)):
        os.makedirs(workdir)
        os.chdir(workdir)
        os.system('ln -sf '+rhalphdir+'/CMSSW_10_2_13')
        os.system('ln -sf '+rhalphdir+'/setup_ralph.sh')
        os.system('ln -sf '+rhalphdir+'/runFit.py')
        os.system('ln -sf '+rhalphdir+'/rhalphalib')
        os.system('ln -sf '+rhalphdir+'/jetmass.py')
        os.system('ln -sf '+rhalphdir+'/plotter')
        os.chdir(rhalphdir)

    names = []
    for name, categories in grids.items():
        print(name)
        config = copy.deepcopy(config_template)
        config['ModelName'] = name.replace('_','')
        config['gridHistFileName'] = rhalphdir+'/../Histograms/grid_'+name+'.root'
        config['histLocation'] = rhalphdir+'/../macros/scaleStudy/Histograms_'+name+'.root'
        json.dump(config,open(workdir+'/'+json_name_template%(name),'w'))
        names.append(name)


    with open(workdir+'/batch_wrapper.sh','w') as wrapper:
        wrapper.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source setup_ralph.sh
python jetmass.py $1
""")
        wrapper.close()
    os.system('chmod +x '+workdir+'/batch_wrapper.sh')
    
    with open(workdir+'/CondorSubmission.jdl','w') as condor_submit_yml:
        condor_submit_yml.write("""Requirements = ( OpSysAndVer == "CentOS7" )
universe          = vanilla
# #Running in local mode with 8 cpu slots
# universe          =  local
# request_cpus      =  8
notification      = Error
notify_user       = steffen.albrecht@desy.de
initialdir        = """+workdir+"""
output            = $(VaryPseudoLike).o$(ClusterId).$(Process)
error             = $(VaryPseudoLike).e$(ClusterId).$(Process)
log               = $(VaryPseudoLike).$(Cluster).log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
RequestMemory     = 4G
RequestDisk       = 4G
#getenv            = True
JobBatchName      = 
executable        = """+workdir+"""/batch_wrapper.sh
arguments         = """+json_name_template%'$(VaryPseudoLike)'+"""
""")
        condor_submit_yml.close()
    print(' '.join(['condor_submit','CondorSubmission.jdl','-a "queue VaryPseudoLike in ('+' '.join(names)+')"']))
    proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(names)+')"'],shell=True, stdout=subprocess.PIPE)
    print(proc_condor_sub.communicate()[0])


if(__name__ == "__main__"):
    grids={
        "oneScale":{
            "ptbins":[0.,100000.],
            "etabins":[0.,9.],
            "pfflavours":["all"]
        }#,
        # "PF_flavours":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # },
        # "2ptbins":{
        #     "ptbins":[0.,10.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["all"]
        # },
        # "2etabins":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0., 1.479, 9.],
        #     "pfflavours":["all"]
        # },
        # "PF_flavours_2ptbins":{
        #     "ptbins":[0.,10.,100000.],
        #     "etabins":[0., 9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # },
        # "PF_flavours_2etabins":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0., 1.479, 9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # }
    }
    config_template = json.load(open('WJets_scaleStudy.json'))
    
    submit_variations(grids,config_template)
    # submit_categories(grids,config_template)
