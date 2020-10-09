from __future__ import print_function
import os, subprocess, copy, json

# scaleStudy_suffix="_0p01QCDMisstag"
scaleStudy_suffix=""

# selection = "WJets"
# selection = "WfromTop"
# selection = "WJets"
selection = "Top"

json_name_template = selection+"_scaleStudy_%s.json"
config_template_name = selection+'_scaleStudy.json'
rhalphdir = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/rhalph/'

def submit_variations(grids,config_template,suffix='',submit=True):
    workdir = 'workdir_scaleStudy_'+selection+suffix
    
    if(not os.path.isdir(workdir)):
        os.makedirs(workdir)
        os.chdir(workdir)
        os.system('ln -sf '+rhalphdir+'/CMSSW_10_2_13')
        os.system('ln -sf '+rhalphdir+'/setup_ralph.sh')
        os.system('ln -sf '+rhalphdir+'/runFit.py')
        os.system('ln -sf '+rhalphdir+'/rhalphalib')
        os.system('ln -sf '+rhalphdir+'/jetmass.py')
        os.system('ln -sf '+rhalphdir+'/fitplotter')
        os.chdir(rhalphdir+'/scaleStudy')

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
                        config['histLocation'] = rhalphdir+'/../macros/scaleStudy%s/Histograms_'%scaleStudy_suffix+name+'.root'
                        print(name,'-> vary pseudo like:','%i_%i_%s_'%(ptbin,etabin,pfflavour),direction)
                        config['Pseudo'].append(category)
                        config['Pseudo'].append(direction)
                        os.system('pwd')
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
    if(submit):
        print(' '.join(['condor_submit','CondorSubmission.jdl','-a "queue VaryPseudoLike in ('+' '.join(variation_names)+')"']))
        proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(variation_names)+')"'],shell=True, stdout=subprocess.PIPE)
        print(proc_condor_sub.communicate()[0])


def submit_categories(grids,config_template,suffix='',submit=True):
    workdir = 'workdir_nominal_'+selection+suffix
    # rhalphdir = os.getcwd()
    
    if(not os.path.isdir(workdir)):
        os.makedirs(workdir)
        os.chdir(workdir)
        os.system('ln -sf '+rhalphdir+'/CMSSW_10_2_13')
        os.system('ln -sf '+rhalphdir+'/setup_ralph.sh')
        os.system('ln -sf '+rhalphdir+'/runFit.py')
        os.system('ln -sf '+rhalphdir+'/rhalphalib')
        os.system('ln -sf '+rhalphdir+'/jetmass.py')
        os.system('ln -sf '+rhalphdir+'/fitplotter')
        os.chdir(rhalphdir+'/scaleStudy')

    names = []
    for name, categories in grids.items():
        print(name)
        config = copy.deepcopy(config_template)
        config['ModelName'] = name.replace('_','')
        config['gridHistFileName'] = rhalphdir+'/../Histograms/grid_'+name+'.root'
        config['histLocation'] = rhalphdir+'/../macros/scaleStudy%s/Histograms_'%scaleStudy_suffix+name+'.root'
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
    if(submit):
        print(' '.join(['condor_submit','CondorSubmission.jdl','-a "queue VaryPseudoLike in ('+' '.join(names)+')"']))
        proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(names)+')"'],shell=True, stdout=subprocess.PIPE)
        print(proc_condor_sub.communicate()[0])


if(__name__ == "__main__"):
    grids={
        # "oneScale":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["all"]
        # },
        "PF_flavours":{
            "ptbins":[0.,100000.],
            "etabins":[0.,9.],
            "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        }#,
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


    # pt_bins = ['_inclusive','_300To500','_300To400','_400To500','_500ToInf']
    # cycles = ["_300To500","_300To350","_350To400","_400To450","_450To500","_500To550","_550To600","_600ToInf"]
    cycles = []

    condor=False
    if len(cycles)>0:
        for pt_bin in cycles:
            config_template = json.load(open(config_template_name.replace('.json','_PtBin.json')))
            old_channel_name = list(config_template['channels'].keys())[0]
            config_template['channels'][old_channel_name+pt_bin.replace('_','').split('To')[0]] = config_template['channels'].pop(old_channel_name)
            
            config_template['channels'][old_channel_name+pt_bin.replace('_','').split('To')[0]]['pt_bin'] = pt_bin.replace('To','to').replace('_','')
            
            
            submit_variations(grids,config_template,pt_bin,submit=condor)
            submit_categories(grids,config_template,pt_bin,submit=condor)
        if(not condor):
            os.system('./submit_in_tmux.sh')
    else:
        config_template = json.load(open(config_template_name))
        submit_variations(grids,config_template,submit=False)
        submit_categories(grids,config_template,submit=False)    
