from __future__ import print_function
import os, subprocess, copy, json, argparse

rhalphdir = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/'

# scaleStudy_suffix="_0p01QCDMisstag"
scaleStudy_suffix="/QCDMisstag_0p05/"
# scaleStudy_suffix=""

# selection = "WJets"
# selection = "WfromTop"
# selection = "WJets"
# selection = "WJets_CHFLowBin"
# selection = "Top"
builder_arguments = ' --combineOptions "--freezeParameters r --preFitValue 0"'
# workdir_prefix = "../DUST/scaleStudy/"
workdir_prefix = ""

# json_name_template = selection+"_scaleStudy_%s.json"
# config_template_name = selection+'_scaleStudy.json'

def plot_scaleStudy_result(grids,workdir_prefix):
    import plot_nuisance_effect as scaleStudyPlotter
    scaleStudyPlotter.prefix = workdir_prefix
    for name,config in grids.items():
        scaleStudyPlotter.plot_nuisance_effect_grid(name, config, workdir_prefix=workdir_prefix)


def submit_scaleStudy(grids, config_template_name="WJets_scaleStudy.json", workdir_prefix="WJets_scaleStudy", submit = False):
    # config_template = json.load(open(config_template_name))
    submit_variations(grids,config_template_name,submit=submit,workdir_prefix=workdir_prefix)
    submit_nominal(grids,config_template_name,submit=submit,workdir_prefix=workdir_prefix)    


def submit_variations(grids,config_template_name,workdir_prefix='',submit=True):
    config_template = json.load(open(config_template_name))
    json_name_template = config_template_name.replace('.json','_%s.json')
    # workdir = workdir_prefix+'workdir_scaleStudy_'+selection+suffix
    workdir = workdir_prefix+'/workdir_variations'
    
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
python jetmass.py $1 """+builder_arguments+"""
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
    if(submit):
        proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(variation_names)+')"'],shell=True, stdout=subprocess.PIPE)
        print(proc_condor_sub.communicate()[0])

def submit_nominal(grids,config_template_name,workdir_prefix='',submit=True):
    config_template = json.load(open(config_template_name))
    json_name_template = config_template_name.replace('.json','_%s.json')
    # workdir = workdir_prefix+'workdir_nominal_'+selection+suffix
    workdir = workdir_prefix+'/workdir_nominal'
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
python jetmass.py $1 """+builder_arguments+"""
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
    if(submit):
        proc_condor_sub = subprocess.Popen(['condor_submit '+workdir+'/CondorSubmission.jdl -a "queue VaryPseudoLike in ('+' '.join(names)+')"'],shell=True, stdout=subprocess.PIPE)
        print(proc_condor_sub.communicate()[0])


def submit_bernstein_scan(grids, config_template, pt_orders, rho_orders, pt_orders_MC=[1], rho_orders_MC=[2], suffix='', submit=True):    

    for qcd_pt_order_MCTF in pt_orders_MC:
        for qcd_rho_order_MCTF in rho_orders_MC:
            for qcd_pt_order in pt_orders:
                for qcd_rho_order in rho_orders:
                    this_configs = copy.deepcopy(config_template)
                    TF_substring = ""
                    if(this_configs.get("InitialQCDFit","False") == "True"):
                        TF_substring = 'TFPt%iRho%iTFPt%iRho%i'%(qcd_pt_order_MCTF,qcd_rho_order_MCTF,qcd_pt_order,qcd_rho_order)
                    else:
                        TF_substring = 'Pt%iRho%i'%(qcd_pt_order,qcd_rho_order)
                        
                    this_configs['ModelName']+= TF_substring
                    this_configs['BernsteinOrders']=[qcd_pt_order,qcd_rho_order]
                    this_configs['InitialQCDFitOrders']=[qcd_pt_order_MCTF,qcd_rho_order_MCTF]
                    # this_configs['gridHistFileName'] = os.path.abspath(rhalphdir+'/../Histograms/grid_'+scaleCategory+'.root')
                    # this_configs['histLocation'] = os.path.abspath(rhalphdir+'/../macros/scaleStudy/QCDMisstag_0p05/Histograms_'+scaleCategory+'.root')
                    
                    submit_variations(grids,this_configs,suffix=suffix+TF_substring,submit=submit)
                    submit_nominal(grids,this_configs,suffix=suffix+TF_substring,submit=submit)
                    
if(__name__ == "__main__"):
    grids={
        "oneScale":{
            "ptbins":[0.,100000.],
            "etabins":[0.,9.],
            "pfflavours":["all"]
        },
        # "PF_flavours":{
        #     "ptbins":[0.,100000.],
        #     "etabins":[0.,9.],
        #     "pfflavours":["chargedH", "neutralH", "gamma", "other"]
        # }#,
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
    
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with base config", default="WJets_CHFHighBin_scaleStudy")
    parser.add_argument("--workdir",type=str,default="")
    parser.add_argument('--submit',action="store_true",help='submit scaleStudy workflow')
    parser.add_argument('--plot',action="store_true",help='plot scaleStudy workflow results')
    parser.add_argument('--debug',action="store_true")
    args = parser.parse_args()

    config_name = args.config.replace('.json','')

    if(args.submit):
        submit_scaleStudy(grids,config_template_name = config_name + '.json', workdir_prefix=args.workdir+config_name, submit = not args.debug)
    if(args.plot):
        plot_scaleStudy_result(grids,args.workdir+config_name)
    exit(0)
    # pt_bins = ['_inclusive','_300To500','_300To400','_400To500','_500ToInf']
    # cycles = ["_300To500","_300To350","_350To400","_400To450","_450To500","_500To550","_550To600","_600ToInf"]
    cycles = []

    bernstein_scan = False
    if(bernstein_scan):
        # submit_bernstein_scan(grids, json.load(open(config_template_name)), range(1,6), range(1,6), suffix='test', submit=True)
        submit_bernstein_scan(grids, json.load(open(config_template_name)), range(1,6), range(1,6), pt_orders_MC=range(1,6), rho_orders_MC=range(1,6), suffix='10GeVBins', submit=True)  
    
        exit(0)
    
    condor=False
    if len(cycles)>0:
        for pt_bin in cycles:
            config_template = json.load(open(config_template_name.replace('.json','_PtBin.json')))
            old_channel_name = list(config_template['channels'].keys())[0]
            config_template['channels'][old_channel_name+pt_bin.replace('_','').split('To')[0]] = config_template['channels'].pop(old_channel_name)
            
            config_template['channels'][old_channel_name+pt_bin.replace('_','').split('To')[0]]['pt_bin'] = pt_bin.replace('To','to').replace('_','')
            
            
            submit_variations(grids,config_template,pt_bin,submit=condor)
            submit_nominal(grids,config_template,pt_bin,submit=condor)
        if(not condor):
            os.system('./submit_in_tmux.sh')
    else:
        cycle = "WJets_CHFHighBin_scaleStudy"
        submit_scaleStudy(grids, config_template_name=cycle+".json", workdir_prefix=cycle, submit = True)
        
