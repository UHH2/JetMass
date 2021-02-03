from jetmass import jet_mass_producer
import argparse,sys,json,copy,os,subprocess

# sys.path.append("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/util/python")
sys.path.append("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/util/python")
from combine_workflow import plot_gof_result

selection = 'WJets'
N_toys = 100
N_toys_per_fit = 100
N_jobs_per_fit = int(N_toys/N_toys_per_fit)

cycle_name =".MassScaleFrozen"
freezeParameters = "massScale"
suffix="NoInitialQCDFit"
# scaleCategory = "oneScale"
scaleCategory = "PF_flavours"

#suffix="InternToys"
def submit_bernstein_orders(configs,pt_orders, rho_orders,args):
    workdir = 'workdir_bernstein_gof_'+suffix+'_'+selection
    rhalphdir = os.getcwd()
    if(not os.path.isdir(workdir)):
        os.makedirs(workdir)
        os.chdir(workdir)
        os.system('ln -sf '+rhalphdir+'/CMSSW_10_2_13')
        os.system('ln -sf '+rhalphdir+'/setup_ralph.sh')
        os.system('ln -sf '+rhalphdir+'/runFit.py')
        os.system('ln -sf '+rhalphdir+'/rhalphalib')
        os.system('ln -sf '+rhalphdir+'/jetmass.py')
        os.system('ln -sf '+rhalphdir+'/fitplotter')
        os.chdir(rhalphdir)


    config_names = []
    for qcd_pt_order in [2,3,4,5,6]:
        for qcd_rho_order in [2,3,4,5,6]:
            this_configs = copy.deepcopy(configs)
            this_configs['ModelName']+='Pt%iRho%i'%(qcd_pt_order,qcd_rho_order)
            this_configs['BernsteinOrders']=[qcd_pt_order,qcd_rho_order]
            this_configs['gridHistFileName'] = os.path.abspath(rhalphdir+'/../Histograms/grid_'+scaleCategory+'.root')
            this_configs['histLocation'] = os.path.abspath(rhalphdir+'/../macros/scaleStudy/Histograms_'+scaleCategory+'.root')

            print(this_configs['ModelName'])
            print(this_configs['BernsteinOrders'])
            if(args.merge):
                # os.chdir(workdir+'/'+this_configs['ModelName'])
                os.chdir(workdir)
                hadd_command = 'hadd -f higgsCombine%s.GoodnessOfFit.mH0.merged.root '%(cycle_name+this_configs['ModelName']) + this_configs['ModelName']+'*/higgsCombine%s.GoodnessOfFit.mH0.*.root'%cycle_name
                # hadd_command = 'hadd  higgsCombine%s.GoodnessOfFit.mH0.merged.root '%cycle_name + ' '.join(['higgsCombine%s.GoodnessOfFit.mH0.merged.%i.root'%(cycle_name,i) for i in range(N_jobs_per_fit)])
                print(hadd_command)
                os.system(hadd_command)
                try:
                    plot_gof_result(this_configs['ModelName']+"0/higgsCombine%sBaseline.GoodnessOfFit.mH0.root"%cycle_name,"higgsCombine%s.GoodnessOfFit.mH0.merged.root"%(cycle_name+this_configs['ModelName']),freezeParameters+"Frozen_"+this_configs['ModelName'])
                except:
                    print('something went wrong. probably when hadding the trees')
                # os.chdir('../../')
                os.chdir('../')
                continue
            json.dump(this_configs,open(workdir+'/'+this_configs['ModelName']+'.json','w'))
            config_names.append(this_configs['ModelName'])
            
    if(args.merge):
        return

        
    util_path = os.path.abspath("../util/python")

    with open(workdir+'/build_wrapper.sh','w') as build_wrapper:
        build_wrapper.write("""#!/bin/bash
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        source setup_ralph.sh
        python jetmass.py $1.json --build --job_index $2 --noMassScales
        """)
        build_wrapper.close()
    os.system('chmod +x '+workdir+'/build_wrapper.sh')
    with open(workdir+'/batch_wrapper.sh','w') as wrapper:
        wrapper.write("""#!/bin/bash
        # source /cvmfs/cms.cern.ch/cmsset_default.sh
        # source setup_ralph.sh
        # python jetmass.py $1.json
        env -i bash build_wrapper.sh $1 $2
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd CMSSW_10_2_13
        eval `scramv1 runtime -sh`
        cd ../
        export PATH=$PATH:"""+util_path+"""
        cd $1$2
        source build.sh
        # combine_workflow.py --method GOF -w $1$2_combined.root --freezeParameters """+freezeParameters+""" --name """+cycle_name+""" -t """+str(N_toys_per_fit)+""" --job_index $2 --skipplots
        combine_workflow.py --method GOF -w $1$2_combined.root --name """+cycle_name+""" -t """+str(N_toys_per_fit)+""" --job_index $2 --skipplots 
        """)
        wrapper.close()
    os.system('chmod +x '+workdir+'/batch_wrapper.sh')
    
    
    argument_strings=[]
    for i in range(N_jobs_per_fit):
        argument_strings += [config_name + ' ' + str(i) for config_name in config_names]
        
    with open(workdir+'/CondorSubmission.jdl','w') as condor_submit_yml:
        condor_submit_yml.write("""Requirements = ( OpSysAndVer == "CentOS7" )
universe          = vanilla
# #Running in local mode with 8 cpu slots
# universe          =  local
# request_cpus      =  8
notification      = Error
notify_user       = steffen.albrecht@desy.de
initialdir        = """+workdir+"""
output            = $(ConfigName).o$(ClusterId).$(Process)
error             = $(ConfigName).e$(ClusterId).$(Process)
log               = $(ConfigName).$(Cluster).log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
RequestMemory     = 4G
RequestDisk       = 4G
#getenv            = True
JobBatchName      = QCDBernsteinGOF"""+suffix+"""
executable        = """+workdir+"""/batch_wrapper.sh
arguments         = $(ConfigName) $(FitJobIndex)
queue ConfigName,FitJobIndex from ("""+'\n'.join(argument_strings)+"""
)""")
        condor_submit_yml.close()
    command = ['condor_submit '+workdir+'/CondorSubmission.jdl']
    print(' '.join(command))
    proc_condor_sub = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE)
    print(proc_condor_sub.communicate()[0])



if(__name__ == "__main__"):
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with base config")
    parser.add_argument('--merge',action="store_true",help='merge results from toys and plot')
    args = parser.parse_args()
    
    try:
        configs = json.load(open(args.config))
    except IndexError:
        print("You must specify a configuration JSON!")
        sys.exit(0)

    submit_bernstein_orders(configs,range(2,7),range(2,7),args)
