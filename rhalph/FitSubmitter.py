from __future__ import print_function
import json,os,copy,subprocess
from CombineWorkflows import CombineWorkflows

class FitSubmitter:
    def __init__(self, config, workdir,dry_run = False):

        self._email = "steffen.albrecht@desy.de"

        self._workdir = os.path.realpath(workdir)
        
        self.rhalphdir = os.getcwd()
        self.CombinePath = self.rhalphdir + "/CMSSW_10_2_13"
            
        self.setup_workdir(self._workdir)

        self._QCDTFOrders = (None,None)
        self._DATATFOrders = (None,None)
        self.config = config
        self._dry_run = dry_run
        self.batchname = ""
        self.extra_jetmass_options = ""
        self.fit_qcd_model = False
        self.postfit_command = "" #will be written into htc wrapper which recieves fit config as parameter $1


    def setup_workdir(self,workdir):
        if(not os.path.isdir(workdir)):
            os.makedirs(workdir)
            scriptdir = os.getcwd()
            os.chdir(workdir)
            os.system('ln -sf '+self.CombinePath)
            os.system('ln -sf '+self.rhalphdir+'/setup_ralph.sh')
            os.system('ln -sf '+self.rhalphdir+'/runFit.py')
            os.system('ln -sf '+self.rhalphdir+'/rhalphalib')
            os.system('ln -sf '+self.rhalphdir+'/jetmass.py')
            os.system('ln -sf '+self.rhalphdir+'/fitplotter')
            os.system('ln -sf '+self.rhalphdir+'/CombineWorkflows.py')
            os.system('ln -sf '+self.rhalphdir+'/CombinePlotter.py')
            os.chdir(scriptdir)
        
    
    @property
    def config(self):
        return self._config
    
    @config.setter
    def config(self,c):
        if(isinstance(c,dict)):
            self._config = c
        else:
            try:
                self._config = json.load(open(c))
            except:
                raise(TypeError("Could not convert fit configuration. Valid input types are a path of json-file or python dict."))
        self._config['gridHistFileName'] = os.path.abspath(self._config['gridHistFileName'])
        self._config['histLocation'] = os.path.abspath(self._config['histLocation'])
    
    def do_fit_diagnostics_plots(self):
        fitdiagnostics_plotting = "python CombinePlotter.py --method plot_fit_diagnostics_results --parameter $1"
        print('adding to postfit_command: ',fitdiagnostics_plotting)
        self.postfit_command += fitdiagnostics_plotting

        
    def display(self):
        print('-> ModelBaseName:\t',self._config['ModelName'])
        print('-> histLocation:\t',self._config['histLocation'])
        print('-> gridHistFileName:\t',self._config['gridHistFileName'])
        print('-> workdir:\t\t',self._workdir)


    def submit_fits(self, job_workdirs, batchname = "JetMassJobs", defaultFit = True):
        
        # config_model_names = [config['ModelName'] for config in configs]
        #write HTCondor submissio jdl
        htc_jdl_name = self._workdir+'/'+batchname+'CondorSubmission.jdl'
        with open(htc_jdl_name, 'w') as condor_submit_jdl:
            condor_submit_jdl.write("""Requirements = ( OpSysAndVer == "CentOS7" )
universe          = vanilla
# #Running in local mode with 8 cpu slots
# universe          =  local
# request_cpus      =  8
notification      = Error
notify_user       = """ +self._email+"""
initialdir        = """+self._workdir+"""
output            = $(ConfigName).o$(ClusterId).$(Process)
error             = $(ConfigName).e$(ClusterId).$(Process)
log               = $(ConfigName).$(Cluster).log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
RequestMemory     = 4G
RequestDisk       = 4G
#getenv            = True
JobBatchName      = """+batchname+"""
executable        = """+self._workdir+"""/"""+batchname+"""_wrapper.sh
arguments         = $(ConfigName) $(FitJobIndex)
queue ConfigName from ("""+'\n'.join(job_workdirs)+"""
)""")
            condor_submit_jdl.close()

        batch_wrapper_name = self._workdir+'/'+batchname+'_wrapper.sh'
        batch_wrapper = open(batch_wrapper_name,'w')
        batch_wrapper.write("#!/bin/bash\n")
        batch_wrapper.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
        batch_wrapper.write("source setup_ralph.sh\n")
        jetmass_options = self.extra_jetmass_options +("" if defaultFit else " --customCombineWrapper")
        batch_wrapper.write('cd $1\n')

        # produce ModelDir with datacards using rhalphalib
        batch_wrapper.write('for i in *.json; do python jetmass.py $i %s;done\n'%jetmass_options) 
        # write wrapper that runs combine workflow if previous command did not do that yet!
        if('--build' in self.extra_jetmass_options):
            batch_wrapper.write('env -i bash $1/'+('qcdmodel/' if self.fit_qcd_model else '')+'wrapper.sh\n')
        batch_wrapper.write(self.postfit_command)
                        
        batch_wrapper.close()
        os.system('chmod +x '+batch_wrapper_name)
        command = ['condor_submit '+htc_jdl_name]
        print(' '.join(command))
        if(not self._dry_run):
            proc_condor_sub = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE)
            print(proc_condor_sub.communicate()[0])

        
    def scan_TF_orders(self,DATAOrders,QCDOrders=(None,None),combine_method=None):
        do_data_scan , do_qcd_scan = True, True
        if(DATAOrders is None):
            DATAOrders = (0,0)
        if(isinstance(DATAOrders[0],int)):
            do_data_scan = False
            DATAOrders = ([DATAOrders[0]],[DATAOrders[1]])
        if(QCDOrders is None):
            do_qcd_scan = False
            QCDOrders = (None,None)
        elif(isinstance(QCDOrders[0],int)):
            do_qcd_scan = False
            QCDOrders = ([QCDOrders[0]],[QCDOrders[1]])
           
        initial_qcd_fit = not (QCDOrders[0] is  None and QCDOrders[1] is None)
        if(not initial_qcd_fit):
            if(self.config.get("InitialQCDFitOrders","False")=="True"):
                QCDOrders = (self.config["InitialQCDFitOrders"][0],self.config["InitialQCDFitOrders"][1])
            else:
                QCDOrders = ([0],[0])
        jobs = []
        for QCD_pt in QCDOrders[0]:
            for QCD_rho in QCDOrders[1]:
                for DATA_pt in DATAOrders[0]:
                    for DATA_rho in DATAOrders[1]:
                        FTest = isinstance(combine_method,CombineWorkflows) and combine_method.method == "FTest" 
                        #print(QCD_pt,QCD_rho,DATA_pt,DATA_rho)
                        QCDOrders_str = 'MCTFPt%iRho%i'
                        DataOrders_str = 'DataResTFPt%iRho%i'
                        Orders_str = (QCDOrders_str%(QCD_pt,QCD_rho)+DataOrders_str%(DATA_pt,DATA_rho)) if initial_qcd_fit else 'TFPt%iRho%i'%(DATA_pt,DATA_rho)
                        
                        this_config = copy.deepcopy(self.config)
                        this_config['ModelName']+= Orders_str
                        this_config['InitialQCDFitOrders'] = [QCD_pt,QCD_rho]
                        this_config['InitialQCDFit'] = "True" if initial_qcd_fit else "False"
                        this_config['BernsteinOrders'] = [DATA_pt,DATA_rho]

                        job_workdir = self._workdir +"/"+this_config["ModelName"]+"/"
                        self.setup_workdir(job_workdir)
                        jobs.append(this_config["ModelName"])
                        
                        json.dump(this_config, open(job_workdir + '/' + this_config['ModelName']+'.json','w'))

                        if(FTest):
                            incrPt_config = copy.deepcopy(self.config)
                            incrPt_config['ModelName']+= (QCDOrders_str%(QCD_pt+(1 if do_qcd_scan else 0),QCD_rho)+DataOrders_str%(DATA_pt+(1 if do_data_scan else 0),DATA_rho)) if initial_qcd_fit else 'TFPt%iRho%i'%(DATA_pt+1,DATA_rho)
                            incrPt_config['InitialQCDFitOrders'] = [QCD_pt+(1 if do_qcd_scan else 0),QCD_rho]
                            incrPt_config['InitialQCDFit'] = "True" if initial_qcd_fit else "False"
                            incrPt_config['BernsteinOrders'] = [DATA_pt+(1 if do_data_scan else 0),DATA_rho]
                            json.dump(incrPt_config, open(job_workdir + '/' + incrPt_config['ModelName']+'.json','w'))

                            incrRho_config = copy.deepcopy(self.config)
                            incrRho_config['ModelName']+= (QCDOrders_str%(QCD_pt,QCD_rho+(1 if do_qcd_scan else 0))+DataOrders_str%(DATA_pt,DATA_rho+(1 if do_data_scan else 0))) if initial_qcd_fit else 'TFPt%iRho%i'%(DATA_pt,DATA_rho+1)
                            incrRho_config['InitialQCDFitOrders'] = [QCD_pt,QCD_rho+(1 if do_qcd_scan else 0)]
                            incrRho_config['InitialQCDFit'] = "True" if initial_qcd_fit else "False"
                            incrRho_config['BernsteinOrders'] = [DATA_pt,DATA_rho+(1 if do_data_scan else 0)]
                            json.dump(incrRho_config, open(job_workdir + '/' + incrRho_config['ModelName']+'.json','w'))

                        if(combine_method is not None):
                            if(isinstance(combine_method,str)):
                                cw = CombineWorkflows()
                                cw.method(combine_method)
                                cw.workspace = job_self._workdir + '/' + this_config['ModelName'] + '/'+this_config['ModelName'] + '_combined.root'
                                from jetmass import build_mass_scale_variations
                                cw.POI = build_mass_scale_variations(this_config['gridHistFileName'])[1]
                                cw.freezeParameters = "r"
                                cw.algo = "saturated"
                                # cw.extraOptions = "--preFitValue 0" 
                                # cw.bernsteinOrders = 'MCTFPt%iRho%iDataResTFPt%iRho%i'%(QCD_pt,QCD_rho,DATA_pt,DATA_rho) if initial_qcd_scan else 'TFPt%iRho%i'%(DATA_pt,DATA_rho)
                                cw.write_wrapper()
                            elif(isinstance(combine_method,CombineWorkflows)):
                                if(self.fit_qcd_model):
                                    # if(FTest):
                                    #     combine_method.workspace = 
                                    combine_method.workspace = job_workdir + '/' + this_config['ModelName'] + '/qcdmodel/qcdmodel_combined.root'
                                else:
                                    combine_method.workspace = job_workdir + '/' + this_config['ModelName'] + '/'+this_config['ModelName'] + '_combined.root'
                                combine_method.write_wrapper()
                            else:
                                raise(TypeError("combine_method should be either str or CombineWorkflow instance. You provided "+type(combine_method)))
                        
        # self.submit_fits(configs,batchname = "2TFScan" if initial_qcd_scan else "1TFScan", fitDiagnostics = (combine_method is None))
        self.submit_fits(jobs,batchname = self.batchname, defaultFit = (combine_method is None))
        
    def scan(self, attributes, scanName, combine_method=None):
        configs=[]
        for attribute, values in attributes.items():
            print(attribute)
            for value in values:
                this_config = copy.deepcopy(self.config)
                this_config['ModelName'] += value[0]
                this_config[attribute] = value[1]
                json.dump(this_config, open(self._workdir + '/' + this_config['ModelName']+'.json','w'))
                configs.append(this_config)
                if(combine_method is not None):
                    if(isinstance(combine_method,str)):
                        cw = CombineWorkflows()
                        cw.method(combine_method)
                        from jetmass import build_mass_scale_variations
                        cw.workspace = self._workdir + '/' + this_config['ModelName'] + '/'+this_config['ModelName'] + '_combined.root'
                        cw.POI = build_mass_scale_variations(this_config['gridHistFileName'])[1]
                        cw.freezeParameters = "r"
                        cw.algo = "saturated"
                        # cw.extraOptions = "--preFitValue 0"
                        cw.lumi = float(this_config.get('Pseudo',['lumiscale:1.0'])[0].split(':')[1]) * 41.8
                        cw.write_wrapper()
                    elif(isinstance(combine_method,CombineWorkflows)):
                        if(self.fit_qcd_model):
                            combine_method.workspace = self._workdir + '/' + this_config['ModelName'] + '/qcdmodel/qcdmodel_combined.root'
                        else:
                            combine_method.workspace = self._workdir + '/' + this_config['ModelName'] + '/'+this_config['ModelName'] + '_combined.root'
                        combine_method.write_wrapper()
                    else:
                        raise(TypeError("combine_method should be either str or CombineWorkflow instance. You provided "+type(combine_method)))


        self.submit_fits(configs,batchname = scanName, fitDiagnostics = (combine_method is None))

if(__name__ == '__main__'):
    # submitter = FitSubmitter("WJets.json","DUST/test_dir",dry_run=False)
    submitter = FitSubmitter("TTbarTopOnlyInclPt.json","DUST/TTbarTopOnlyPtIncl",dry_run=False)
    # submitter = FitSubmitter("TTbarWOnlyInclPt.json","DUST/TTbarWOnlyPtIncl",dry_run=False)
    submitter.display()
    # print('1TF:')
    # submitter.scan_TF_orders( (range(1,3),range(1,3)) , combine_method = "GOF")

    LumiScan = {'Pseudo':[
        ("L4fb",['lumiScale:0.1']),
        ("L41fb",['lumiScale:1.0']),
        ("L100fb",['lumiScale:2.3889154324']),
        ("L200fb",['lumiScale:4.7778308648'])]}

    # submitter.scan(LumiScan,"LumiScan",combine_method = "GOF")
    submitter.scan(LumiScan,"LumiScan")
    # print('2TF:')
    # submitter.scan_TF_orders( (range(1,3),range(1,3)) , (range(1,3),range(1,3)) )
