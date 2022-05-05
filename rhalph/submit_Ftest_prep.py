from jetmass import build_mass_scale_variations
import FitSubmitter
from CombineWorkflows import CombineWorkflows
import json

#dirty hack to keep this primarily python3 compatible, but still be able to run in py2 wihtout changes
try:
    range = xrange
except:
    pass
    
def submit_ftest_prep(config, ntoys, njobs, QCDOrders=None, DataOrders=None, algo = "saturated", dry_run=False,PseudoData=False):
    for i in range(njobs):
        do_qcd_test = False
        do_2TF = QCDOrders is not None
        
        if((QCDOrders is not None and (isinstance(QCDOrders[0],range) or isinstance(QCDOrders[0],list))) or DataOrders is None):
            do_qcd_test = True
        do_data_test = False
        if(not do_qcd_test and DataOrders is not None and  (isinstance(DataOrders[0],range) or isinstance(DataOrders[0],list))):
            do_data_test = True
        if(do_qcd_test and do_data_test):
            raise ValueError("You provided fixed orders for both QCD and Data TF. You have to provide a lists of orders that should be scanned for either TF!")
        if(not do_qcd_test and not do_data_test):
            raise ValueError("You provided running orders for both QCD and Data TF. You have to provide lists of orders that should be scanned for just one TF!")
        
        if(PseudoData or do_qcd_test):
            config["Pseudo"]=[""]
        else:
            config["Pseudo"]=[]

        #config["QCDSigmaScale"] = 0.01
            
        cw = CombineWorkflows()
        cw.combineCMSSW = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/CombineTestground/CMSSW_10_2_13/'
        cw.method = "FTestBatch"
        cw.seed = 42 + i
        cw.POI = "r"
        cw.freezeParameters = "r"
        
        cw.extraOptions = ("--toysFrequentist --setParameters r=1 " if do_data_test else  "") + ("--toysNoSystematics --setParameters r=0 " if do_qcd_test else "") +  " --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01" + " --redefineSignalPOIs massScale_pt0_eta0_all_Z_PT800,massScale_pt0_eta0_all_W_PT500,massScale_pt0_eta0_all_W_PT650,massScale_pt0_eta0_all_Z_PT650,massScale_pt0_eta0_all_Z_PT500,massScale_pt0_eta0_all_W_PT800"

        cw.algo = algo
        # cw.toysOptions = ("--toysFrequentist --expectSignal 1" if do_data_test else  "") + ("--toysNoSystematics --expectSignal 0" if do_qcd_test else "")
        cw.toysOptions = ""
        cw.toys = int(ntoys/njobs)

        workdir = "FTest_" + ("QCDTFScan" if do_qcd_test else "") + ("DataTFScan"+ ("_QCDOrder%ix%i"%QCDOrders if do_2TF else "") if do_data_test else "") + "_%s_Seed%i"%(cw.algo,cw.seed)

        fs = FitSubmitter.FitSubmitter(config, "DUST/"+workdir, dry_run = dry_run)
        fs.batchname = workdir.replace("_","")
        fs.CombinePath = cw.combineCMSSW
        # fs.extra_jetmass_options = "--build --noMassScales"
        fs.extra_jetmass_options = "--build"
        fs.fit_qcd_model = do_qcd_test
        fs.scan_TF_orders(DataOrders,QCDOrders,combine_method=cw)
 
def ftest_QCDMC_TF_fitDiagnostics_closure(config,QCDOrders = (range(0,7),range(0,7))):
    cw = CombineWorkflows()
    cw.set_method("diagnostics")
    cw.POI = "r"
    cw.freezeParameters = ""
    workdir = "FitDiagClosure_1TFScan"

    fs = FitSubmitter.FitSubmitter(config,"DUST/"+workdir,dry_run=dry_run)
    fs.batchname = workdir.replace("_","")
    fs.extra_jetmass_options = "--build" #just build workspace and not perform fit with main workspace
    fs.fit_qcd_model = True

    fs.do_fit_diagnostics_plots()
    
    fs.scan_TF_orders(([0],[0]),QCDOrders,combine_method = cw)


def ftest_Data_TF_fitDiagnostics_closure(config,Orders = (range(0,7),range(0,7)), QCDOrders=([3],[1]),PseudoData=False):
    if(PseudoData):
        config["Pseudo"]=[""]
    else:
        config["Pseudo"]=[]

    cw = CombineWorkflows()
    cw.set_method("diagnostics")
    cw.POI = "r"
    cw.freezeParameters = ""
    workdir = "FitDiagClosure_2TFScan_%sData_QCDOrders_%ix%i"%("Pseudo" if PseudoData else "",QCDOrders[0][0],QCDOrders[1][0])

    fs = FitSubmitter.FitSubmitter(config,workdir,dry_run=dry_run)
    fs.batchname = workdir.replace("_","")
    fs.extra_jetmass_options = "--build --noMassScales" #just build workspace and not perform fit with main workspace

    fs.do_fit_diagnostics_plots()
    
    fs.scan_TF_orders(Orders,QCDOrders,combine_method = cw)

if(__name__ == '__main__'):
    NToys = 500
    NToysPerJob=50
    NJobs=int(NToys/NToysPerJob)
    

    config = json.load(open("VJets.json"))
    # config = json.load(open("Zbb.json"))
    Orders = (range(0,7)),(range(0,7))
    dry_run = False
    for algo in ["saturated"]:#,"KS"]:#,"AD"]:
    # for algo in ["saturated"]:
        if(algo == "saturated"):
            # submit_ftest_prep(config, NToys, NJobs, QCDOrders=(0,4), DataOrders = Orders, algo=algo, dry_run=dry_run)
            submit_ftest_prep(config, NToys, NJobs, QCDOrders=None, DataOrders = Orders, algo=algo, dry_run=dry_run)
        # submit_ftest_prep(config, NToys, NJobs, QCDOrders=Orders,  algo=algo, dry_run=dry_run)
