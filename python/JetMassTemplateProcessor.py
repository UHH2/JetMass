#!/nfs/dust/cms/user/albrechs/miniconda3/envs/coffea/bin/python
import awkward as ak
import numpy as np
from coffea import processor
from coffea import hist as coffeahist
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea.analysis_tools import PackedSelection
import coffea.lookup_tools#.dense_lookup import dense_lookup
#from coffea.lookup_tools.extractor import extractor
import hist

ddtmaps_path = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/Histograms/ddtmaps.npy'
kfactor_path = '/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/NLOweights/'

class JMSTemplates(processor.ProcessorABC):
    def __init__(self,year="2017", selection = 'vjets'):
        self._year = year
        self._selection = selection

        # dataset_ax = hist.Cat("dataset", "Dataset")
        # pT_ax  = hist.Bin("pt" , "$p_{T}$ [GeV]", 100, 0., 1500.)

        #bins: int, start: float, stop: float, *, name: str = '', label: str = '', 
        dataset_ax = hist.axis.StrCategory([], name="dataset", growth=True)
        shift_ax = hist.axis.StrCategory([], name="shift", growth=True)
        jec_applied_ax = hist.axis.StrCategory([], name='jecAppliedOn',label='JEC applied on',growth=True)
        mJ_ax  = hist.axis.Regular(50, 0., 500.,name="mJ" , label=r"$m_{SD}$ [GeV]")
        pT_ax  = hist.axis.Regular(300, 0., 3000., name="pt" , label=r"$p_{T}$ [GeV]")
        mJ_fit_ax  = hist.axis.Regular(250, 50., 300.,name="mJ" , label=r"$m_{SD}$ [GeV]")
        rho_ax  = hist.axis.Regular(100, -10., 0,name="rho" , label=r"$\rho$")

        self._pT_fit_ax = {
            'vjets':hist.axis.Variable(np.array([500, 650, 800, 1200,np.inf]), name="pt" , label=r"$p_{T}$ [GeV]"),
            'ttbar':hist.axis.Variable(np.array([200, 300, 400, 500, 650, np.inf]), name="pt" , label=r"$p_{T}$ [GeV]"),
        }

        # mJ_ax  = hist.Bin("mJ" , "$m_{SD}$ [GeV]", 50, 0., 500.)
        # dataset_ax = hist.Cat("dataset", "Dataset")
        # pT_ax  = hist.Bin("pt" , "$p_{T}$ [GeV]", 300, 0., 3000.)
        # mJ_fit_ax  = hist.Bin("mJ" , "$m_{SD}$ [GeV]", 250, 50., 300.)
        # rho_ax  = hist.Bin("rho" , "$\\rho$", 100, -10., 0.)

        # self._pT_fit_ax = {
        #     'vjets':hist.Bin("pt" , "$p_{T}$ [GeV]", np.array([500, 650, 800, 1200])),
        #     'ttbar':hist.Bin("pt" , "$p_{T}$ [GeV]", np.array([200, 300, 400, 500, 650, 100000])),
        # }
        
        hists = {}

        #create dense_lookup from custom n2ddt map
        n2ddtmap = np.load(ddtmaps_path,allow_pickle=True).item()
        corrected_str={
            'none':'',
            'pt':'_corrected_pt',
            'pt&mJ':'_corrected_pt_mass'

        }
        self._n2ddtmaps = {
            jec_applied_on:coffea.lookup_tools.dense_lookup.dense_lookup(
                n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}'][0],
                dims = (n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}'][1],
                        n2ddtmap[f'n2ddt_map_{year}_smooth_4_0p05{corrected_str[jec_applied_on]}'][2]) )
            for jec_applied_on in ['none','pt','pt&mJ']
        }
        #get some corrections and pack them into dense_lookups
        corrections_extractor = coffea.lookup_tools.extractor()
        for boson in ['W','Z']:
            fname = f'{kfactor_path}/{boson}JetsCorr.root'
            corrections_extractor.import_file(fname)
            corrections_extractor.add_weight_sets([
                f'{boson}_kfactor kfactor {fname}',
                f'{boson}_ewcorr ewcorr {fname}',
            ])
        corrections_extractor.finalize()
        
        self.corrections = corrections_extractor.make_evaluator()


        self._selections = ['vjets','ttbar']

        #these are selections that are linked as one specifies later in self._regions!
        self._cuts = {
            'vjets':['n2ddt','rhocut'],
            'ttbar':['tau32','tau21']
        }

        #these are sample specific matching criteria
        #these are handled via 'any' of PackedSelection, so in case of multiple requirements,
        #be aware, that they are linked with OR!!
        self._matching_mappings = {
            'WJetsMatched':{'IsMergedWZ':1},
            'WJetsUnmatched':{'IsMergedWZ':0},
            'ZJetsMatched':{'IsMergedWZ':1},
            'ZJetsUnmatched':{'IsMergedWZ':0},
            "TTbar_semilep_mergedTop":{'IsMergedTop':1},
            "TTbar_semilep_mergedW":{'IsMergedWZ':1},
            "TTbar_semilep_mergedQB":{'IsMergedQB':1},
            "TTbar_semilep_semiMergedTop":{'IsMergedWZ':1,'IsMergedQB':1},
            "TTbar_semilep_notMerged":{'IsNotMerged':1}
        }

        #common control-plots
        hists.update({
            
            'pt':hist.Hist(pT_ax, dataset_ax, jec_applied_ax,storage=hist.storage.Weight()),
            # 'pt_raw':hist.Hist(dataset_ax, pT_ax,storage=hist.storage.Weight()),
            'mjet':hist.Hist(mJ_ax, dataset_ax, jec_applied_ax,storage=hist.storage.Weight()),
            # 'mjet_raw':hist.Hist(dataset_ax, mJ_ax,storage=hist.storage.Weight()),
            'rho':hist.Hist(rho_ax, dataset_ax, jec_applied_ax,storage=hist.storage.Weight()),
            # 'rho_raw':hist.Hist(dataset_ax, rho_ax,storage=hist.storage.Weight()),
            'npv':hist.Hist(dataset_ax,hist.axis.Regular(80,0,80,name='npv',label=r'$N_{PV}$'),storage=hist.storage.Weight()),
            'ntrueint':hist.Hist(dataset_ax,hist.axis.Regular(80,0,80,name='ntrueint',label=r'$N_{TrueInt}$'),storage=hist.storage.Weight()),

        # 'pt':hist.Hist("Events",dataset_ax, pT_ax),
        #     'pt_raw':hist.Hist("Events",dataset_ax, pT_ax),
        #     'mjet':hist.Hist("Events",dataset_ax, mJ_ax),
        #     'mjet_raw':hist.Hist("Events",dataset_ax, mJ_ax),
        #     'rho':hist.Hist("Events",dataset_ax, rho_ax),
        #     'rho_raw':hist.Hist("Events",dataset_ax, rho_ax),
        #     'npv':hist.Hist("Events",dataset_ax,hist.axis.Regular(80,0,80,name='npv',label=r'$N_{PV}$')),
        #     'ntrueint':hist.Hist("Events",dataset_ax,hist.Bin(80,0,80,name='ntrueint',label=r'$N_{TrueInt}$')),
        })

        self._regions = {
            'vjets':{
                'inclusive':{'rhocut':True,'pt500cut':True},
                'pass':{'n2ddt':True,'rhocut':True,'pt500cut':True},
                'fail':{'n2ddt':False,'rhocut':True,'pt500cut':True},
            },
            'ttbar':{
                'inclusive':{'pt200cut':True},
                'pass' :{'tau32':True,'pt200cut':True},
                'passW':{'tau32':False,'tau21':True,'pt200cut':True},
                'fail' :{'tau32':False, 'tau21':False,'pt200cut':True},
                
            }
        }

        self._variations = ['0_0_all', '0_0_chargedH', '0_0_gamma', '0_0_neutralH', '0_0_other']
        for selection in self._selections:
            for region in self._regions[selection].keys():
                hists.update({
                    f"{selection}_mjet_{region}" : hist.Hist(mJ_ax,self._pT_fit_ax[selection],dataset_ax,jec_applied_ax,shift_ax,storage=hist.storage.Weight()),
                })
                for variation in self._variations:
                    hists.update({
                        f"{selection}_mjet_variation_{variation}_{region}__up" : hist.Hist(mJ_ax,self._pT_fit_ax[selection],dataset_ax,jec_applied_ax,shift_ax,storage=hist.storage.Weight()),
                    })
                    hists.update({
                        f"{selection}_mjet_variation_{variation}_{region}__down" : hist.Hist(mJ_ax,self._pT_fit_ax[selection],dataset_ax,jec_applied_ax,shift_ax,storage=hist.storage.Weight()),
                    })

        self._hists = lambda:hists#processor.dict_accumulator(hists)
    
    @property
    def accumulator(self):
        return self._hists

    def n2ddt(self,pt,rho,n2,corrected='none'):
        quantile =self._n2ddtmaps[corrected](rho,pt)
        return n2 - quantile

    def postprocess(self, accumulator):
        return accumulator
    
    def process(self, events):
        #print(self.accumulator)
        out = self.accumulator()

        
        dataset = events.metadata['dataset']

        isMC = 'data' not in dataset.lower()
        
        #evaluate matching criteria and created new masked events dataframe
        matching_mask = np.ones(len(events), dtype='bool')
        if(dataset in self._matching_mappings.keys()):
            matching_selection =  PackedSelection()
            for branch_name, branch_value in self._matching_mappings[dataset].items():
                matching_selection.add(branch_name, events[branch_name] == branch_value)
            #BE AWARE OF THE FOLLOWING OR!!!!
            matching_mask = matching_selection.any(*self._matching_mappings[dataset].keys())
        events = events[matching_mask]

        genpt = events.genjetpt
        if('WJets' in dataset):
            events.weight = events.weight*self.corrections['W_kfactor'](genpt)*self.corrections[f'W_ewcorr'](genpt)
        if('ZJets' in dataset):
            genpt = events.genjetpt
            events.weight = events.weight*self.corrections['Z_kfactor'](genpt)*self.corrections[f'Z_ewcorr'](genpt)

        
        nevents = len(events)
        _true = np.ones(nevents, dtype='bool')
        _false = np.zeros(nevents, dtype='bool')
        
        jecfactor = events.jecfactor
        
        pt = events.pt
        pt_raw = pt/jecfactor
        
        mjet_raw = events.mjet
        mjet = mjet_raw*jecfactor
        
        rho = 2*np.log(mjet/pt)
        rho_raw = 2*np.log(mjet_raw/pt_raw)
        rho_corrected_pt = 2*np.log(mjet_raw/pt)
        rho_corrected_mJ = 2*np.log(mjet/pt_raw)

        out['pt'].fill(dataset=dataset, jecAppliedOn='pt', pt=pt, weight = events.weight)
        out['pt'].fill(dataset=dataset, jecAppliedOn='none', pt=pt_raw, weight = events.weight)
        # out['pt_raw'].fill(dataset=dataset, pt=pt_raw)
        out['mjet'].fill(dataset=dataset, jecAppliedOn='mJ', mJ=mjet, weight = events.weight)
        out['mjet'].fill(dataset=dataset, jecAppliedOn='none', mJ=mjet_raw, weight = events.weight)
        # out['mjet_raw'].fill(dataset=dataset, mJ=mjet_raw)
        out['rho'].fill(dataset=dataset,jecAppliedOn='pt&mJ', rho=rho, weight = events.weight)
        out['rho'].fill(dataset=dataset,jecAppliedOn='pt', rho=rho_corrected_pt, weight = events.weight)
        out['rho'].fill(dataset=dataset,jecAppliedOn='mJ', rho=rho_corrected_mJ, weight = events.weight)
        out['rho'].fill(dataset=dataset,jecAppliedOn='none', rho=rho_raw, weight = events.weight)
        # out['rho_raw'].fill(dataset=dataset, rho=rho_raw)

        out['npv'].fill(dataset=dataset,npv = events.n_pv, weight = events.weight)
        if(isMC):
            out['ntrueint'].fill(dataset=dataset,ntrueint = events.n_trueint, weight = events.weight)

        
        ddb = events["MIDeepDoubleBHbbprob"]/(events["MIDeepDoubleBHbbprob"]+events["MIDeepDoubleBQCDprob"])
        for jec_applied_on in ['none','pt','pt&mJ']:
            selections = PackedSelection()
            pt_ = pt_raw
            if('pt' in jec_applied_on):
                pt_ = pt
            mJ_ = mjet_raw
            if('mJ' in jec_applied_on):
                mJ_ = mjet

            rho_ = 2*np.log(mJ_/pt_)
            
            selections.add("pt500cut",(pt_>500))
            selections.add("pt200cut",(pt_>200))

            selections.add("n2ddt",
                           (events.N2>0)      # make sure to fail on events with default values for topjets (do we need to remove those from fail region as well?)
                           & (self.n2ddt(pt_, rho_, events.N2, corrected=jec_applied_on)<0) # actual N2-DDT tagger
                           )
            selections.add("rhocut",
                           (rho_<-2.1)
                           &(rho_>-6.0))

            selections.add("tau21",events.tau21<0.45)
            selections.add("tau32",events.tau32<0.5)

        
            for selection in self._selections:
                for region in self._regions[selection].keys():
                    smask = selections.require(**self._regions[selection][region])
                    out[f"{selection}_mjet_{region}"].fill(
                        dataset=dataset,shift='nominal',
                        jecAppliedOn=jec_applied_on,
                        pt = pt_[smask],
                        mJ = mJ_[smask],
                        weight = events.weight[smask]
                    )
                    
                    if(isMC):
                        for variation in self._variations:
                            mJVar_ = events[f'mjet_{variation}']
                            if('mJ' in jec_applied_on):
                                mJVar_ = mJVar_ * jecfactor
                            out[f"{selection}_mjet_variation_{variation}_{region}__up"].fill(
                                dataset=dataset,shift='nominal',
                                jecAppliedOn=jec_applied_on,
                                pt = pt_[smask],
                                mJ = mJVar_[:,0][smask],
                                weight = events.weight[smask]
                            )
                            out[f"{selection}_mjet_variation_{variation}_{region}__down"].fill(
                                dataset=dataset,shift='nominal',
                                jecAppliedOn=jec_applied_on,
                                pt = pt_[smask],
                                mJ = mJVar_[:,1][smask],
                                weight = events.weight[smask]
                            )
            
        return out

    def postprocess(self, accumulator):
        return accumulator


    
if(__name__ == "__main__"):
    from coffea.util import load, save
    import os,glob
    from coffea_util import CoffeaWorkflow

    workflow  = CoffeaWorkflow("JMSTemplates")

    workflow.parser.add_argument("--output","-o",type=str,default='jms_templates.coffea')
    workflow.parser.add_argument("--year", default='2017')
    args = workflow.parse_args()

    workflow.processor_instance = JMSTemplates(args.year,args.selection)
    workflow.processor_schema = BaseSchema

    sample_pattern = '/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/{SELECTION}Trees/workdir_{SELECTION}_{YEAR}/*{SAMPLE}*.root'
    sample_names = ['Data','WJets','ZJets','QCD']
    samples = {sample:glob.glob(sample_pattern.format(SELECTION=args.selection,YEAR=args.year,SAMPLE=sample)) for sample in sample_names}

    if('WJets' in samples):
        samples['WJetsMatched'] = samples['WJets']
        samples['WJetsUnmatched'] = samples['WJets']
        del samples['WJets']
    if('ZJets' in samples):
        samples['ZJetsMatched'] = samples['ZJets']
        samples['ZJetsUnmatched'] = samples['ZJets']
        del samples['ZJets']

    
    for k,v in samples.items():
        print(k,len(v))

    output_file_path = os.path.join(os.getcwd(),args.output)
    print('changing into /tmp dir')
    os.chdir(os.environ['TMPDIR'])
    if(args.scaleout>0):
        print('init dask client')
        workflow.init_dask_htcondor_client(1,10,5)

    print('starting coffea runner')
    output = workflow.run(samples)

            
    if(not output_file_path.endswith('.coffea')):
        output_file_path += '.coffea'
    
    save(output, output_file_path)
