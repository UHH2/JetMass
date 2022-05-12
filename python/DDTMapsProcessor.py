#!/nfs/dust/cms/user/albrechs/miniconda3/envs/coffea/bin/python
import json
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea import processor,hist
import mplhep as hep
hep.style.use('CMS')
import numpy as np

class DDTMapPrep(processor.ProcessorABC):
    def __init__(self):
        hists = {}
        
        year_axis = hist.Cat("year", "Year")

        pT_ax  = hist.Bin("pt" , "$p_{T}$ [GeV]", 100, 0., 1500.)
        rho_ax  = hist.Bin("rho" , "$\\rho$ [GeV]", 50, -10., 0.)
        N2_ax  = hist.Bin("n2" , "N2", 300, -1.5, 1.5)        
        
        hists.update({'leadingJet':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})
        hists.update({'leadingJet_corrected_pt_mass':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})
        hists.update({'leadingJet_corrected_pt':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})

        hists.update({'leadingJet_msd':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})
        hists.update({'leadingJet_msd_corrected_pt_mass':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})
        hists.update({'leadingJet_msd_corrected_pt':hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)})
        
        self._hists = processor.dict_accumulator(hists)

    @property
    def accumulator(self):
        return self._hists

    def process(self, events):
        out = self.accumulator.identity()
        
        year = events.metadata['dataset']
        
        jecfactor = events.jecfactor
        
        
        pt = events.pt
        pt_raw = pt/jecfactor
        
        mjet_raw = events.mjet
        mjet = mjet_raw*jecfactor
        
        rho = 2*np.log(mjet_raw/pt_raw)
        out['leadingJet'].fill(year=year, pt = pt_raw, rho = rho, n2 = events.N2, weight = events.weight)
        
        rho_corrected_pt = 2*np.log(mjet_raw/pt)
        out['leadingJet_corrected_pt'].fill(year=year, pt = pt, rho = rho_corrected_pt, n2 = events.N2, weight = events.weight)

        rho_corrected_pt_mass = 2*np.log(mjet/pt)
        out['leadingJet_corrected_pt_mass'].fill(year=year, pt = pt, rho = rho_corrected_pt_mass, n2 = events.N2, weight = events.weight)

        rho_SD = 2*np.log(mjet_SD_raw/pt_raw)
        out['leadingJet'].fill(year=year, pt = pt_raw, rho = rho_SD, n2 = events.N2, weight = events.weight)
        
        rho_SD_corrected_pt = 2*np.log(mjet_SD_raw/pt)
        out['leadingJet_corrected_pt'].fill(year=year, pt = pt, rho = rho_SD_corrected_pt, n2 = events.N2, weight = events.weight)

        rho_SD_corrected_pt_mass = 2*np.log(mjet_SD/pt)
        out['leadingJet_corrected_pt_mass'].fill(year=year, pt = pt, rho = rho_SD_corrected_pt_mass, n2 = events.N2, weight = events.weight)

        return out

    def postprocess(self, accumulator):
        return accumulator

if(__name__ == "__main__"):
    from coffea.util import load, save
    import os,glob
    import argparse
    from dask.distributed import Client
    from dask_jobqueue import HTCondorCluster

    parser = argparse.ArgumentParser()    
    parser.add_argument("--selection",default='vjets')
    parser.add_argument("--chunk",type=int,default=-1)
    parser.add_argument("--maxchunks",type=int,default=-1)
    parser.add_argument("--output","-o",type=str,default='qcd_pt_v_rho_v_n2.coffea')
    parser.add_argument("--debug",action="store_true")
    parser.add_argument("--years",nargs='+', default=['2017','UL17'])
    
    args = parser.parse_args()

    path = '/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/workdir_{SELECTION}_{YEAR}/'

    sample_pattern = os.path.join(path,'*QCD*.root')
    
    samples = {y:glob.glob(sample_pattern.format(SELECTION=args.selection,YEAR=y)) for y in args.years}

    for k,v in samples.items():
        print(k,len(v))
    
    processor_args = {}
    if(args.chunk > 0): processor_args['chunksize']=args.chunk
    if(args.maxchunks > 0): processor_args['maxchunks']=args.maxchunks

    processor_instance = DDTMapPrep()

    output_file = os.path.join(os.getcwd(),args.output)

    
    cluster = HTCondorCluster(cores=8, memory="24GB", disk="5GB")

    cluster.scale(20)
    client = Client(cluster)
    client.wait_for_workers(2)
    
    print(f"{client}\n{client.dashboard_link}")

    os.chdir(os.environ['TMPDIR'])

    ouput = processor.run_uproot_job(samples,
                                     treename = "AnalysisTree",
                                     processor_instance=processor_instance,
                                     executor = processor.dask_executor,
                                     executor_args = {
                                         'client':client,
                                         'schema':BaseSchema,
                                     },
                                     **processor_args
                                     )
    
    save(ouput, output_file)
