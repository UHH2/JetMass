#!/usr/bin/env python3

import os
from collections import OrderedDict


uhh2datasets_path="/nfs/dust/cms/user/albrechs/UHH2/UHH2-Datasets"
import sys
sys.path.append(uhh2datasets_path)
from CrossSectionHelper import MCSampleValuesHelper


selection_names = {
    "vjets":"vjets",
    "ttbar":"ttbar",
}


sample_lists ={
    "vjets":{
        "DATA":[
            "JetHT_RunB",
            "JetHT_RunC",
            "JetHT_RunD",
            "JetHT_RunE",
            "JetHT_RunF",
        ],
        "MC":[
            "QCD_HT50to100",
            "QCD_HT100to200",
            "QCD_HT200to300",
            "QCD_HT300to500",
            "QCD_HT500to700",
            "QCD_HT700to1000",
            "QCD_HT1000to1500",
            "QCD_HT1500to2000",
            "QCD_HT2000toInf",
            "WJetsToQQ_HT400to600",
            "WJetsToQQ_HT600to800",
            "WJetsToQQ_HT800toInf",
            "ZJetsToQQ_HT400to600",
            "ZJetsToQQ_HT600to800",
            "ZJetsToQQ_HT800toInf",           
            "TTToHadronic",
            "TTToSemiLeptonic",
            "ST_tW_antitop_5f_inclusiveDecays",
            "ST_tW_top_5f_inclusiveDecays",
        ]
    },
    "ttbar":{
        "DATA":[
            "SingleMuon_RunB",
            "SingleMuon_RunC",
            "SingleMuon_RunD",
            "SingleMuon_RunE",
            "SingleMuon_RunF",
        ],
        "MC":[
            "ST_s-channel_4f_leptonDecays",
            "ST_t-channel_antitop_4f_InclusiveDecays",
            "ST_t-channel_top_4f_InclusiveDecays",
            "ST_tW_antitop_5f_NoFullyHadronicDecays_PDFWeights",
            "ST_tW_top_5f_NoFullyHadronicDecays_PDFWeights",
            "QCD_Pt-15To20_MuEnrichedPt5",
            "QCD_Pt-20To30_MuEnrichedPt5",
            "QCD_Pt-30To50_MuEnrichedPt5",
            "QCD_Pt-50To80_MuEnrichedPt5",
            "QCD_Pt-80To120_MuEnrichedPt5",
            "QCD_Pt-120To170_MuEnrichedPt5",
            "QCD_Pt-170To300_MuEnrichedPt5",
            "QCD_Pt-300To470_MuEnrichedPt5",
            "QCD_Pt-470To600_MuEnrichedPt5",
            "QCD_Pt-600To800_MuEnrichedPt5",
            "QCD_Pt-800To1000_MuEnrichedPt5",
            "QCD_Pt-1000_MuEnrichedPt5",
            "DYJetsToLL_M-50_HT-100to200",
            "DYJetsToLL_M-50_HT-1200to2500",
            "DYJetsToLL_M-50_HT-200to400",
            "DYJetsToLL_M-50_HT-2500toInf",
            "DYJetsToLL_M-50_HT-400to600",
            "DYJetsToLL_M-50_HT-600to800",
            "DYJetsToLL_M-50_HT-70to100",
            "DYJetsToLL_M-50_HT-800to1200",
            "WJetsToLNu_HT-70to100",
            "WJetsToLNu_HT-100to200",
            "WJetsToLNu_HT-200to400",
            "WJetsToLNu_HT-400to600",
            "WJetsToLNu_HT-600to800",
            "WJetsToLNu_HT-800to1200",
            "WJetsToLNu_HT-1200to2500",
            "WJetsToLNu_HT-2500toInf",
        ]
    }
}


lumi_files = {
    # "UL16":"Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_normtag.root",
    "UL16postVFP":"Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root",
    "UL16preVFP":"Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root",
    "UL17":"Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root",
    "UL18":"Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root",
}

class Sample(object):
    def __init__(self,xml,lumi=1.0):
        self.name=""
        self.xml=xml
        self.lumi=float(lumi)
    def __str__(self):
        return f"{self.name}({self.lumi}:{self.xml})"
    def __repr__(self):
        return f"{self.name}({self.lumi}:{self.xml})"

class SFrameConfig(object):

    def __init__(self,sample,year):
        self.samples={}
        self.setup_samples(sample_lists[sample],year)
        self.selection = selection_names.get(sample,sample)
        self.year = year
        
        self.AnalysisName = "JetMass" #i.e. name of the directory of the analysis (.../UHH2/<AnalysisName>/)
        
        self.output_directory= f"/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/{self.selection}Trees/{self.year}/"
        if(not os.path.isdir(self.output_directory)):
            os.makedirs(self.output_directory)
        self.config_dir="."
        if(not os.path.isdir(self.config_dir)):
            os.makedirs(self.config_dir)

        self.target_lumi = self.get_target_lumi()
        
        self.user_config = OrderedDict({
            "selection":self.selection,
            "pileup_directory_data":f"common/UHH2-data/{year}/MyDataPileupHistogram_{year}.root",
            "pileup_directory":f"common/UHH2-data/{year}/MyMCPileupHistogram_{year}.root",
            "lumi_file":f"common/UHH2-data/{year}/{lumi_files[year]}",
        })

        self.mail = "steffen.albrecht@desy.de"
        self.file_split = 5
        self.memory = 4
        self.disk = 4

    def get_target_lumi(self):
        brilcalc_output = os.path.join(os.environ['CMSSW_BASE'],
                                       "src/UHH2/common/UHH2-data/",
                                       self.year,
                                       f"TempOut{self.year}.csv")
        with open(brilcalc_output) as brilcalc_result:
            summary = brilcalc_result.read().split("\n")
            
            return float(summary[-4].split(",")[5])/1e6
            
            
        
    def setup_samples(self,samples,year):
        samples_db = MCSampleValuesHelper()
        self.samples = {
            "DATA":{
                s:Sample(samples_db.get_xml(s,"13TeV",year))
                for s in samples["DATA"]
            },
            "MC":{
                s:Sample(samples_db.get_xml(s,"13TeV",year), samples_db.get_lumi(s,"13TeV",year))
                for s in samples["MC"]
            }
        }

    def build_xml(self):
        import lxml.etree as ET
        import lxml.builder as B
        F = B.ElementMaker()

        entities = "\n".join([f"<!ENTITY \t {name} \t SYSTEM \t \"{os.path.join(uhh2datasets_path,sample.xml)}\">" for name,sample in self.samples["DATA"].items()])
        entities += "\n\n"
        entities += "\n".join([f"<!ENTITY \t {name} \t SYSTEM \t \"{os.path.join(uhh2datasets_path,sample.xml)}\">" for name,sample in self.samples["MC"].items()])
        entities += "\n\n"
        entities += f"<!ENTITY \t OUTdir \t \"{self.output_directory}\">\n"
        doctype_string=f"<!DOCTYPE JobConfiguration PUBLIC \"\" \"JobConfig.dtd\" [\n{entities}\n]>"


        doctype_string += ('\n\n<!--\n'
                           f'<ConfigParse NEventsBreak="0" LastBreak="0" FileSplit="{self.file_split}" />\n'
                           f'<ConfigSGE RAM ="{self.memory}" DISK ="{self.disk}" Mail="{self.mail}" Notification="as" Workdir="workdir_{self.selection}_{self.year}"/>\n'
                           '-->\n\n')


        config = F.JobConfiguration(
            F.Library(Name=f"libSUHH2{self.AnalysisName}"),
            F.Package(Name=f"SUHH2{self.AnalysisName}.par"),
            F.Cycle(
                *(
                    F.InputData(
                    f"&{name};",
                    F.InputTree(Name="AnalysisTree"),
                    F.OutputTree(Name="AnalysisTree"),
                    Version=f'{name}_{self.year}', Lumi=str(sample.lumi), Type="Data", NEventsMax="-1", Cacheable="False")
                    for name,sample in self.samples['DATA'].items()
                ),
                *(
                    F.InputData(
                    f"&{name};",
                    F.InputTree(Name="AnalysisTree"),
                    F.OutputTree(Name="AnalysisTree"),
                    Version=f'{name}_{self.year}', Lumi=str(sample.lumi), Type="MC", NEventsMax="-1", Cacheable="False")
                    for name,sample in self.samples['MC'].items()
                ),

                F.UserConfig(
                    *{F.Item(Name=k, Value=v) for k,v in self.user_config.items()}
                ),
                Name="uhh2::AnalysisModuleRunner", OutputDirectory="&OUTdir;",PostFix="",TargetLumi=str(self.target_lumi)),
            JobName=f"{self.selection}Job", OutputLevel="INFO")
        s = ET.tostring(config, pretty_print=True,xml_declaration=True,encoding="UTF-8",doctype=doctype_string).decode("utf-8").encode("utf-8").decode("utf-8")

        with open(f"{self.config_dir}/{self.selection}_{self.year}.xml","w") as fout:
            fout.write(s.replace('&amp;','&'))        


if(__name__ == '__main__'):

    import argparse

    years = ["UL16preVFP","UL16postVFP","UL17","UL18"]
    samples = ["vjets","ttbar"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--year",default = "UL17", choices = years)
    parser.add_argument("--sample", default = "vjets", choices = samples )
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()


    all_years = years if args.all else [args.year]
    all_samples = samples if args.all else [args.sample]
    
    for year in years:
        for sample in samples:
            sfconfig = SFrameConfig(sample,year)
            sfconfig.user_config.update(OrderedDict({
                "GridFile":"JetMass/Histograms/grid_AllPFFlavours.root",
                "doGenStudies":"false",
                "PrimaryVertexCollection":"offlineSlimmedPrimaryVertices",
                "GenParticleCollection":"GenParticles",
                "ElectronCollection":"slimmedElectronsUSER",
                "MuonCollection":"slimmedMuonsUSER",
                "JetCollection":"jetsAk4CHS",
                "GenJetCollection":"slimmedGenJets",
                "TopJetCollection":"jetsAk8PuppiSubstructure_SoftDropPuppi",
                "GenTopJetCollection":"genjetsAk8Substructure",
                "PFParticleCollection":"PFParticles",
                "METName":"slimmedMETs",
                "genMETName":"slimmedMETs_GenMET",
                "additionalBranches":"jetsAk8CHSSubstructure_SoftDropCHSgenjetsAk8SubstructureSoftDrop",
                "chsjets":"jetsAk8CHSSubstructure_SoftDropCHS",
                "sdgenjets":"genjetsAk8SubstructureSoftDrop",
                "use_sframe_weight":"false",
                "dometfilters":"true",
                "dopvfilter":"true",
                "doSpikeKiller":"true",
                "NLOweightsDir":"UHHNtupleConverter/NLOweights",
                "AnalysisModule":"JetMassModule",
                
            }))
            sfconfig.build_xml()
            
    
    # helper = MCSampleValuesHelper()
    # helper.get_lumi("TTTo2L2Nu","13TeV","2018")
    # print(helper.get_xs("TTTo2L2Nu","13TeV","UL18"))
    # os.system(f'tail -n 10 {uhh2datasets_path}/{helper.get_xml("TTTo2L2Nu","13TeV","UL18")}')
