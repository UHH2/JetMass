#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h" //mainly for runPeriods vectors

#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MuonHists.h"

#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/TopJetCorrections.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/common/include/MCLargeWeightKiller.h"
#include "UHH2/JetMass/include/WriteOutput.h"


#include <unistd.h>

#include "TFile.h"
#include "TH1F.h"

#define EXTRAOUT false


using namespace std;
using namespace uhh2;

Particle GetGeneratorW(Event & event);

class PreSelModule: public AnalysisModule {
public:

  explicit PreSelModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;
  std::unique_ptr<TopJetCorrections> topjetCorr;

  std::unique_ptr<JetCleaner> ak4cleaner, ak4cleaner15;
  std::unique_ptr<TopJetCleaner> ak8cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> ak8cleaner_dRlep;

  std::unique_ptr<Selection> TRIGGER_sel, MET_sel;
  std::unique_ptr<Selection> N_AK8_300_sel, N_AK8_500_sel;
  std::unique_ptr<Selection> N_MUON_sel, N_ELEC_sel, TwoD_sel, bjetCloseToLepton_sel;
  std::unique_ptr<AndSelection> LEP_VETO_sel;
  
  std::vector<std::unique_ptr<uhh2::Hists>> hists;

  std::unique_ptr<AnalysisModule> writer;


  bool is_mc, is_QCD, matchV, is_WSample, is_ZSample, is_buggyPU;
  bool isTopSel = false;
  bool isWSel = false;
  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

};


PreSelModule::PreSelModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";
  is_QCD = ctx.get("dataset_version").find("QCD") != std::string::npos;
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;
  std::cout << "This sample has "<< (is_buggyPU ? "buggyPU" : "normalPU") << std::endl;
  const std::string& channel = ctx.get("channel", "");
  if     (channel == "top") isTopSel = true;
  else if(channel == "W")   isWSel = true;
  else throw runtime_error("PreSelModule: Select 'top' or 'W' channel");

  // common modules
  MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Fall17_medium_noIso, PtEtaCut(55., 2.4));
  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));


  if(is_buggyPU){
    TString pu_file_path = (TString) ctx.get("dataset_version");
    pu_file_path = pu_file_path.ReplaceAll("MC","");
    pu_file_path = pu_file_path.ReplaceAll("_buggyPU","");
    pu_file_path = pu_file_path.ReplaceAll("_test","");
    pu_file_path = "common/data/2017/Pileup_QCD_PtBinned/MyMCPileupHistogram"+pu_file_path+".root";
    ctx.set("pileup_directory",(std::string) pu_file_path);
  }
  std::cout << "reweighting mc pileup using " << ctx.get("pileup_directory")<<" as mc profile dir" <<std::endl;

  common.reset(new CommonModules());
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->set_jet_id(jetid);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection(true);
  common->switch_jetPtSorter();
  common->set_HTjetid(jetid);
  if(is_mc) common->disable_metfilters();
  common->init(ctx);

  if(is_mc && isWSel && is_QCD){
    mcSpikeKiller.reset(new MCLargeWeightKiller(ctx,
                                                2, // maximum allowed ratio of leading reco jet pT / generator HT
                                                2, // maximum allowed ratio of leading gen jet pT / generator HT
                                                2, // maximum allowed ratio of leading reco jet pT / Q scale
                                                2, // maximum allowed ratio of PU maximum pTHat / gen HT (ensures scale of PU < scale of hard interaction)
                                                2, // maximum allowed ratio of leading reco jet pT / pTHat
                                                2, // maximum allowed ratio of leading gen jet pT / pTHat
                                                "jets", // name of jet collection to be used
                                                "genjets" // name of genjet collection to be used
                          ));
  }


  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);

  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;
  // ak4cleaner15.reset(new JetCleaner(ctx, 15.0, AK4_Clean_eta));
  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  ak8cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(0.8));

  // SELECTIONS
  N_AK8_300_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(300.,100000.))));
  N_AK8_500_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(500.,100000.))));
  N_MUON_sel.reset(new NMuonSelection(1,1));
  N_ELEC_sel.reset(new NElectronSelection(0,0));
  TwoD_sel.reset(new TwoDCut(0.4, 25));
  MET_sel.reset(new METCut(50., 100000.));
  if(isTopSel)    TRIGGER_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  else if(isWSel) TRIGGER_sel.reset(new AndSelection(ctx));
  bjetCloseToLepton_sel.reset(new NMuonBTagSelection(1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM) ));

  LEP_VETO_sel.reset(new AndSelection(ctx,"lepton-veto"));
  LEP_VETO_sel->add<NElectronSelection>("ele-veto",0,0);
  LEP_VETO_sel->add<NMuonSelection>("muon-veto",0,0);
  
  // HISTOGRAMS
  hists.emplace_back(new ElectronHists(ctx, "ElectronHists"));
  hists.emplace_back(new EventHists(ctx, "EventHists"));
  hists.emplace_back(new MuonHists(ctx, "MuonHists"));
  hists.emplace_back(new JetHists(ctx, "JetHists"));
  hists.emplace_back(new TopJetHists(ctx, "TopJetHists"));

  ctx.undeclare_all_event_output();

  writer.reset(new WriteOutput(ctx));
}

bool PreSelModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "PreSelModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }
  // MATCHING

  // Throw away events with NVTX in buggy area for buggy samples
  if(is_mc && is_buggyPU){
    float n_true = event.genInfo->pileup_TrueNumInteractions();
    if(n_true < 10. || n_true > 72.) return false;
  }

  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;

  //remove MC Events with very large unphysical weights
  if(is_mc && isWSel && is_QCD){
    if(!mcSpikeKiller->passes(event)) return false;
  }

  // AK8 JEC
  topjetCorr->process(event);

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

  // CLEANER
  // ak4cleaner15->process(event);
  ak4cleaner->process(event);
  ak8cleaner->process(event);
  ak8cleaner_dRlep->process(event);

  // TRIGGER
  if(!TRIGGER_sel->passes(event)) return false;

  // SELECTIONS
  bool passTOP = true;
  bool passW = true;

  if(isTopSel) passW = false;
  else if(isWSel) passTOP = false;

  if(isTopSel){
  if(!N_AK8_300_sel->passes(event)) passTOP = false;
  if(!N_MUON_sel->passes(event)) passTOP = false;
  if(!N_ELEC_sel->passes(event)) passTOP = false;
  if(!MET_sel->passes(event)) passTOP = false;
  if(!TwoD_sel->passes(event)) passTOP = false;
  if(!bjetCloseToLepton_sel->passes(event)) passTOP = false;
  }


  if(isWSel){
  if(!LEP_VETO_sel->passes(event)) passW = false;
  if(!N_AK8_500_sel->passes(event)) passW = false;
  }


  // FILL HISTS
  if(passTOP && isTopSel)  for(auto & h: hists) h->fill(event);
  else if(passW && isWSel) for(auto & h: hists) h->fill(event);

  // DECIDE TO STORE EVENT
  if(!passTOP && !passW) return false;
  else{
    writer->process(event);
    return true;
  }
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the PreSelModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(PreSelModule)
