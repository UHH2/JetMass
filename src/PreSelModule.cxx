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
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
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
  std::unique_ptr<Selection> N_AK8_300_sel, N_AK8_500_sel, RHO_sel;
  std::unique_ptr<Selection> N_MUON_sel, N_ELEC_sel, TwoD_sel, bjetCloseToLepton_sel;

  std::vector<std::unique_ptr<uhh2::Hists>> hists;

  std::unique_ptr<AnalysisModule> pfparticles_jec_corrector,pf_applyPUPPI;

  std::unique_ptr<AnalysisModule> writer;


  bool is_mc, matchV, is_WSample, is_ZSample;
  bool isTopSel = false;
  bool isWSel = false;
  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

};


PreSelModule::PreSelModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";

  const std::string& channel = ctx.get("channel", "");
  if     (channel == "top") isTopSel = true;
  else if(channel == "W")   isWSel = true;
  else throw runtime_error("PreSelModule: Select 'top' or 'W' channel");

  // common modules
  MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Fall17_medium_noIso, PtEtaCut(55., 2.4));
  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));

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

  if(is_mc){
    mcSpikeKiller.reset(new MCLargeWeightKiller(ctx,
                                                infinity, // maximum allowed ratio of leading reco jet pT / generator HT
                                                infinity, // maximum allowed ratio of leading gen jet pT / generator HT
                                                infinity, // maximum allowed ratio of leading reco jet pT / Q scale
                                                2.5, // maximum allowed ratio of PU maximum pTHat / gen HT (ensures scale of PU < scale of hard interaction)
                                                3, // maximum allowed ratio of leading reco jet pT / pTHat
                                                3, // maximum allowed ratio of leading gen jet pT / pTHat
                                                "jets", // name of jet collection to be used
                                                "genjets" // name of genjet collection to be used
                          ));
  }


  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);

  // Application of Puppi weights onto pfparticles
  pf_applyPUPPI.reset(new ApplyPuppiToPF());

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
  RHO_sel.reset(new RhoSelection(-6.0, -2.1));
  N_MUON_sel.reset(new NMuonSelection(1,1));
  N_ELEC_sel.reset(new NElectronSelection(0,0));
  TwoD_sel.reset(new TwoDCut(0.4, 25));
  MET_sel.reset(new METCut(50., 100000.));
  if(isTopSel)    TRIGGER_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  else if(isWSel) TRIGGER_sel.reset(new AndSelection(ctx));
  bjetCloseToLepton_sel.reset(new NMuonBTagSelection(1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM) ));

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


  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;

  //remove MC Events with very large unphysical weights
  if(is_mc && isWSel){
    if(!mcSpikeKiller->passes(event)) return false;
  }

  // AK8 JEC
  topjetCorr->process(event);

  pf_applyPUPPI->process(event);

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

  if(!N_AK8_300_sel->passes(event)) passTOP = false;
  if(!N_MUON_sel->passes(event)) passTOP = false;
  if(!N_ELEC_sel->passes(event)) passTOP = false;
  if(!MET_sel->passes(event)) passTOP = false;
  if(!TwoD_sel->passes(event)) passTOP = false;
  if(!bjetCloseToLepton_sel->passes(event)) passTOP = false;


  if(!N_AK8_500_sel->passes(event)) passW = false;
  if(!RHO_sel->passes(event)) passW = false;


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
