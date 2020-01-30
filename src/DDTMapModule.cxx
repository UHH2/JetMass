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
#include "UHH2/JetMass/include/H3DDTHist.h"

#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/TopJetCorrections.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"


#include <unistd.h>


#define EXTRAOUT false


using namespace std;
using namespace uhh2;


class DDTMapModule: public AnalysisModule {
public:

  explicit DDTMapModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;
  std::unique_ptr<TopJetCorrections> topjetCorr;

  std::unique_ptr<JetCleaner> ak4cleaner, ak4cleaner15;
  std::unique_ptr<TopJetCleaner> ak8cleaner;

  std::unique_ptr<Selection> TRIGGER_sel, MET_sel;
  std::unique_ptr<Selection> N_AK8_300_sel, N_AK8_500_sel, RHO_sel;
  std::unique_ptr<Selection> N_MUON_sel, N_ELEC_sel, TwoD_sel;
  std::unique_ptr<MatchingSelection> MatchLeadingGenW_sel;

  // std::unique_ptr<uhh2::Hists> maps_raw, maps_raw_pfmass;
  std::unique_ptr<uhh2::Hists> maps_cleaner, maps_cleaner_pfmass;
  std::unique_ptr<uhh2::Hists> maps_nak8_pt500, maps_nak8_pt500_pfmass;

  std::unique_ptr<AnalysisModule> pfparticles_jec_corrector,pf_applyPUPPI;

  bool is_mc,matchW,is_WSample;
  bool isTopSel = false;
  bool isWSel = false;

  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;
};


DDTMapModule::DDTMapModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";

  std::string version = ctx.get("dataset_version");
  is_WSample = version.find("WJets") != std::string::npos;
  matchW = version.find("WMatched") != std::string::npos;

  const std::string& channel = ctx.get("channel", "");
  if     (channel == "top") isTopSel = true;
  else if(channel == "W")   isWSel = true;
  else throw runtime_error("DDTMapModule: Select 'top' or 'W' channel");

  // common modules
  MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Fall17_medium_noIso, PtEtaCut(55., 2.4));
  common.reset(new CommonModules());
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->init(ctx);

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

  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);

  // PF correctors
  if(is_mc || !isWSel){
    pfparticles_jec_corrector.reset(new CorrectParticles());
    pf_applyPUPPI.reset(new ApplyPuppiToPF());
  }

  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.5;
  AK8_Clean_pT = 200.0;
  AK8_Clean_eta = 2.5;
  ak4cleaner15.reset(new JetCleaner(ctx, 15.0, AK4_Clean_eta));
  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));

  // SELECTIONS
  N_AK8_300_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(300.,100000.))));
  N_AK8_500_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(500.,100000.))));
  RHO_sel.reset(new RhoSelection(-6.0, -2.1));
  N_MUON_sel.reset(new NMuonSelection(1,1));
  N_ELEC_sel.reset(new NElectronSelection(0,0));
  TwoD_sel.reset(new TwoDCut(0.4, 40));
  MET_sel.reset(new METCut(50., 100000.));
  if(isTopSel)    TRIGGER_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  else if(isWSel) TRIGGER_sel.reset(new AndSelection(ctx));
  MatchLeadingGenW_sel.reset(new MatchingSelection(ctx, MatchingSelection::oIsLeadingGenW));

  // HISTOGRAMS
  // maps_raw.reset(new H3DDTHist(ctx, "maps_raw"));
  maps_cleaner.reset(new H3DDTHist(ctx, "maps_cleaner"));
  maps_nak8_pt500.reset(new H3DDTHist(ctx, "maps_NAK8_Pt500"));

  // maps_raw_pfmass.reset(new H3DDTHist(ctx, "maps_raw_PFMass"));
  maps_cleaner_pfmass.reset(new H3DDTHist(ctx, "maps_cleaner_PFMass"));
  maps_nak8_pt500_pfmass.reset(new H3DDTHist(ctx, "maps_NAK8_Pt500_PFMass"));
}

bool DDTMapModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "DDTMapModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }
  // MATCHING


  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;

  //Spikekiller
  if (is_mc){
    if (!mcSpikeKiller->passes(event)) return false;
  }

  
  // AK8 JEC
  topjetCorr->process(event);

  // if(is_mc || !isWSel) pfparticles_jec_corrector->process(event);
  if(is_mc || !isWSel) pf_applyPUPPI->process(event);

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

  // maps_raw->fill(event);
  // maps_raw_pfmass->fill(event);

  // CLEANER
  ak4cleaner15->process(event);
  ak8cleaner->process(event);

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
  if(!TwoD_sel->passes(event)) passTOP = false;
  if(!MET_sel->passes(event)) passTOP = false;

  if(!N_AK8_500_sel->passes(event)) passW = false;
  if(!RHO_sel->passes(event)) passW = false;

  if(isWSel && is_WSample && passW){
    bool Wmatched = MatchLeadingGenW_sel->passes_matching(event,event.topjets->at(0));
    if(matchW) passW = Wmatched;
    else passW = !Wmatched;
  }

  // RUN CLEANER AGAIN with pt > 30 (has to be after TwoD sel)
  ak4cleaner->process(event);
  maps_cleaner->fill(event);
  maps_cleaner_pfmass->fill(event);
  
  // FILL HISTS
  if(passTOP && isTopSel){
    maps_nak8_pt500->fill(event);
    maps_nak8_pt500_pfmass->fill(event);
  }
  else if(passW && isWSel){
    maps_nak8_pt500->fill(event);
    maps_nak8_pt500_pfmass->fill(event);
  }

  // DECIDE TO STORE EVENT
  if(!passTOP && !passW) return false;
  else return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the DDTMapModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(DDTMapModule)
