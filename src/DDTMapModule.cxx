#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h" //mainly for runPeriods vectors

#include "UHH2/JetMass/include/H3DDTHist.h"

#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/TopJetCorrections.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
#include "UHH2/common/include/MCLargeWeightKiller.h"


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

  // std::unique_ptr<Selection> N_AK8_500_sel, RHO_sel;

  std::unique_ptr<uhh2::Hists> maps_cleaner, maps_cleaner_pfmass;
  std::unique_ptr<uhh2::Hists> maps_nak8_pt500, maps_nak8_pt500_pfmass;

  std::unique_ptr<AnalysisModule> pf_applyPUPPI;

  bool is_mc;

  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;
};


DDTMapModule::DDTMapModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";

  std::string version = ctx.get("dataset_version");

  // common modules
  common.reset(new CommonModules());
  common->init(ctx);

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

  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);

  // PF correctors
  pf_applyPUPPI.reset(new ApplyPuppiToPF());

  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));

  // SELECTIONS
  // N_AK8_500_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(500.,100000.))));
  // RHO_sel.reset(new RhoSelection(-6.0, -2.1));

  // HISTOGRAMS
  maps_cleaner.reset(new H3DDTHist(ctx, "maps_cleaner"));
  // maps_nak8_pt500.reset(new H3DDTHist(ctx, "maps_NAK8_Pt500"));

  maps_cleaner_pfmass.reset(new H3DDTHist(ctx, "maps_cleaner_PFMass"));
  // maps_nak8_pt500_pfmass.reset(new H3DDTHist(ctx, "maps_NAK8_Pt500_PFMass"));

}

bool DDTMapModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "DDTMapModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }
  // MATCHING
  assert(event.pfparticles);
  assert(event.genparticles);

  
  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;
  //Spikekiller
  if (is_mc){
    if (!mcSpikeKiller->passes(event)) return false;
  }

  
  // AK8 JEC
  topjetCorr->process(event);

  pf_applyPUPPI->process(event);

  sort_by_pt<TopJet>(*event.topjets);

  // CLEANER
  ak8cleaner->process(event);

  // SELECTIONS
  // bool passSelections = true;

  // if(!N_AK8_500_sel->passes(event)) passSelections = false;
  // if(!RHO_sel->passes(event)) passSelections = false;

  
  // FILL HISTS
  maps_cleaner->fill(event);
  maps_cleaner_pfmass->fill(event);

  // if(passSelections){
  //   maps_nak8_pt500->fill(event);
  //   maps_nak8_pt500_pfmass->fill(event);
  // }

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the DDTMapModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(DDTMapModule)
