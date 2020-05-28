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

  std::unique_ptr<uhh2::Hists> maps_cleaner, maps_cleaner_pfmass;

  bool is_mc;

  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;
  bool is_buggyPU;

};


DDTMapModule::DDTMapModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";

  std::string version = ctx.get("dataset_version");
  
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;

  
  if(is_buggyPU){
    TString pu_file_path = (TString) ctx.get("dataset_version");
    pu_file_path = pu_file_path.ReplaceAll("MC","");
    pu_file_path = pu_file_path.ReplaceAll("_buggyPU","");
    pu_file_path = pu_file_path.ReplaceAll("_test","");
    pu_file_path = "common/data/2017/Pileup_QCD_PtBinned/MyMCPileupHistogram"+pu_file_path+".root";
    ctx.set("pileup_directory",(std::string) pu_file_path);
  }
  std::cout << "reweighting mc pileup using " << ctx.get("pileup_directory")<<" as mc profile dir" <<std::endl;

  
  // common modules
  common.reset(new CommonModules());
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


  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));

  // HISTOGRAMS
  maps_cleaner.reset(new H3DDTHist(ctx, "maps_cleaner"));

  maps_cleaner_pfmass.reset(new H3DDTHist(ctx, "maps_cleaner_PFMass"));

}

bool DDTMapModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "DDTMapModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }
  // MATCHING
  assert(event.pfparticles);
  assert(event.genparticles);


// Throw away events with NVTX in buggy area for buggy samples
  if(is_buggyPU){
    float n_true = event.genInfo->pileup_TrueNumInteractions();
    if(n_true < 10. || n_true > 72.) return false;
  }  

  
  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;
  //Spikekiller
  if (is_mc){
    if(!mcSpikeKiller->passes(event)) return false;
  }

  
  // AK8 JEC
  topjetCorr->process(event);

  sort_by_pt<TopJet>(*event.topjets);

  // CLEANER
  ak8cleaner->process(event);
  
  // FILL HISTS
  maps_cleaner->fill(event);
  maps_cleaner_pfmass->fill(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the DDTMapModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(DDTMapModule)
