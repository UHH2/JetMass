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
#include "UHH2/JetMass/include/PFHists.h"

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

  std::unique_ptr<AndSelection> full_selection;

  std::vector<std::unique_ptr<uhh2::Hists>> hists;
  std::unique_ptr<uhh2::Hists> h_pfhists_200to500, h_pfhists_500to1000, h_pfhists_1000to2000, h_pfhists_2000to3000, h_pfhists_3000to4000, h_pfhists_4000to5000;
  std::unique_ptr<uhh2::Hists> h_pfhists_inclusive, h_pfhists_500to550, h_pfhists_550to600, h_pfhists_600to675, h_pfhists_675to800, h_pfhists_800to1200, h_pfhists_1200toInf;
  
  std::unique_ptr<AnalysisModule> writer;


  bool is_mc, is_QCD, matchV, is_WSample, is_ZSample, is_buggyPU;
  bool isTTbarSel = false;
  bool isVJetsSel = false;
  bool do_genStudies = false;
  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

};


PreSelModule::PreSelModule(Context & ctx){

  // Set some boolians
  is_mc = ctx.get("dataset_type") == "MC";
  is_QCD = ctx.get("dataset_version").find("QCD") != std::string::npos;
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;
  std::cout << "This sample has "<< (is_buggyPU ? "buggyPU" : "normalPU") << std::endl;
  const std::string& selection_ = ctx.get("selection", "");
  if     (selection_ == "ttbar") isTTbarSel = true;
  else if(selection_ == "vjets")   isVJetsSel = true;
  else throw runtime_error("PreSelModule: Select 'ttbar' or 'vjets' selection");

  do_genStudies = string2bool(ctx.get("doGenStudies", "true"));
  
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

  if(is_mc && isVJetsSel && is_QCD){
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

  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  ak8cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(0.8));

  full_selection.reset(new AndSelection(ctx,"full_selection"));
  if(isTTbarSel){
    full_selection->add<TriggerSelection>("Trigger selection","HLT_Mu50_v*");
    full_selection->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 200 GeV", 1,-1,TopJetId(PtEtaCut(200.,100000.)));
    full_selection->add<NMuonSelection>("N_{#mu} #geq 1", 1,1);
    full_selection->add<NElectronSelection>("N_{e} = 0", 0,0);
    full_selection->add<METCut>("MET > 50 GeV", 50.,100000.);
    full_selection->add<TwoDCut>("2D-Cut",0.4,25.);
    full_selection->add<NMuonBTagSelection>("b-jet in muon hemisphere", 1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM));
  }else if(isVJetsSel){
    full_selection->add<NElectronSelection>("ele-veto",0,0);
    full_selection->add<NMuonSelection>("muon-veto",0,0);
    full_selection->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 500 GeV", 1,-1,TopJetId(PtEtaCut(500.,100000.)));
    full_selection->add<HTCut>("H_{T} > 1000 GeV",ctx, 1000.);  
  }

  
  // HISTOGRAMS
  hists.emplace_back(new ElectronHists(ctx, "ElectronHists"));
  hists.emplace_back(new EventHists(ctx, "EventHists"));
  hists.emplace_back(new MuonHists(ctx, "MuonHists"));
  hists.emplace_back(new JetHists(ctx, "JetHists"));
  hists.emplace_back(new TopJetHists(ctx, "TopJetHists"));

  h_pfhists_200to500.reset(new PFHists(ctx, "PFHists_200to500"));
  h_pfhists_500to1000.reset(new PFHists(ctx, "PFHists_500to1000"));
  h_pfhists_1000to2000.reset(new PFHists(ctx, "PFHists_1000to2000"));
  h_pfhists_2000to3000.reset(new PFHists(ctx, "PFHists_2000to3000"));
  h_pfhists_3000to4000.reset(new PFHists(ctx, "PFHists_3000to4000"));
  h_pfhists_4000to5000.reset(new PFHists(ctx, "PFHists_4000to5000"));

  h_pfhists_inclusive.reset(new PFHists(ctx, "PFHists_inclusive"));
  h_pfhists_500to550.reset(new PFHists(ctx, "PFHists_500to550"));
  h_pfhists_550to600.reset(new PFHists(ctx, "PFHists_550to600"));
  h_pfhists_600to675.reset(new PFHists(ctx, "PFHists_600to675"));
  h_pfhists_675to800.reset(new PFHists(ctx, "PFHists_675to800"));
  h_pfhists_800to1200.reset(new PFHists(ctx, "PFHists_800to1200"));
  h_pfhists_1200toInf.reset(new PFHists(ctx, "PFHists_1200toInf"));
  
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
  if(is_mc && isVJetsSel && is_QCD){
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

  //PFHists
  if(event.topjets->size()>0){
    float AK8_pt = event.topjets->at(0).pt();
    h_pfhists_inclusive->fill(event);
    if(AK8_pt>200 && AK8_pt<500)h_pfhists_200to500->fill(event);
    if(AK8_pt>500 && AK8_pt<1000)h_pfhists_500to1000->fill(event);
    if(AK8_pt>1000 && AK8_pt<2000)h_pfhists_1000to2000->fill(event);
    if(AK8_pt>2000 && AK8_pt<3000)h_pfhists_2000to3000->fill(event);
    if(AK8_pt>3000 && AK8_pt<4000)h_pfhists_3000to4000->fill(event);
    if(AK8_pt>4000 && AK8_pt<5000)h_pfhists_4000to5000->fill(event);
    if(AK8_pt>500 && AK8_pt<550)h_pfhists_500to550->fill(event);
    if(AK8_pt>550 && AK8_pt<600)h_pfhists_550to600->fill(event);
    if(AK8_pt>600 && AK8_pt<675)h_pfhists_600to675->fill(event);
    if(AK8_pt>675 && AK8_pt<800)h_pfhists_675to800->fill(event);
    if(AK8_pt>800 && AK8_pt<1200)h_pfhists_800to1200->fill(event);
    if(AK8_pt>1200)h_pfhists_1200toInf->fill(event);
  }
  
  bool   pass_full_selection = full_selection->passes(event);
  
  // make sure there is a closest gentopjet to topjet is closer than dR<0.6  
  if(is_mc && do_genStudies && pass_full_selection){
    auto dR = numeric_limits<double>::infinity();
    if(event.gentopjets->size()>0){
      const GenTopJet * closest_gentopjet = closestParticle(event.topjets->at(0), *event.gentopjets);
      dR = deltaR(event.topjets->at(0),*closest_gentopjet);
    }
    if(dR>0.6){
      pass_full_selection = false;
    }
  }
  if(pass_full_selection){  
    // FILL HISTS
    for(auto & h: hists) h->fill(event);
    // STORE EVENT
    writer->process(event);
    return true;
  }else return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the PreSelModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(PreSelModule)
