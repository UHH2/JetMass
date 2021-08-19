#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/TopPtReweight.h"
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

#include "UHH2/common/include/TTbarGen.h"
#include <unistd.h>

#include "TFile.h"
#include "TH1F.h"

#define EXTRAOUT false


using namespace std;
using namespace uhh2;

Particle GetGeneratorW(Event & event);

class JetMassModule: public AnalysisModule {
public:

  explicit JetMassModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> topPtReweighting;
  std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;
  std::unique_ptr<TopJetCorrections> topjetCorr;

  std::unique_ptr<JetCleaner> ak4cleaner, ak4cleaner15;
  std::unique_ptr<TopJetCleaner> ak8cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> ak8cleaner_dRlep;

  std::unique_ptr<AndSelection> reco_selection;

  std::vector<std::unique_ptr<uhh2::Hists>> hists;
  std::unique_ptr<uhh2::Hists> h_pfhists_200to500, h_pfhists_500to1000, h_pfhists_1000to2000, h_pfhists_2000to3000, h_pfhists_3000to4000, h_pfhists_4000to5000;
  std::unique_ptr<uhh2::Hists> h_pfhists_inclusive, h_pfhists_500to550, h_pfhists_550to600, h_pfhists_600to675, h_pfhists_675to800, h_pfhists_800to1200, h_pfhists_1200toInf;
  
  std::unique_ptr<AnalysisModule> writer;

  TString version;
  bool is_mc, is_QCD, matchV, is_WSample, is_ZSample, is_buggyPU;
  bool isTTbarSel = false;
  bool isVJetsSel = false;
  bool do_genStudies = false;
  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

  // uhh2::Event::Handle<double>h_mtt_gen;
  uhh2::Event::Handle<double>h_weight_pre_ttbar_reweight;

  std::string NLOWeightsDir = "UHHNtupleConverter/NLOweights";
  TH1F *h_kfactor, *h_ewcorr;
  std::unique_ptr<MatchingSelection> matching_selection;
  
};


JetMassModule::JetMassModule(Context & ctx){

  // Set some boolians
  version = ctx.get("dataset_version");
  is_mc = ctx.get("dataset_type") == "MC";
  is_QCD = ctx.get("dataset_version").find("QCD") != std::string::npos;
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;
  std::cout << "This sample has "<< (is_buggyPU ? "buggyPU" : "normalPU") << std::endl;
  const std::string& selection_ = ctx.get("selection", "");
  if     (selection_ == "ttbar") isTTbarSel = true;
  else if(selection_ == "vjets")   isVJetsSel = true;
  else throw runtime_error("JetMassModule: Select 'ttbar' or 'vjets' selection");

  matching_selection.reset(new MatchingSelection(ctx));
  if(isVJetsSel && (version.Contains("WJets") || version.Contains("ZJets"))){
    std::string NLOWeightsFilename = NLOWeightsDir + (std::string)(version.Contains("W") ? "/WJets" : "/ZJets") + "Corr.root";
    TFile * NLOWeightsFile = new TFile(locate_file(NLOWeightsFilename).c_str());
    h_kfactor = (TH1F*) NLOWeightsFile->Get("kfactor");
    h_ewcorr = (TH1F*) NLOWeightsFile->Get("ewcorr");
  }
  
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


  ctx.undeclare_all_event_output();
  
  h_weight_pre_ttbar_reweight = ctx.declare_event_output<double>("weight_pre_ttbar_reweight");  
  topPtReweighting.reset(new TopPtReweight(ctx, 0.0615, -0.0005,"","",true));
  
  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);
  topjetCorr->disable_jersmear();

  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;

  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  ak8cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(0.8));

  reco_selection.reset(new AndSelection(ctx,"reco_selection"));
  if(isTTbarSel){
    reco_selection->add<TriggerSelection>("Trigger selection","HLT_Mu50_v*");
    reco_selection->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 200 GeV", 1,-1,TopJetId(PtEtaCut(200.,100000.)));
    reco_selection->add<NMuonSelection>("N_{#mu} #geq 1", 1,1);
    reco_selection->add<NElectronSelection>("N_{e} = 0", 0,0);
    reco_selection->add<METCut>("MET > 50 GeV", 50.,100000.);
    reco_selection->add<TwoDCut>("2D-Cut",0.4,25.);
    reco_selection->add<NMuonBTagSelection>("b-jet in muon hemisphere", 1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM));
  }else if(isVJetsSel){
    reco_selection->add<NElectronSelection>("ele-veto",0,0);
    reco_selection->add<NMuonSelection>("muon-veto",0,0);
    reco_selection->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 500 GeV", 1,-1,TopJetId(PtEtaCut(500.,100000.)));
    reco_selection->add<HTCut>("H_{T} > 1000 GeV",ctx, 1000.);  
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
  

  writer.reset(new WriteOutput(ctx));
}

bool JetMassModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "JetMassModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }

  // if(version.Contains("TTbar_") && (version.Contains("had") || version.Contains("semilep") || version.Contains("dilep")) ){
  // float mtt = -1.0;
  // if(version.Contains("TT")){
  //   TTbarGen ttgen(*event.genparticles);
  //   mtt = (ttgen.Top().v4()+ttgen.Antitop().v4()).M();
  //   if(version.Contains("had") || version.Contains("semilep") || version.Contains("dilep")){
  //     if(mtt>700) return false;
  //   }
  // }
  // event.set(h_mtt_gen,mtt);
  
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

  event.set(h_weight_pre_ttbar_reweight,event.weight);
  if(is_mc)topPtReweighting->process(event);
  
  // AK8 JEC
  topjetCorr->process(event);

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

  // CLEANER
  // ak4cleaner15->process(event);
  ak4cleaner->process(event);
  ak8cleaner->process(event);
  ak8cleaner_dRlep->process(event);

  if(event.topjets->size()>0){
    //k-factors
    if(isVJetsSel){
      matching_selection->init(event);
      bool IsMergedWZ  = matching_selection->passes_matching(event.topjets->at(0),MatchingSelection::oIsMergedV);
      const GenJet * closest_genjet_1 = closestParticle(event.topjets->at(0), *event.genjets);
      const GenJet * closest_genjet_2 = event.topjets->size() > 1 ? closestParticle(event.topjets->at(1), *event.genjets) : closest_genjet_1;
      
      float gen_pt_1 = closest_genjet_1 ? closest_genjet_1->pt() : -9999;
      float gen_pt_2 = closest_genjet_2 ? closest_genjet_2->pt() : -9999;
      float genjetpt = IsMergedWZ ? gen_pt_1 : gen_pt_2;
      
      double kfactor_pt = genjetpt;
      double ewk_pt = genjetpt;
      
      if( kfactor_pt > 3000 ) kfactor_pt = 2800;
      if( kfactor_pt < 200 ) kfactor_pt = 205;
      
      float kfactor_bin = h_kfactor->GetXaxis()->FindBin(kfactor_pt);
      
      float w= h_kfactor->GetBinContent(kfactor_bin);
      
      if( ewk_pt > 1205 ) ewk_pt = 1205;
      if( ewk_pt < 160 ) ewk_pt = 165;
      
      float ewk_bin = h_ewcorr->GetXaxis()->FindBin(ewk_pt);

      float w_ew= h_ewcorr->GetBinContent(ewk_bin);
      float nlo_weight = w * w_ew;
     
      event.weight *= nlo_weight ;
    }
    
    //PFHists
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
  
  bool   pass_reco_selection = reco_selection->passes(event);
  
  // make sure there is a closest gentopjet to topjet is closer than dR<0.6  
  if(is_mc && do_genStudies && pass_reco_selection){
    auto dR = numeric_limits<double>::infinity();
    if(event.gentopjets->size()>0){
      const GenTopJet * closest_gentopjet = closestParticle(event.topjets->at(0), *event.gentopjets);
      dR = deltaR(event.topjets->at(0),*closest_gentopjet);
    }
    if(dR>0.6){
      pass_reco_selection = false;
    }
  }
  if(pass_reco_selection){  
    // FILL HISTS
    for(auto & h: hists) h->fill(event);
    // STORE EVENT
    writer->process(event);
    return true;
  }else return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the JetMassModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(JetMassModule)
