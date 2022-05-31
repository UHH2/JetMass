#include "UHH2/JetMass/include/TriggerHists.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;


TriggerHists::TriggerHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";

  book<TH1F>("HT"        , "H_{T}^{PFJets,/JETHT/} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT"        , "p_{T}^{PFJet,/JETHT/} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT"        , "p_{T}^{PFJet,/JETHT/} [GeV/c]", 300 ,0 ,3000);

  //lower threshold trigger as reference
  book<TH1F>("HT_320"    , "H_{T}^{PFJets,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_450_320", "H_{T}^{PFJets,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_450_320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_450_320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_500_320", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500_320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500_320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_550_320", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_550_320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_550_320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  //
  book<TH1F>("HT_450"    , "H_{T}^{PFJets,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_500_450", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500_450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500_450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);

  //  
  book<TH1F>("HT_500"    , "H_{T}^{PFJets,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_550_500", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_550_500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_550_500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);


  //offline selection as reference
  book<TH1F>("HT_o450"    , "H_{T}^{PFJets,offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o450"    , "p_{T}^{PFJet,offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o450"    , "p_{T}^{PFJet,offlinePT450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_450_o450", "H_{T}^{PFJets,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_450_o450", "p_{T}^{PFJet,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_450_o450", "p_{T}^{PFJet,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_o500"    , "H_{T}^{PFJets,offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o500"    , "p_{T}^{PFJet,offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o500"    , "p_{T}^{PFJet,offlinePT500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_500_o500", "H_{T}^{PFJets,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500_o500", "p_{T}^{PFJet,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500_o500", "p_{T}^{PFJet,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_o550"    , "H_{T}^{PFJets,offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o550"    , "p_{T}^{PFJet,offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o550"    , "p_{T}^{PFJet,offlinePT550} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_550_o550", "H_{T}^{PFJets,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_550_o550", "p_{T}^{PFJet,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_550_o550", "p_{T}^{PFJet,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  

  // "triggerObjects_hltSinglePFJet500"
  // "triggerObjects_hltSinglePFJet450"
  //handle_trigger220_HF = ctx.declare_event_input< vector< FlavorParticle > >("triggerObjects_hltDiPFJetAve220ForHFJEC" );
  // h_lowTrObjects = ctx.declare_event_input<vector<FlavorParticle>>("triggerObjects_hltSinglePFJet450");
  // h_highTrObjects = ctx.declare_event_input<vector<FlavorParticle>>("triggerObjects_hltSinglePFJet500");
  h_ht = ctx.get_handle<double>("HT");}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void TriggerHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  
  auto PFJET320 = event.get_trigger_index("HLT_PFJet320_v*");//prescale = 20
  auto PFJET450 = event.get_trigger_index("HLT_PFJet450_v*");//prescale = 1 ff.
  auto PFJET500 = event.get_trigger_index("HLT_PFJet500_v*");		
  auto PFJET550 = event.get_trigger_index("HLT_PFJet550_v*");
  // std::cout << "550 trigger prescale: " << event.trigger_prescale(PFJET550) << std::endl;

  float HT(-1.0f), AK8_PT(-1.0f), AK4_PT(-1.0f);
  if(event.is_valid(h_ht)) HT = event.get(h_ht);
  if(event.topjets->size()>0) AK8_PT = event.topjets->at(0).pt();
  if(event.jets->size()>0) AK4_PT = event.jets->at(0).pt();

  hist("HT")->Fill(HT,weight);
  hist("AK8_PT")->Fill(AK8_PT,weight);
  hist("AK4_PT")->Fill(AK4_PT,weight);

  //PFJet320 as reference selection
  if(event.passes_trigger(PFJET320)){
    hist("HT_320")->Fill(HT,weight);
    hist("AK8_PT_320")->Fill(AK8_PT,weight);
    hist("AK4_PT_320")->Fill(AK4_PT,weight);
  }
  if(event.passes_trigger(PFJET320) && event.passes_trigger(PFJET450) ){
    hist("HT_450_320")->Fill(HT,weight);
    hist("AK8_PT_450_320")->Fill(AK8_PT,weight);
    hist("AK4_PT_450_320")->Fill(AK4_PT,weight);
  }
  if(event.passes_trigger(PFJET320) && event.passes_trigger(PFJET500) ){
    hist("HT_500_320")->Fill(HT,weight);
    hist("AK8_PT_500_320")->Fill(AK8_PT,weight);
    hist("AK4_PT_500_320")->Fill(AK4_PT,weight);
  }
  if(event.passes_trigger(PFJET320) && event.passes_trigger(PFJET550) ){
    hist("HT_550_320")->Fill(HT,weight);
    hist("AK8_PT_550_320")->Fill(AK8_PT,weight);
    hist("AK4_PT_550_320")->Fill(AK4_PT,weight);
  }

  //PFJet450 as reference selection
  if(event.passes_trigger(PFJET450)){
    hist("HT_450")->Fill(HT,weight);
    hist("AK8_PT_450")->Fill(AK8_PT,weight);
    hist("AK4_PT_450")->Fill(AK4_PT,weight);
  }
  if(event.passes_trigger(PFJET450) && event.passes_trigger(PFJET500) ){
    hist("HT_500_450")->Fill(HT,weight);
    hist("AK8_PT_500_450")->Fill(AK8_PT,weight);
    hist("AK4_PT_500_450")->Fill(AK4_PT,weight);
  }

  //PFJet500 as reference selection
  if(event.passes_trigger(PFJET500)){
    hist("HT_500")->Fill(HT,weight);
    hist("AK8_PT_500")->Fill(AK8_PT,weight);
    hist("AK4_PT_500")->Fill(AK4_PT,weight);
  }
  if(event.passes_trigger(PFJET500) && event.passes_trigger(PFJET550) ){
    hist("HT_550_500")->Fill(HT,weight);
    hist("AK8_PT_550_500")->Fill(AK8_PT,weight);
    hist("AK4_PT_550_500")->Fill(AK4_PT,weight);
  }


  //pT_offline>450 as reference selection
  if((AK8_PT>450.)){
    hist("HT_o450")->Fill(HT,weight);
    hist("AK8_PT_o450")->Fill(AK8_PT,weight);
    hist("AK4_PT_o450")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>450.) && event.passes_trigger(PFJET450) ){
    hist("HT_450_o450")->Fill(HT,weight);
    hist("AK8_PT_450_o450")->Fill(AK8_PT,weight);
    hist("AK4_PT_450_o450")->Fill(AK4_PT,weight);
  }

  //pT_offline>500 as reference selection
  if((AK8_PT>500.)){
    hist("HT_o500")->Fill(HT,weight);
    hist("AK8_PT_o500")->Fill(AK8_PT,weight);
    hist("AK4_PT_o500")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>500.) && event.passes_trigger(PFJET500) ){
    hist("HT_500_o500")->Fill(HT,weight);
    hist("AK8_PT_500_o500")->Fill(AK8_PT,weight);
    hist("AK4_PT_500_o500")->Fill(AK4_PT,weight);
  }

  //pT_offline>550 as reference selection
  if((AK8_PT>550.)){
    hist("HT_o550")->Fill(HT,weight);
    hist("AK8_PT_o550")->Fill(AK8_PT,weight);
    hist("AK4_PT_o550")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>550.) && event.passes_trigger(PFJET550) ){
    hist("HT_550_o550")->Fill(HT,weight);
    hist("AK8_PT_550_o550")->Fill(AK8_PT,weight);
    hist("AK4_PT_550_o550")->Fill(AK4_PT,weight);
  }

}


TriggerHists::~TriggerHists(){}
