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
  // AK8 PFJet HLT Triggers
  book<TH1F>("HT_PFJET320"    , "H_{T}^{PFJets,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET450_PFJET320", "H_{T}^{PFJets,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET450_PFJET320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET450_PFJET320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET500_PFJET320", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET500_PFJET320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET500_PFJET320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET550_PFJET320", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET550_PFJET320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET550_PFJET320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);

  //
  book<TH1F>("HT_PFJET450"    , "H_{T}^{PFJets,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET500_PFJET450", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET500_PFJET450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET500_PFJET450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);

  //  
  book<TH1F>("HT_PFJET500"    , "H_{T}^{PFJets,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET550_PFJET500", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET550_PFJET500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET550_PFJET500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);

  // AK8 PFJet HLT Triggers
  book<TH1F>("HT_AK8PFJET320"    , "H_{T}^{AK8PFJets,hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET320"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET320"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET450_AK8PFJET320", "H_{T}^{AK8PFJets,hltSingleAK8PFJet450|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET450_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet450|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET450_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet450|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET500_AK8PFJET320", "H_{T}^{AK8PFJets,hltSingleAK8PFJet500|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET500_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET500_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET550_AK8PFJET320", "H_{T}^{AK8PFJets,hltSingleAK8PFJet550|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET550_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET550_AK8PFJET320", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|hltSingleAK8PFJet320} [GeV/c]", 300 ,0 ,3000);

  //
  book<TH1F>("HT_AK8PFJET450"    , "H_{T}^{AK8PFJets,hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET450"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET450"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET500_AK8PFJET450", "H_{T}^{AK8PFJets,hltSingleAK8PFJet500|hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET500_AK8PFJET450", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET500_AK8PFJET450", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|hltSingleAK8PFJet450} [GeV/c]", 300 ,0 ,3000);

  //  
  book<TH1F>("HT_AK8PFJET500"    , "H_{T}^{AK8PFJets,hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET500"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET500"    , "p_{T}^{AK8PFJet,hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET550_AK8PFJET500", "H_{T}^{AK8PFJets,hltSingleAK8PFJet550|hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET550_AK8PFJET500", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET550_AK8PFJET500", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|hltSingleAK8PFJet500} [GeV/c]", 300 ,0 ,3000);


  //offline selection as reference
  book<TH1F>("HT_o450"    , "H_{T}^{PFJets,offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o450"    , "p_{T}^{PFJet,offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o450"    , "p_{T}^{PFJet,offlinePT450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_o500"    , "H_{T}^{PFJets,offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o500"    , "p_{T}^{PFJet,offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o500"    , "p_{T}^{PFJet,offlinePT500} [GeV/c]", 300 ,0 ,3000);
  
  book<TH1F>("HT_o550"    , "H_{T}^{PFJets,offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_o550"    , "p_{T}^{PFJet,offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_o550"    , "p_{T}^{PFJet,offlinePT550} [GeV/c]", 300 ,0 ,3000);

  
  book<TH1F>("HT_PFJET450_o450", "H_{T}^{PFJets,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET450_o450", "p_{T}^{PFJet,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET450_o450", "p_{T}^{PFJet,hltSinglePFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET500_o500", "H_{T}^{PFJets,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET500_o500", "p_{T}^{PFJet,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET500_o500", "p_{T}^{PFJet,hltSinglePFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_PFJET550_o550", "H_{T}^{PFJets,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_PFJET550_o550", "p_{T}^{PFJet,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_PFJET550_o550", "p_{T}^{PFJet,hltSinglePFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET450_o450", "H_{T}^{AK8PFJets,hltSingleAK8PFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET450_o450", "p_{T}^{AK8PFJet,hltSingleAK8PFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET450_o450", "p_{T}^{AK8PFJet,hltSingleAK8PFJet450|offlinePT450} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET500_o500", "H_{T}^{AK8PFJets,hltSingleAK8PFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET500_o500", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET500_o500", "p_{T}^{AK8PFJet,hltSingleAK8PFJet500|offlinePT500} [GeV/c]", 300 ,0 ,3000);

  book<TH1F>("HT_AK8PFJET550_o550", "H_{T}^{AK8PFJets,hltSingleAK8PFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_AK8PFJET550_o550", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_AK8PFJET550_o550", "p_{T}^{AK8PFJet,hltSingleAK8PFJet550|offlinePT550} [GeV/c]", 300 ,0 ,3000);


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

  bool pass_PFJET320 = event.lookup_trigger_index(PFJET320) ? event.passes_trigger(PFJET320) : false;
  bool pass_PFJET450 = event.lookup_trigger_index(PFJET450) ? event.passes_trigger(PFJET450) : false;
  bool pass_PFJET500 = event.lookup_trigger_index(PFJET500) ? event.passes_trigger(PFJET500) : false;
  bool pass_PFJET550 = event.lookup_trigger_index(PFJET550) ? event.passes_trigger(PFJET550) : false;
  
  auto AK8PFJET320 = event.get_trigger_index("HLT_AK8PFJet320_v*"); // prescales?
  auto AK8PFJET450 = event.get_trigger_index("HLT_AK8PFJet450_v*");
  auto AK8PFJET500 = event.get_trigger_index("HLT_AK8PFJet500_v*");
  auto AK8PFJET550 = event.get_trigger_index("HLT_AK8PFJet550_v*");

  bool pass_AK8PFJET320 = event.lookup_trigger_index(AK8PFJET320) ? event.passes_trigger(AK8PFJET320) : false;
  bool pass_AK8PFJET450 = event.lookup_trigger_index(AK8PFJET450) ? event.passes_trigger(AK8PFJET450) : false;
  bool pass_AK8PFJET500 = event.lookup_trigger_index(AK8PFJET500) ? event.passes_trigger(AK8PFJET500) : false;
  bool pass_AK8PFJET550 = event.lookup_trigger_index(AK8PFJET550) ? event.passes_trigger(AK8PFJET550) : false;


  float HT(-1.0f), AK8_PT(-1.0f), AK4_PT(-1.0f);
  if(event.is_valid(h_ht)) HT = event.get(h_ht);
  if(event.topjets->size()>0) AK8_PT = event.topjets->at(0).pt();
  if(event.jets->size()>0) AK4_PT = event.jets->at(0).pt();

  hist("HT")->Fill(HT,weight);
  hist("AK8_PT")->Fill(AK8_PT,weight);
  hist("AK4_PT")->Fill(AK4_PT,weight);

  // AK8 PFJet HLT trigger as reference
  //PFJet320 as reference selection
  if(pass_PFJET320){
    hist("HT_PFJET320")->Fill(HT,weight);
    hist("AK8_PT_PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_PFJET320 && pass_PFJET450 ){
    hist("HT_PFJET450_PFJET320")->Fill(HT,weight);
    hist("AK8_PT_PFJET450_PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET450_PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_PFJET320 && pass_PFJET500 ){
    hist("HT_PFJET500_PFJET320")->Fill(HT,weight);
    hist("AK8_PT_PFJET500_PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET500_PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_PFJET320 && pass_PFJET550 ){
    hist("HT_PFJET550_PFJET320")->Fill(HT,weight);
    hist("AK8_PT_PFJET550_PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET550_PFJET320")->Fill(AK4_PT,weight);
  }

  //PFJet450 as reference selection
  if(pass_PFJET450){
    hist("HT_PFJET450")->Fill(HT,weight);
    hist("AK8_PT_PFJET450")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET450")->Fill(AK4_PT,weight);
  }
  if(pass_PFJET450 && pass_PFJET500 ){
    hist("HT_PFJET500_PFJET450")->Fill(HT,weight);
    hist("AK8_PT_PFJET500_PFJET450")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET500_PFJET450")->Fill(AK4_PT,weight);
  }

  //PFJet500 as reference selection
  if(pass_PFJET500){
    hist("HT_PFJET500")->Fill(HT,weight);
    hist("AK8_PT_PFJET500")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET500")->Fill(AK4_PT,weight);
  }
  if(pass_PFJET500 && pass_PFJET550 ){
    hist("HT_PFJET550_PFJET500")->Fill(HT,weight);
    hist("AK8_PT_PFJET550_PFJET500")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET550_PFJET500")->Fill(AK4_PT,weight);
  }


  // AK8 PFJet HLT trigger as reference

  if(pass_AK8PFJET320){
    hist("HT_AK8PFJET320")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_AK8PFJET320 && pass_AK8PFJET450 ){
    hist("HT_AK8PFJET450_AK8PFJET320")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET450_AK8PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET450_AK8PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_AK8PFJET320 && pass_AK8PFJET500 ){
    hist("HT_AK8PFJET500_AK8PFJET320")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET500_AK8PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET500_AK8PFJET320")->Fill(AK4_PT,weight);
  }
  if(pass_AK8PFJET320 && pass_AK8PFJET550 ){
    hist("HT_AK8PFJET550_AK8PFJET320")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET550_AK8PFJET320")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET550_AK8PFJET320")->Fill(AK4_PT,weight);
  }

  //PFJet450 as reference selection
  if(pass_AK8PFJET450){
    hist("HT_AK8PFJET450")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET450")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET450")->Fill(AK4_PT,weight);
  }
  if(pass_AK8PFJET450 && pass_AK8PFJET500 ){
    hist("HT_AK8PFJET500_AK8PFJET450")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET500_AK8PFJET450")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET500_AK8PFJET450")->Fill(AK4_PT,weight);
  }

  //PFJet500 as reference selection
  if(pass_AK8PFJET500){
    hist("HT_AK8PFJET500")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET500")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET500")->Fill(AK4_PT,weight);
  }
  if(pass_AK8PFJET500 && pass_AK8PFJET550 ){
    hist("HT_AK8PFJET550_AK8PFJET500")->Fill(HT,weight);
    hist("AK8_PT_AK8PFJET550_AK8PFJET500")->Fill(AK8_PT,weight);
    hist("AK4_PT_AK8PFJET550_AK8PFJET500")->Fill(AK4_PT,weight);
  }


  // offline cuts as reference!

  //pT_offline>450 as reference selection
  if((AK8_PT>450.)){
    hist("HT_o450")->Fill(HT,weight);
    hist("AK8_PT_o450")->Fill(AK8_PT,weight);
    hist("AK4_PT_o450")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>450.) && pass_PFJET450 ){
    hist("HT_PFJET450_o450")->Fill(HT,weight);
    hist("AK8_PT_PFJET450_o450")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET450_o450")->Fill(AK4_PT,weight);
  }

  //pT_offline>500 as reference selection
  if((AK8_PT>500.)){
    hist("HT_o500")->Fill(HT,weight);
    hist("AK8_PT_o500")->Fill(AK8_PT,weight);
    hist("AK4_PT_o500")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>500.) && pass_PFJET500 ){
    hist("HT_PFJET500_o500")->Fill(HT,weight);
    hist("AK8_PT_PFJET500_o500")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET500_o500")->Fill(AK4_PT,weight);
  }

  //pT_offline>550 as reference selection
  if((AK8_PT>550.)){
    hist("HT_o550")->Fill(HT,weight);
    hist("AK8_PT_o550")->Fill(AK8_PT,weight);
    hist("AK4_PT_o550")->Fill(AK4_PT,weight);
  }
  if( (AK8_PT>550.) && pass_PFJET550 ){
    hist("HT_PFJET550_o550")->Fill(HT,weight);
    hist("AK8_PT_PFJET550_o550")->Fill(AK8_PT,weight);
    hist("AK4_PT_PFJET550_o550")->Fill(AK4_PT,weight);
  }

}


TriggerHists::~TriggerHists(){}
