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

  book<TH1F>("AK4_PT"        , "p_{T}^{PFJet,/JETHT/} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_450_320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500_320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_550_320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500_450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK4_PT_550_500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  
  
  book<TH1F>("AK8_PT"        , "p_{T}^{PFJet,/JETHT/} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_320"    , "p_{T}^{PFJet,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_450_320", "p_{T}^{PFJet,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500_320", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_550_320", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_450"    , "p_{T}^{PFJet,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500_450", "p_{T}^{PFJet,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_500"    , "p_{T}^{PFJet,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("AK8_PT_550_500", "p_{T}^{PFJet,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  
  book<TH1F>("HT"        , "H_{T}^{PFJets,/JETHT/} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_320"    , "H_{T}^{PFJets,hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_450_320", "H_{T}^{PFJets,hltSinglePFJet450|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_500_320", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_550_320", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet320} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_450"    , "H_{T}^{PFJets,hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_500_450", "H_{T}^{PFJets,hltSinglePFJet500|hltSinglePFJet450} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_500"    , "H_{T}^{PFJets,hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  book<TH1F>("HT_550_500", "H_{T}^{PFJets,hltSinglePFJet550|hltSinglePFJet500} [GeV/c]", 300 ,0 ,3000);
  

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
  
  auto PFJET320 = event.get_trigger_index("HLT_PFJet320_v*");
  auto PFJET450 = event.get_trigger_index("HLT_PFJet450_v*");
  auto PFJET500 = event.get_trigger_index("HLT_PFJet500_v*");		
  auto PFJET550 = event.get_trigger_index("HLT_PFJet550_v*");


}


TriggerHists::~TriggerHists(){}
