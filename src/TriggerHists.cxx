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

  //IsoMuon as reference
  probe_triggers = {"PFJet400","PFJet450","PFJet500","PFJet550","AK8PFJet400","AK8PFJet450","AK8PFJet500","AK8PFJet550","PFHT800","PFHT890","PFHT900","PFHT1050"};
  for(unsigned int i=0;i<probe_triggers.size();i++){
    TString common_name = "IsoMuonReference_"+probe_triggers[i]+"_";
    book<TH1F>("HT_"+common_name+"denom", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("HT_"+common_name+"num", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("AK8_PT_"+common_name+"denom", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK8_PT_"+common_name+"num", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"denom", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"num", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    common_name += "prescale_";
    book<TH1F>("HT_"+common_name+"denom", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("HT_"+common_name+"num", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("AK8_PT_"+common_name+"denom", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK8_PT_"+common_name+"num", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"denom", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"num", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
  }
  for(unsigned int i=0;i<probe_triggers.size();i++){
    TString common_name = "IsoMuonReference_"+probe_triggers[i]+"_lookup_";
    book<TH1F>("HT_"+common_name+"denom", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("HT_"+common_name+"num", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("AK8_PT_"+common_name+"denom", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK8_PT_"+common_name+"num", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"denom", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"num", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    common_name += "prescale_";
    book<TH1F>("HT_"+common_name+"denom", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("HT_"+common_name+"num", "H_{T}^{AK8PFJets} [GeV/c]", 300,0,3000);
    book<TH1F>("AK8_PT_"+common_name+"denom", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK8_PT_"+common_name+"num", "p_{T}^{AK8PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"denom", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
    book<TH1F>("AK4_PT_"+common_name+"num", "p_{T}^{AK4PFJets} [GeV/c]", 300,0,1500);
  }

  h_ht = ctx.get_handle<double>("HT");
  jet_tight_id = JetPFID(JetPFID::WP_TIGHT_CHS);
  is_SingleMuon = ctx.get("dataset_version").find("SingleMuon") != std::string::npos;
  nmuon_sel.reset(new NMuonSelection(1,-1, muid_fortriggereff));
}  

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void TriggerHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  if(! event.muons) return;
  if(is_SingleMuon && event.muons->size() < 1) return;
  
  
    
  std::vector<bool> trigger_lookups = {};
  std::vector<bool> trigger_bits = {};
  std::vector<int> trigger_prescales = {};
  
  for(unsigned int i=0; i<probe_triggers.size();i++){
    std::string TriggerName("HLT_"+probe_triggers[i]+"_v*");
    auto trigger_index = event.get_trigger_index(TriggerName);
    bool trigger_lookup = event.lookup_trigger_index(trigger_index);
    trigger_lookups.push_back(trigger_lookup);
    trigger_bits.push_back(trigger_lookup ? event.passes_trigger(trigger_index) : false);
    // std::cout << TriggerName << " lookup" << (trigger_lookup ? "True" : "False") << " bit:" <<  trigger_bits[i] << std::endl;
    int prescale(1);
    if(!isMC){
      prescale = event.lookup_trigger_index(trigger_index) ? event.trigger_prescale(trigger_index) : 0;
    }
    trigger_prescales.push_back(prescale);
  }
  
  float min_dR = 0.8;
  bool lepton_veto = true;
  for(const auto& muo : *event.muons){
    if(uhh2::deltaR(topjets->at(0), muo) < min_dR) lepton_veto = false;
  }
  if(event.electrons){
    for(const auto& ele : *event.electrons){
      if(uhh2::deltaR(topjets->at(0), ele) < min_dR) lepton_veto = false;
    }
  }
  bool jet_criteria = JetPFID(JetPFID::WP_TIGHT_PUPPI)(topjets->at(0), event) && lepton_veto; 
  // single muon reference trigger
  auto IsoMu27 = event.get_trigger_index("HLT_IsoMu27_v*");
  bool pass_IsoMu27 = event.lookup_trigger_index(IsoMu27) ? event.passes_trigger(IsoMu27) : false;  
  bool passing_IsoMu_reference = pass_IsoMu27 && jet_criteria && nmuon_sel->passes(event);
    
  // std::cout << "passing ref probe " <<  (passing_IsoMu_reference ? "True" : "False") <<  " " << (pass_AK8PFJET450 ? "True" : "False")  << std::endl;

  float HT(-1.0f), AK8_PT(-1.0f), AK4_PT(-1.0f);
  if(event.is_valid(h_ht)) HT = event.get(h_ht);
  if(event.topjets->size()>0) AK8_PT = event.topjets->at(0).pt();
  if(event.jets->size()>0) AK4_PT = event.jets->at(0).pt();

  hist("HT")->Fill(HT,weight);
  hist("AK8_PT")->Fill(AK8_PT,weight);
  hist("AK4_PT")->Fill(AK4_PT,weight);


 //IsoMuon as reference
  for(unsigned int i=0;i<probe_triggers.size();i++){
    TString common_name = "IsoMuonReference_"+probe_triggers[i] + "_";
    if(passing_IsoMu_reference){ 
      hist("HT_"+common_name+"denom")->Fill(HT, weight);
      hist("AK8_PT_"+common_name+"denom")->Fill(AK8_PT, weight);
      hist("AK4_PT_"+common_name+"denom")->Fill(AK4_PT, weight);
      hist("HT_"+common_name+"prescale_denom")->Fill(HT, weight*trigger_prescales[i]);
      hist("AK8_PT_"+common_name+"prescale_denom")->Fill(AK8_PT, weight*trigger_prescales[i]);
      hist("AK4_PT_"+common_name+"prescale_denom")->Fill(AK4_PT, weight*trigger_prescales[i]);
      if(trigger_bits[i]){
        hist("HT_"+common_name+"num")->Fill(HT, weight);
        hist("AK8_PT_"+common_name+"num")->Fill(AK8_PT, weight);
        hist("AK4_PT_"+common_name+"num")->Fill(AK4_PT, weight);
        hist("HT_"+common_name+"prescale_num")->Fill(HT, weight*trigger_prescales[i]);
        hist("AK8_PT_"+common_name+"prescale_num")->Fill(AK8_PT, weight*trigger_prescales[i]);
        hist("AK4_PT_"+common_name+"prescale_num")->Fill(AK4_PT, weight*trigger_prescales[i]);
      }
    }
  }
 //IsoMuon as reference with check on existence of probetrigger
   for(unsigned int i=0;i<probe_triggers.size();i++){
    TString common_name = "IsoMuonReference_"+probe_triggers[i]+"_lookup_";
    if(trigger_lookups[i]){
      if(passing_IsoMu_reference){ 
        hist("HT_"+common_name+"denom")->Fill(HT, weight);
        hist("AK8_PT_"+common_name+"denom")->Fill(AK8_PT, weight);
        hist("AK4_PT_"+common_name+"denom")->Fill(AK4_PT, weight);
        hist("HT_"+common_name+"prescale_denom")->Fill(HT, weight*trigger_prescales[i]);
        hist("AK8_PT_"+common_name+"prescale_denom")->Fill(AK8_PT, weight*trigger_prescales[i]);
        hist("AK4_PT_"+common_name+"prescale_denom")->Fill(AK4_PT, weight*trigger_prescales[i]);
        if(trigger_bits[i]){
          hist("HT_"+common_name+"num")->Fill(HT, weight);
          hist("AK8_PT_"+common_name+"num")->Fill(AK8_PT, weight);
          hist("AK4_PT_"+common_name+"num")->Fill(AK4_PT, weight);
          hist("HT_"+common_name+"prescale_num")->Fill(HT, weight*trigger_prescales[i]);
          hist("AK8_PT_"+common_name+"prescale_num")->Fill(AK8_PT, weight*trigger_prescales[i]);
          hist("AK4_PT_"+common_name+"prescale_num")->Fill(AK4_PT, weight*trigger_prescales[i]);
        }
      }
    }
  }


}


TriggerHists::~TriggerHists(){}
