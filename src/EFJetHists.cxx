#include "UHH2/JetMass/include/EFJetHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

EFJetHists::EFJetHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_chargedH_frac = book<TH1F>("chargedH_frac", "charged hadron energy fraction", 50, 0, 1);
  h_neutralH_frac = book<TH1F>("neutralH_frac", "neutral hadron energy fraction", 50, 0, 1);
  h_gamma_frac = book<TH1F>("gamma_frac", "gamma energy fraction", 50, 0, 1);
  h_chargedEM_frac = book<TH1F>("chargedEM_frac", "charged EM energy fraction", 50, 0, 1);
  h_chf_nhf = book<TH2F>("chf_nhf", "x=chf y=nhf", 50, 0, 1, 50, 0, 1);

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void EFJetHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;

  h_chargedH_frac->Fill(topjets->at(0).chargedHadronEnergyFraction(), weight);
  h_neutralH_frac->Fill(topjets->at(0).neutralHadronEnergyFraction(), weight);
  h_gamma_frac->Fill(topjets->at(0).photonEnergyFraction(), weight);
  h_chargedEM_frac->Fill(topjets->at(0).chargedEmEnergyFraction(), weight);
  h_chf_nhf->Fill(topjets->at(0).chargedHadronEnergyFraction(), topjets->at(0).neutralHadronEnergyFraction(), weight);

}



EFJetHists::~EFJetHists(){}
