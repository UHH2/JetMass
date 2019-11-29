#include "UHH2/JetMass/include/JetMassHists_central_noConst.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassHists_central_noConst::JetMassHists_central_noConst(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // setup jetmass hists
  TString xtitleMass = "jet mass";
  double minMass = 50;
  double maxMass = 250;
  int nMass = 100;
  TString xtitleRho = "rho";
  double minRho = -10;
  double maxRho = 0;
  int nRho = 150;

  // book some hists without variations
  h_mass = book<TH1F>("Mass_central", xtitleMass, nMass, minMass, maxMass);
  h_rho = book<TH1F>("Rho_central", xtitleRho, nRho, minRho, maxRho);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void JetMassHists_central_noConst::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  double mass = event.topjets->at(0).softdropmass();
  double rho = 2*TMath::Log(event.topjets->at(0).softdropmass()/event.topjets->at(0).pt());
  h_mass->Fill(mass, weight);
  h_rho->Fill(rho, weight);
}

JetMassHists_central_noConst::~JetMassHists_central_noConst(){}
