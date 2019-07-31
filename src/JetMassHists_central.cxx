#include "UHH2/JetMass/include/JetMassHists_central.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassHists_central::JetMassHists_central(Context & ctx, const string & dirname): Hists(ctx, dirname){
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
void JetMassHists_central::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  vector<PFParticle>* allparticles = event.pfparticles;

  // get the PFParticles inside the first topjet
  // if mode is soft drop, get particles from subjets
  vector<PFParticle> particles;
  vector<Jet> subjets = topjets->at(0).subjets();
  for(auto subjet: subjets){
    for(const auto candInd : subjet.pfcand_indexs()){
      particles.push_back(allparticles->at(candInd));
    }
  }

  // fill central histograms
  double mass = CalculateMJet(particles);
  double rho = CalculateRho(particles);
  h_mass->Fill(mass, weight);
  h_rho->Fill(rho, weight);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// calculate jet mass from vector of PF particles
double JetMassHists_central::CalculateMJet(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto p:Particles){
    jet_v4 += p.v4() * p.puppiWeight();
  }
  double mjet = jet_v4.M();
  return mjet;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// calculate jet rho from vector of PF particles
double JetMassHists_central::CalculateRho(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto p:Particles){
    jet_v4 += p.v4() * p.puppiWeight();
  }
  double mjet = jet_v4.M();
  double pt = jet_v4.Pt();
  double rho = TMath::Log(mjet*mjet/(pt*pt));
  return rho;
}


JetMassHists_central::~JetMassHists_central(){}
