#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassHists::JetMassHists(Context & ctx, const string & dirname, const vector<double> ptbins_, const vector<double> etabins_, double variation_, TString mode ): Hists(ctx, dirname){

  // should SD be used
  use_SD = false;
  if(mode == "SD") use_SD = true;

  // get binnings and size of variation
  etabins = etabins_;
  ptbins = ptbins_;
  Nbins_eta = etabins.size()-1;
  Nbins_pt = ptbins.size()-1;
  variation = variation_;

  // setup jetmass hists
  TString xtitle = "jet mass";
  double min = 0;
  double max = 300;
  int Nbins = 50;

  TH1F* dummy = new TH1F("dummy", xtitle, Nbins, min, max);
  h_variationsUP.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  h_variationsDOWN.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  delete dummy;

  h_central = book<TH1F>("JetMass_central", xtitle, Nbins, min, max);
  h_central_mjet = book<TH1F>("JetMass_central_mjet", xtitle, Nbins, min, max);
  h_particle_pt = book<TH1F>("particle_pt", "p_T", 50, 0, 50);
  h_particle_eta = book<TH1F>("particle_eta", "eta", 50, -5, 5);
  h_weights = book<TH1F>("weights", "weight", 100, -1, 2);

  for(unsigned int i=0; i<Nbins_pt; i++){
    for(unsigned int j=0; j<Nbins_eta; j++){
      TString mjet_name = "mjet";
      TString bin_name = "_" + to_string(i) + "_" + to_string(j);
      h_variationsUP[i][j] = book<TH1F>(mjet_name+bin_name+"_up", xtitle, Nbins, min, max);
      h_variationsDOWN[i][j] = book<TH1F>(mjet_name+bin_name+"_down", xtitle, Nbins, min, max);
    }
  }
}


void JetMassHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  vector<PFParticle>* allparticles = event.pfparticles;

  h_weights->Fill(weight);

  // get the PFParticles inside the first topjet
  // if mode is soft drop, get particles from subjets
  vector<PFParticle> particles;
  if(use_SD){
    vector<Jet> subjets = topjets->at(0).subjets();
    for(auto subjet: subjets){
      for(const auto candInd : subjet.pfcand_indexs()){
        particles.push_back(allparticles->at(candInd));
      }
    }
  }
  else{
    for(const auto candInd : topjets->at(0).pfcand_indexs()){
      particles.push_back(allparticles->at(candInd));
    }
  }

  // fill central histograms
  double mass = CalculateMJet(particles);
  double mjet;
  if(use_SD) mjet = topjets->at(0).softdropmass();
  else       mjet = topjets->at(0).v4().M();
  h_central->Fill(mass, weight);
  h_central_mjet->Fill(mjet, weight);


  // fill some particle histograms
  for(auto p: particles){
    h_particle_pt->Fill(p.pt(), weight);
    h_particle_eta->Fill(p.eta(), weight);
  }

  // loop over every bin in pt and eta
  // in every bin, vary the grid, apply to pf particles, compute jetmas
  // and fill up/down histograms
  for(unsigned int i=0; i<Nbins_pt; i++){
    for(unsigned int j=0; j<Nbins_eta; j++){
      vector<vector<double>> sf_up = GetSF(i,j,"up");
      vector<vector<double>> sf_down = GetSF(i,j,"down");
      vector<PFParticle> new_particles_up = VaryParticles(particles, sf_up);
      vector<PFParticle> new_particles_down = VaryParticles(particles, sf_down);
      h_variationsUP[i][j]->Fill(CalculateMJet(new_particles_up), weight);
      h_variationsDOWN[i][j]->Fill(CalculateMJet(new_particles_down), weight);
    }
  }
}

// construct a grid with all entries set to 1
// for up/down variations one particular bin is variead
vector<vector<double>> JetMassHists::GetSF(unsigned int ptbin, unsigned int etabin, TString direction){
  vector<vector<double>> sf;
  for(unsigned int i=0; i<Nbins_pt; i++){
    vector<double> sf_oneptbin;
    for(unsigned int j=0; j<Nbins_eta; j++){
      if(i==ptbin && j==etabin){
        if(direction == "up") sf_oneptbin.push_back(1.0 + variation);
        else                  sf_oneptbin.push_back(1.0 - variation);
      }
      else sf_oneptbin.push_back(1.0);
    }
    sf.push_back(sf_oneptbin);
  }
  return sf;
}

// this applies the factor to PF particles according to the grid
vector<PFParticle> JetMassHists::VaryParticles(vector<PFParticle> oldParticles, vector<vector<double>> sf){
  vector<PFParticle> newParticles;
  for(auto p:oldParticles){
    double pt = p.pt();
    double eta = fabs(p.eta());
    unsigned int ptbin = 0;
    unsigned int etabin = 0;
    for(unsigned int i=0; i<Nbins_pt;i++){
      if(pt>ptbins[i]) ptbin = i;
    }
    for(unsigned int i=0; i<Nbins_eta;i++){
      if(eta>etabins[i]) etabin = i;
    }
    PFParticle newP;
    newP.set_v4(sf[ptbin][etabin] * p.v4());
    newParticles.push_back(newP);
  }
  return newParticles;
}

// calculate jet mass from vector of PF particles
double JetMassHists::CalculateMJet(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto p:Particles){
    jet_v4 += p.v4() * p.puppiWeight();
  }
  double mjet = jet_v4.M();
  return mjet;
}


JetMassHists::~JetMassHists(){}
