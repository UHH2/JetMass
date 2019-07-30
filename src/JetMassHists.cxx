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
  TString xtitleMass = "jet mass";
  double minMass = 40;
  double maxMass = 300;
  int nMass = 130;
  TString xtitleRho = "rho";
  double minRho = -10;
  double maxRho = 0;
  int nRho = 150;

  TH1F* dummy = new TH1F("dummy", "dummy", 1, 0, 1);
  h_mass_UP.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  h_mass_DOWN.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  h_rho_UP.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  h_rho_DOWN.resize(Nbins_pt, std::vector<TH1F*>(Nbins_eta-1, dummy));
  delete dummy;

  h_mass = book<TH1F>("Mass_central", xtitleMass, nMass, minMass, maxMass);
  h_mass_jet = book<TH1F>("Mass_central_jet", xtitleMass, nMass, minMass, maxMass);
  h_rho = book<TH1F>("Rho_central", xtitleRho, nRho, minRho, maxRho);
  h_rho_jet = book<TH1F>("Rho_central_jet", xtitleRho, nRho, minRho, maxRho);
  h_particle_pt = book<TH1F>("particle_pt", "p_T", 50, 0, 50);
  h_particle_eta = book<TH1F>("particle_eta", "eta", 50, -5, 5);
  h_weights = book<TH1F>("weights", "weight", 100, -1, 2);

  for(unsigned int i=0; i<Nbins_pt; i++){
    for(unsigned int j=0; j<Nbins_eta; j++){
      TString mass_name = "Mass";
      TString rho_name = "Rho";
      TString bin_name = "_" + to_string(i) + "_" + to_string(j);
      h_mass_UP[i][j] = book<TH1F>(mass_name+bin_name+"_up", xtitleMass, nMass, minMass, maxMass);
      h_mass_DOWN[i][j] = book<TH1F>(mass_name+bin_name+"_down", xtitleMass, nMass, minMass, maxMass);
      h_rho_UP[i][j] = book<TH1F>(rho_name+bin_name+"_up", xtitleRho, nRho, minRho, maxRho);
      h_rho_DOWN[i][j] = book<TH1F>(rho_name+bin_name+"_down", xtitleRho, nRho, minRho, maxRho);
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
  double rho = CalculateRho(particles);
  double mjet;
  if(use_SD) mjet = topjets->at(0).softdropmass();
  else       mjet = topjets->at(0).v4().M();
  double ptjet = topjets->at(0).v4().Pt();
  double rhojet = TMath::Log(mjet*mjet/(ptjet*ptjet));
  h_mass->Fill(mass, weight);
  h_mass_jet->Fill(mjet, weight);
  h_rho->Fill(rho, weight);
  h_rho_jet->Fill(rhojet, weight);


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
      h_mass_UP[i][j]->Fill(CalculateMJet(new_particles_up), weight);
      h_mass_DOWN[i][j]->Fill(CalculateMJet(new_particles_down), weight);
      h_rho_UP[i][j]->Fill(CalculateRho(new_particles_up), weight);
      h_rho_DOWN[i][j]->Fill(CalculateRho(new_particles_down), weight);
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

// calculate jet rho from vector of PF particles
double JetMassHists::CalculateRho(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto p:Particles){
    jet_v4 += p.v4() * p.puppiWeight();
  }
  double mjet = jet_v4.M();
  double pt = jet_v4.Pt();
  double rho = TMath::Log(mjet*mjet/(pt*pt));
  return rho;
}


JetMassHists::~JetMassHists(){}
