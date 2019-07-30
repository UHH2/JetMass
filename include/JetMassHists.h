#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>


using namespace std;

class JetMassHists: public uhh2::Hists {
public:
  JetMassHists(uhh2::Context & ctx, const std::string & dirname, double variation_, TString mode = "SD");

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JetMassHists();

private:
  vector<PFParticle> VaryParticles(vector<PFParticle> oldParticles, int i, int j, TString var);
  double CalculateMJet(vector<PFParticle> Particles);
  double CalculateRho(vector<PFParticle> Particles);
  bool isMC;
  bool use_SD;
  double variation;
  int Nbins_pt, Nbins_eta;
  TH1F *h_rho, *h_mass, *h_mass_jet, *h_rho_jet;
  vector<vector<TH1F*>> h_mass_UP;
  vector<vector<TH1F*>> h_mass_DOWN;
  vector<vector<TH1F*>> h_rho_UP;
  vector<vector<TH1F*>> h_rho_DOWN;
  TH1F* h_particle_pt;
  TH1F* h_particle_eta;
  TH1F* h_weights;
  TH2F* grid;
};
