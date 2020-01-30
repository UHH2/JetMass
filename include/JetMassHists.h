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
  JetMassHists(uhh2::Context & ctx, const std::string & dirname, TString mode = "SD");
  JetMassHists(uhh2::Context & ctx, const std::string & dirname, double variation_, TString mode = "SD");

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JetMassHists();

private:
  vector<PFParticle> VaryParticles(vector<PFParticle>, int, int, TString, TString);
  void ConstructOtherIDs();
  bool inCategory(PFParticle p, TString cat);
  double CalculateMJet(vector<PFParticle>, float JEC_factor_raw=1.);
  double CalculateRho(vector<PFParticle>, float JEC_factor_raw=1.);
  vector<double> CalculateMJetVariation(vector<PFParticle>, int, int, TString);
  bool isMC;
  bool use_SD;
  bool use_constituents;
  bool do_variations;
  bool no_JEC_on_Mass;
  double variation;
  int Nbins_pt, Nbins_eta, Nbins_cat;
  TH1F *h_rho, *h_mass, *h_mass_jet, *h_rho_jet;
  vector<vector<vector<TH1F*>>> h_mass_UP;
  vector<vector<vector<TH1F*>>> h_mass_DOWN;
  vector<vector<vector<TH2F*>>> h_pt_v_mass_UP;
  vector<vector<vector<TH2F*>>> h_pt_v_mass_DOWN;
  TH1F* h_particle_pt;
  TH1F* h_particle_eta;
  TH1F* h_weights;
  TH2F* grid;
  TH2F* h_pt_v_mass;
  vector<TString> categories;
  vector<int> otherIDs;

  /*
  Particle IDs
  0 - undefined
  1 - charged hadron
  2 - elec
  3 - muon
  4 - gamma
  5 - neutral hadron
  6 - HF tower identified as a hadron
  7 - HF tower identified as an EM particle
  */
  vector<TString> all_cat = {"undef", "chargedH", "elec", "muon", "gamma", "neutralH", "had", "em"};
};
