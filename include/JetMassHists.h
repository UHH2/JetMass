#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>


using namespace std;

class JetMassHists: public uhh2::Hists {
public:
  JetMassHists(uhh2::Context & ctx, const std::string & dirname, const vector<double> ptbins_, const vector<double> etabins_, double variation_, TString mode = "SD");

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JetMassHists();

private:
  vector<vector<double>> GetSF(unsigned int ptbin, unsigned int etabin, TString direction);
  vector<PFParticle> VaryParticles(vector<PFParticle> oldParticles, vector<vector<double>> sf);
  double CalculateMJet(vector<PFParticle> Particles);
  bool use_SD;
  double variation;
  vector<double> ptbins, etabins;
  unsigned int Nbins_pt, Nbins_eta;
  TH1F *h_central, *h_central_mjet;
  vector<vector<TH1F*>> h_variationsUP;
  vector<vector<TH1F*>> h_variationsDOWN;
  TH1F* h_particle_pt;
  TH1F* h_particle_eta;
  TH1F* h_weights;
};
