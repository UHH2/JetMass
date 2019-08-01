#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TAxis.h>


using namespace std;

class PFHists: public uhh2::Hists {
public:
  PFHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~PFHists();

private:
  vector<double> ExtractBinning(TH2F*, TString);
  // pt
  TH1F* h_pt_all;
  TH1F* h_pt_chargedH;
  TH1F* h_pt_neutralH;
  TH1F* h_pt_gamma;
  TH1F* h_pt_other;

  // eta
  TH1F* h_eta_all;
  TH1F* h_eta_chargedH;
  TH1F* h_eta_neutralH;
  TH1F* h_eta_gamma;
  TH1F* h_eta_other;

  // eta weighted by energy
  TH1F* h_eta_Eweight_all;
  TH1F* h_eta_Eweight_chargedH;
  TH1F* h_eta_Eweight_neutralH;
  TH1F* h_eta_Eweight_gamma;
  TH1F* h_eta_Eweight_other;

  // count
  TH2F* h_count_all;
  TH2F* h_count_chargedH;
  TH2F* h_count_neutralH;
  TH2F* h_count_gamma;
  TH2F* h_count_other;

  // number of pf in jet
  TH1F* h_particle_injet;

};
