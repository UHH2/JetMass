#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TAxis.h>


using namespace std;

class EFJetHists: public uhh2::Hists {
public:
  EFJetHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~EFJetHists();

private:
  TH1F* h_chargedH_frac;
  TH1F* h_neutralH_frac;
  TH1F* h_gamma_frac;
  TH1F* h_chargedEM_frac;

};
