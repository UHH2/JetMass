#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>


using namespace std;

class JetMassHists_central: public uhh2::Hists {
public:
  JetMassHists_central(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JetMassHists_central();

private:
  double CalculateMJet(vector<PFParticle>);
  double CalculateRho(vector<PFParticle>);
  TH1F *h_rho, *h_mass;

};
