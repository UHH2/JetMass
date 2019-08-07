#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include "UHH2/JetMass/include/EnergyCorrelations.h"
#include "fastjet/PseudoJet.hh"


using namespace std;

class ECFHists: public uhh2::Hists {
public:
  ECFHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~ECFHists();

private:
  vector<fastjet::PseudoJet> ConvertParticles(vector<PFParticle>);
  TH1F *h_comparison;
  TH1F *h_N2_MIT;
  TH1F *h_N2_UHH;
  TH1F *h_LessThan3;

};
