#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include <vector>
#include <TString.h>
#include <TFile.h>

#include <TH3.h>
#include <TH2.h>
#include "UHH2/common/include/TTbarGen.h"

using namespace std;

class JetMassGenHists: public uhh2::Hists {
public:
  JetMassGenHists(uhh2::Context & ctx, const std::string & dirname, const std::string & ttbargen_handlename, const std::string & genHT_handlename, const std::string & matched_gentopjet_handlename="");

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~JetMassGenHists();

private:
  bool isMC;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<double> h_genHT;
  uhh2::Event::Handle<const GenTopJet*> h_matched_gentopjet;

  TH2D *h_m_vs_pt,*h_m_vs_eta;
  TH2D *h_m_vs_dRreco,*h_m_vs_dRV;
};
