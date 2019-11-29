#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include <vector>
#include <TString.h>
#include <TFile.h>

#include <TH3.h>
#include <TH2.h>

using namespace std;

class SubstructureHists: public uhh2::Hists {
public:
  SubstructureHists(uhh2::Context & ctx, const std::string & dirname, TH2F* ddtMap,bool useRho);

  virtual void fill(const uhh2::Event & ev) override;
  static TH2F* computeDDTMap(double working_point, TH3D* n2_msd_pt);
  virtual ~SubstructureHists();

private:
  TH2F* ddtMap;
  TH2F* rho_v_pt;

  bool useRho,isMC;

};
