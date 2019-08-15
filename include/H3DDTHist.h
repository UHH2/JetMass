#pragma once

#include "UHH2/core/include/Hists.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH3.h>


using namespace std;

class H3DDTHist: public uhh2::Hists {
public:
  H3DDTHist(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~H3DDTHist();

private:
  bool isMC;
  TH3D * N2_v_mSD_v_pT;
  TH3D * N2_v_rho_v_pT;
  TH2D * N2_v_rho;
};
