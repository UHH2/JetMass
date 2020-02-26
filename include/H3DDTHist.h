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
  bool use_PFMass;
  bool use_SD;

  // TH3D * N2_v_pT_v_rho_oldBinning; 

  TH3D *  N2_v_pT_v_rho;
  // TH3D *  N2_beta2_v_pT_v_rho;
  TH3D *  DeepBoosted_WvsQCD_v_pT_v_rho;

  /* TH3D *  N2_v_pT_v_mSD; */
  /* TH3D *  N2_beta2_v_pT_v_mSD; */
  /* TH3D *  DeepBoosted_WvsQCD_v_pT_v_mSD;  */
};
