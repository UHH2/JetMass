#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/common/include/Utils.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include "TUnfoldBinning.h"
#include "UHH2/JetMass/include/MatchingSelections.h"

using namespace std;

class UnfoldingHists: public uhh2::Hists {
public:
  UnfoldingHists(uhh2::Context & ctx, const std::string & dirname, const std::vector<double> msd_edges, const std::vector<double> pt_edges, const string & reco_sel_handle_name, const string & gen_sel_handle_name, const string & matching_sel_handle_name, const string & recotopjet_handle_name="", const string & gentopjet_handle_name="");

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~UnfoldingHists();

protected:
  bool isMC;
  bool debug;
  TString dirname_;
  uhh2::Event::Handle<bool> reco_sel_handle;
  uhh2::Event::Handle<bool> gen_sel_handle;
  uhh2::Event::Handle<MatchingSelection> matching_sel_handle;
  uhh2::Event::Handle<const TopJet*> recotopjet_handle;
  uhh2::Event::Handle<const GenTopJet*> gentopjet_handle;
  
  //without under-/over-flow bins
  TUnfoldBinning* detector_binning_msd_pt;
  TUnfoldBinning* generator_binning_msd_pt;

  TH2D *h_pt_msd_response;
  TH1D *h_msd_detector, *h_msd_generator;

  TH2D *h_pt_msd_response_matched_cat1;
  TH1D *h_msd_detector_matched_cat1, *h_msd_generator_matched_cat1;

  TH2D *h_pt_msd_response_matched_cat2;
  TH1D *h_msd_detector_matched_cat2, *h_msd_generator_matched_cat2;
  
  TH2D *h_pt_msd_response_unmatched;
  TH1D *h_msd_detector_unmatched, *h_msd_generator_unmatched;

  //with under-/over-flow bins
  TUnfoldBinning* detector_binning_msd_pt_flow;
  TUnfoldBinning* generator_binning_msd_pt_flow;

  TH2D *h_pt_msd_response_flow;
  TH1D *h_msd_detector_flow, *h_msd_generator_flow;

  TH2D *h_pt_msd_response_matched_cat1_flow;
  TH1D *h_msd_detector_matched_cat1_flow, *h_msd_generator_matched_cat1_flow;

  TH2D *h_pt_msd_response_matched_cat2_flow;
  TH1D *h_msd_detector_matched_cat2_flow, *h_msd_generator_matched_cat2_flow;
  
  TH2D *h_pt_msd_response_unmatched_flow;
  TH1D *h_msd_detector_unmatched_flow, *h_msd_generator_unmatched_flow;

  TH1D * h_control_gen_tau21,* h_control_reco_tau21;
  TH1D * h_control_gen_tau32,* h_control_reco_tau32;
  TH1D * h_control_gen_N2,* h_control_reco_N2;

  bool isTTbarSel,isVJetsSel;
  
  TH1D * copy_book_th1d(TH1 * h, const std::string & newName);
  TH2D * copy_book_th2d(TH2 * h, const std::string & newName);

};
