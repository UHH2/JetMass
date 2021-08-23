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
  UnfoldingHists(uhh2::Context & ctx, const std::string & dirname, const string & reco_sel_handle_name, const string & gen_sel_handle_name, const string & matching_sel_handle_name);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~UnfoldingHists();

protected:
  bool isMC;
  TH2D* h_N2DDT;  
  
  uhh2::Event::Handle<bool> reco_sel_handle;
  uhh2::Event::Handle<bool> gen_sel_handle;
  uhh2::Event::Handle<MatchingSelection> matching_sel_handle;
  
  TUnfoldBinning* detector_binning_msd_pt;
  TUnfoldBinning* generator_binning_msd_pt;
  TH2D *h_pt_msd_response;
  TH1D *h_msd_detector, *h_msd_generator;

  TH2D *h_pt_msd_response_matched;
  TH1D *h_msd_detector_matched, *h_msd_generator_matched;
  
  TH2D *h_pt_msd_response_unmatched;
  TH1D *h_msd_detector_unmatched, *h_msd_generator_unmatched;
  bool isTTbarSel,isVJetsSel;

  
  
  const std::vector<double> pt_edges = {200.00,250.00,300.00,350.00,400.00,450.00,500.00,550.00,600.00,650.00,700.00,750.00,800.00,850.00,900.00,950.00,1000.00,1050.00,1100.00,1150.00,1200.00};

  const std::vector<double> msd_edges = {50.00,55.00,60.00,65.00,70.00,75.00,80.00,85.00,90.00,95.00,100.00,105.00,110.00,115.00,120.00,125.00,130.00,135.00,140.00,145.00,150.00,155.00,160.00,165.00,170.00,175.00,180.00,185.00,190.00,195.00,200.00,205.00,210.00,215.00,220.00,225.00,230.00,235.00,240.00,245.00,250.00,255.00,260.00,265.00,270.00,275.00,280.00,285.00,290.00,295.00,300.00};
  
  TH1D * copy_book_th1d(TH1 * h, const std::string & newName);
  TH2D * copy_book_th2d(TH2 * h, const std::string & newName);

};
