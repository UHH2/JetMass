#include "UHH2/JetMass/include/UnfoldingHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

UnfoldingHists::UnfoldingHists(Context & ctx, const string & dirname, const string & reco_sel_handle_name, const string & gen_sel_handle_name):Hists(ctx, dirname){
  
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");
  std::cout << "setting up Unfolding Hists" << std::endl;

  // if (isMC) {
  //   gen_sel_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  // }
  // reco_sel_handle = ctx.get_handle<bool> (reco_sel_handle_name);
  
  std::string N2DDT_file_path = "JetMass/Histograms/ddtmaps/QCD_2017_PFMass_smooth_gaus4p00sigma.root";
  std::string N2DDT_hist_name = "N2_v_pT_v_rho_0p05_smooth_gaus4p00sigma_maps_cleaner_PFMass";  
  TFile* f_N2DDT = new TFile(locate_file(N2DDT_file_path).c_str(),"READ");  
  h_N2DDT = (TH2D*)f_N2DDT->Get(N2DDT_hist_name.c_str());
  
  //TUnfoldBinning * a = new TUnfoldBinning("a");
   detector_binning_msd_pt = new TUnfoldBinning("detector");
  
  // detector_binning_msd_pt->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),true,true);
  // detector_binning_msd_pt->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),true,true);
  
  // generator_binning_msd_pt = new TUnfoldBinning("generator");
  // generator_binning_msd_pt->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),true,true);
  // generator_binning_msd_pt->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),true,true);
  
  // TH2D* h_pt_msd_response_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_binning_msd_pt, detector_binning_msd_pt, "hist_msd_pt_response_tmp");
  // h_pt_msd_response = copy_book_th2d(h_pt_msd_response_tmp, "hist_msd_pt_response");
  
  // TH1* h_msd_detector_tmp = detector_binning_msd_pt->CreateHistogram("hist_msd_detector_tmp");
  // h_msd_detector = copy_book_th1d(h_msd_detector,"hist_msd_detector");
  
  // TH1* h_msd_generator_tmp = generator_binning_msd_pt->CreateHistogram("hist_msd_generator_tmp");
  // h_msd_generator = copy_book_th1d(h_msd_generator,"hist_msd_generator");
}

void UnfoldingHists::fill(const Event & event){
  auto weight = event.weight;
  if(event.topjets->size() <1) return;
  auto jet = event.topjets->at(0);

  // bool pass_reco = event.get(reco_sel_handle);
  // bool pass_gen(false);
  // if(isMC){
  //   pass_gen = event.get(gen_sel_handle);
  // }
  
  // if(pass_reco || pass_gen){
  //   float gen_pt(-1.0),genmsd(-1.0);
  //   float reco_pt(-1.0),reco_msd(-1.0);
  //   if(pass_reco){
  //     reco_pt = event.topjets->at(0).pt();
  //     reco_msd = event.topjets->at(0).softdropmass();
  //   }
  //   if(pass_gen){
  //     gen_pt = event.gentopjets->at(0);
  //     gen_msd = event.gentopjets->at(0).softdropmass();
  //   }
  //}
    
}
  


UnfoldingHists::~UnfoldingHists(){}

//stolen from https://github.com/raggleton/QGAnalysis/blob/master/src/QGAnalysisUnfoldHists.cxx
TH1D * UnfoldingHists::copy_book_th1d(TH1 * h, const std::string & newName) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1D>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(),
                    h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
}


TH2D * UnfoldingHists::copy_book_th2d(TH2 * h, const std::string & newName) {
  return book<TH2D>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}
