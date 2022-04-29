#include "UHH2/JetMass/include/UnfoldingHists.h"
#include "UHH2/JetMass/include/JetMassUtils.h"
#include "UHH2/core/include/Event.h"

// #include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

UnfoldingHists::UnfoldingHists(Context & ctx, const string & dirname, const std::vector<double> msd_edges, const std::vector<double> pt_edges, const string & reco_sel_handle_name, const string & gen_sel_handle_name, const string & matching_sel_handle_name, const string & recotopjet_handle_name, const string & gentopjet_handle_name):Hists(ctx, dirname){
  
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");
  dirname_=dirname;
  
  isTTbarSel = false;
  isVJetsSel = false;

  if(dirname == "unfolding_hists")debug=true;
  else debug=false;

  const std::string& selection_ = ctx.get("selection", "");
  if     (selection_ == "ttbar") isTTbarSel = true;
  else if(selection_ == "vjets")   isVJetsSel = true;
  else if(selection_ == "none") skip_filling = true;
  else throw runtime_error("UnfoldingHists: Select 'ttbar' or 'vjets' selection");

  if (isMC) {
    gen_sel_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }
  reco_sel_handle = ctx.get_handle<bool> (reco_sel_handle_name);

  matching_sel_handle = ctx.get_handle<MatchingSelection>(matching_sel_handle_name);

  recotopjet_handle = ctx.get_handle<const TopJet*>(recotopjet_handle_name);
  gentopjet_handle = ctx.get_handle<const GenTopJet*>(gentopjet_handle_name);
  
  //Binning without Underflow/Overflow
  detector_binning_msd_pt = new TUnfoldBinning("detector");  
  detector_binning_msd_pt->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),false,false);
  detector_binning_msd_pt->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),false,false);
  
  generator_binning_msd_pt = new TUnfoldBinning("generator");
  generator_binning_msd_pt->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),false,false);
  generator_binning_msd_pt->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),false,false);
  
  TH2D* h_pt_msd_response_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_binning_msd_pt, detector_binning_msd_pt, "hist_msd_pt_response_tmp");
  h_pt_msd_response = copy_book_th2d(h_pt_msd_response_tmp, "hist_msd_pt_response");
  h_pt_msd_response_matched_cat1 = copy_book_th2d(h_pt_msd_response_tmp, "hist_msd_pt_response_matched_cat1");
  h_pt_msd_response_matched_cat2 = copy_book_th2d(h_pt_msd_response_tmp, "hist_msd_pt_response_matched_cat2");
  h_pt_msd_response_unmatched = copy_book_th2d(h_pt_msd_response_tmp, "hist_msd_pt_response_unmatched");
  delete h_pt_msd_response_tmp;
  
  TH1* h_msd_detector_tmp = detector_binning_msd_pt->CreateHistogram("hist_msd_detector_tmp");
  h_msd_detector = copy_book_th1d(h_msd_detector_tmp,"hist_msd_detector");
  h_msd_detector_matched_cat1 = copy_book_th1d(h_msd_detector_tmp,"hist_msd_detector_matched_cat1");
  h_msd_detector_matched_cat2 = copy_book_th1d(h_msd_detector_tmp,"hist_msd_detector_matched_cat2");
  h_msd_detector_unmatched = copy_book_th1d(h_msd_detector_tmp,"hist_msd_detector_unmatched");
  delete h_msd_detector_tmp;
   
  TH1* h_msd_generator_tmp = generator_binning_msd_pt->CreateHistogram("hist_msd_generator_tmp");
  h_msd_generator = copy_book_th1d(h_msd_generator_tmp,"hist_msd_generator");
  h_msd_generator_matched_cat1 = copy_book_th1d(h_msd_generator_tmp,"hist_msd_generator_matched_cat1");
  h_msd_generator_matched_cat2 = copy_book_th1d(h_msd_generator_tmp,"hist_msd_generator_matched_cat2");
  h_msd_generator_unmatched = copy_book_th1d(h_msd_generator_tmp,"hist_msd_generator_unmatched");
  delete h_msd_generator_tmp;

  //Binning with Underflow/Overflow
  detector_binning_msd_pt_flow = new TUnfoldBinning("detector_flow");
  
  detector_binning_msd_pt_flow->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),true,true);
  detector_binning_msd_pt_flow->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),true,true);
  
  generator_binning_msd_pt_flow = new TUnfoldBinning("generator_flow");
  generator_binning_msd_pt_flow->AddAxis("msd",msd_edges.size()-1,msd_edges.data(),true,true);
  generator_binning_msd_pt_flow->AddAxis("pt",pt_edges.size()-1,pt_edges.data(),true,true);
  
  TH2D* h_pt_msd_response_flow_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_binning_msd_pt_flow, detector_binning_msd_pt_flow, "hist_msd_pt_response_flow_tmp");
  h_pt_msd_response_flow = copy_book_th2d(h_pt_msd_response_flow_tmp, "hist_msd_pt_response_flow");
  h_pt_msd_response_matched_cat1_flow = copy_book_th2d(h_pt_msd_response_flow_tmp, "hist_msd_pt_response_matched_cat1_flow");
  h_pt_msd_response_matched_cat2_flow = copy_book_th2d(h_pt_msd_response_flow_tmp, "hist_msd_pt_response_matched_cat2_flow");
  h_pt_msd_response_unmatched_flow = copy_book_th2d(h_pt_msd_response_flow_tmp, "hist_msd_pt_response_unmatched_flow");
  delete h_pt_msd_response_flow_tmp;
  
  TH1* h_msd_detector_flow_tmp = detector_binning_msd_pt_flow->CreateHistogram("hist_msd_detector_flow_tmp");
  h_msd_detector_flow = copy_book_th1d(h_msd_detector_flow_tmp,"hist_msd_detector_flow");
  h_msd_detector_matched_cat1_flow = copy_book_th1d(h_msd_detector_flow_tmp,"hist_msd_detector_matched_cat1_flow");
  h_msd_detector_matched_cat2_flow = copy_book_th1d(h_msd_detector_flow_tmp,"hist_msd_detector_matched_cat2_flow");
  h_msd_detector_unmatched_flow = copy_book_th1d(h_msd_detector_flow_tmp,"hist_msd_detector_unmatched_flow");
  delete h_msd_detector_flow_tmp;
   
  TH1* h_msd_generator_flow_tmp = generator_binning_msd_pt_flow->CreateHistogram("hist_msd_generator_flow_tmp");
  h_msd_generator_flow = copy_book_th1d(h_msd_generator_flow_tmp,"hist_msd_generator_flow");
  h_msd_generator_matched_cat1_flow = copy_book_th1d(h_msd_generator_flow_tmp,"hist_msd_generator_matched_cat1_flow");
  h_msd_generator_matched_cat2_flow = copy_book_th1d(h_msd_generator_flow_tmp,"hist_msd_generator_matched_cat2_flow");
  h_msd_generator_unmatched_flow = copy_book_th1d(h_msd_generator_flow_tmp,"hist_msd_generator_unmatched_flow");
  delete h_msd_generator_flow_tmp;


  h_control_gen_tau21 = book<TH1D>("control_gen_tau21","control_gen_tau21",100,-1,1);
  h_control_reco_tau21 = book<TH1D>("control_reco_tau21","control_reco_tau21",100,-1,1);

  h_control_gen_tau32 = book<TH1D>("control_gen_tau32","control_gen_tau32",100,-1,1);
  h_control_reco_tau32 = book<TH1D>("control_reco_tau32","control_reco_tau32",100,-1,1);

  h_control_gen_N2 = book<TH1D>("control_gen_N2","control_gen_N2",100,-1,1);
  h_control_reco_N2 = book<TH1D>("control_reco_N2","control_reco_N2",100,-1,1);

  h_control_gen_pt = book<TH1D>("control_gen_pt","control_gen_pt",200,0,3000);
  h_control_reco_pt = book<TH1D>("control_reco_pt","control_reco_pt",200,0,3000);

  h_control_dRVJets = book<TH2D>("control_dRVJets","control_dRVJets",200,0,6,200,0,6);

  h_control_msd_detector_bin0 = book<TH1D>("control_msd_detector_bin0","control_msd_detector_bin0",msd_edges.size()-1,msd_edges.data());
  h_control_msd_detector_bin1 = book<TH1D>("control_msd_detector_bin1","control_msd_detector_bin1",msd_edges.size()-1,msd_edges.data());
  h_control_msd_detector_bin2 = book<TH1D>("control_msd_detector_bin2","control_msd_detector_bin2",msd_edges.size()-1,msd_edges.data());
  h_control_msd_detector_bin3 = book<TH1D>("control_msd_detector_bin3","control_msd_detector_bin3",msd_edges.size()-1,msd_edges.data());

  h_control_msd_generator_bin0 = book<TH1D>("control_msd_generator_bin0","control_msd_generator_bin0",msd_edges.size()-1,msd_edges.data());
  h_control_msd_generator_bin1 = book<TH1D>("control_msd_generator_bin1","control_msd_generator_bin1",msd_edges.size()-1,msd_edges.data());
  h_control_msd_generator_bin2 = book<TH1D>("control_msd_generator_bin2","control_msd_generator_bin2",msd_edges.size()-1,msd_edges.data());
  h_control_msd_generator_bin3 = book<TH1D>("control_msd_generator_bin3","control_msd_generator_bin3",msd_edges.size()-1,msd_edges.data());


}

void UnfoldingHists::fill(const Event & event){
  //don't fill for Data
  if(!isMC)return;
  if(skip_filling)return;
  auto weight = event.weight;
  if(event.topjets->size() <1) return;
  auto jet = event.topjets->at(0);
  MatchingSelection matching_selection = event.get(matching_sel_handle);

  bool pass_reco = event.get(reco_sel_handle);
  bool pass_gen = event.get(gen_sel_handle);

  //get gentopjet from handle if valid and fallback to leading gentopjet in event if exists.
  const GenTopJet *gentopjet(NULL);
  if(event.is_valid(gentopjet_handle))gentopjet = event.get(gentopjet_handle);
  else{
    if(event.gentopjets->size()>0)gentopjet = &event.gentopjets->at(0);
  }
  //get gentopjet from handle if valid and fallback to leading topjet in event if exists.
  const TopJet *recotopjet(NULL);
  if(event.is_valid(recotopjet_handle))recotopjet = event.get(recotopjet_handle);
  else{
    if(event.topjets->size()>0)recotopjet = &event.topjets->at(0);
  }
  

  if(gentopjet==nullptr && recotopjet==nullptr){
    throw std::runtime_error("UnfoldingHists: No valid jet handles were provided and there are no leading gen- & reco-topjets in the event. Skipping UnfoldingHist filling!");
    return;
  }
  
  int reco_bin(0), gen_bin(0);
  int reco_bin_flow(0), gen_bin_flow(0);
  
  if(pass_reco || pass_gen){
    float gen_pt(-1.0),gen_msd(-1.0);
    float reco_pt(-1.0),reco_msd(-1.0);
    if(pass_reco && recotopjet){
      reco_pt = recotopjet->pt();
      reco_msd = recotopjet->softdropmass();
      reco_bin = detector_binning_msd_pt->GetGlobalBinNumber(reco_msd, reco_pt);
      reco_bin_flow = detector_binning_msd_pt_flow->GetGlobalBinNumber(reco_msd, reco_pt);

      h_control_reco_tau21->Fill(safe_tau21(recotopjet), weight);
      h_control_reco_tau32->Fill(safe_tau32(recotopjet), weight);
      h_control_reco_N2->Fill(recotopjet->ecfN2_beta1(), weight);
    }
    if(pass_gen && gentopjet){
      gen_pt = gentopjet->pt();
      gen_msd = gentopjet->softdropmass();
      // std::cout << "genjet mass: " << gen_msd << std::endl;
      gen_bin = generator_binning_msd_pt->GetGlobalBinNumber(gen_msd, gen_pt);
      gen_bin_flow = generator_binning_msd_pt_flow->GetGlobalBinNumber(gen_msd, gen_pt);

      h_control_gen_tau21->Fill(safe_tau21(gentopjet), weight);
      h_control_gen_tau32->Fill(safe_tau32(gentopjet), weight);
      h_control_gen_N2->Fill(gentopjet->ecfN2_beta1(), weight);
    }

    h_pt_msd_response->Fill(gen_bin, reco_bin, weight);
    h_pt_msd_response_flow->Fill(gen_bin_flow, reco_bin_flow, weight);
    
    if(pass_gen&&pass_reco){
      h_msd_detector->Fill(reco_bin, weight);
      h_msd_generator->Fill(gen_bin, weight);

      h_msd_detector_flow->Fill(reco_bin_flow, weight);
      h_msd_generator_flow->Fill(gen_bin_flow, weight);
      
      h_control_gen_pt->Fill(gen_pt,weight);
      h_control_reco_pt->Fill(reco_pt,weight);
    }    

    bool pass_reco_parton_matched_cat1(false),pass_reco_parton_matched_cat2(false);
    bool pass_gen_parton_matched_cat1(false),pass_gen_parton_matched_cat2(false);
    if(isVJetsSel){
      if(recotopjet)pass_reco_parton_matched_cat1 = matching_selection.passes_matching(*recotopjet,MatchingSelection::oIsMergedV) && pass_reco;
      
      if(gentopjet)pass_gen_parton_matched_cat1 = matching_selection.passes_matching(*gentopjet,MatchingSelection::oIsMergedV) && pass_gen;
      if(pass_reco_parton_matched_cat1 && pass_gen_parton_matched_cat1 ){
        const GenParticle * genV = matching_selection.get_genV();        
        float dRVreco = deltaR(*recotopjet,*genV);
        float dRVgen = deltaR(*gentopjet,*genV);
        h_control_dRVJets->Fill(dRVreco,dRVgen,weight);
      }
    }
    if(isTTbarSel){
      if(recotopjet)pass_reco_parton_matched_cat1 = matching_selection.passes_matching(*recotopjet,MatchingSelection::oIsMergedTop) && pass_reco;
      if(recotopjet)pass_reco_parton_matched_cat2 = ( matching_selection.passes_matching(*recotopjet,MatchingSelection::oIsMergedV) || matching_selection.passes_matching(*recotopjet,MatchingSelection::oIsMergedQB) ) && pass_reco;

      if(gentopjet)pass_gen_parton_matched_cat1 = matching_selection.passes_matching(*gentopjet,MatchingSelection::oIsMergedTop) && pass_gen;
      if(gentopjet)pass_gen_parton_matched_cat2 = ( matching_selection.passes_matching(*gentopjet,MatchingSelection::oIsMergedV) || matching_selection.passes_matching(*gentopjet,MatchingSelection::oIsMergedQB) ) && pass_gen;
    }
    
    if(pass_reco_parton_matched_cat1 || pass_gen_parton_matched_cat1){
      int reco_bin_matched(0),gen_bin_matched(0);
      int reco_bin_flow_matched(0),gen_bin_flow_matched(0);
      if(pass_reco_parton_matched_cat1){
        reco_bin_matched = reco_bin;
        reco_bin_flow_matched = reco_bin_flow;

      }
      if(pass_gen_parton_matched_cat1){
        gen_bin_matched = gen_bin;
        gen_bin_flow_matched = gen_bin_flow;
      }

      h_pt_msd_response_matched_cat1->Fill(gen_bin_matched, reco_bin_matched, weight);
      h_pt_msd_response_matched_cat1_flow->Fill(gen_bin_flow_matched, reco_bin_flow_matched, weight);

      if(pass_reco_parton_matched_cat1 && pass_gen_parton_matched_cat1){
        h_msd_detector_matched_cat1->Fill(reco_bin, weight);        
        h_msd_generator_matched_cat1->Fill(gen_bin_matched, weight);        

        h_msd_detector_matched_cat1_flow->Fill(reco_bin_flow, weight);
        h_msd_generator_matched_cat1_flow->Fill(gen_bin_flow, weight);


        if(reco_pt>500 && reco_pt<650)h_control_msd_detector_bin0->Fill(reco_msd,weight);
        if(reco_pt>650 && reco_pt<800)h_control_msd_detector_bin1->Fill(reco_msd,weight);
        if(reco_pt>800 && reco_pt<1200)h_control_msd_detector_bin2->Fill(reco_msd,weight);
        if(reco_pt>1200)h_control_msd_detector_bin3->Fill(reco_msd,weight);

        if(gen_pt>500 && gen_pt<650)h_control_msd_generator_bin0->Fill(gen_msd,weight);
        if(gen_pt>650 && gen_pt<800)h_control_msd_generator_bin1->Fill(gen_msd,weight);
        if(gen_pt>800 && gen_pt<1200)h_control_msd_generator_bin2->Fill(gen_msd,weight);
        if(gen_pt>1200)h_control_msd_generator_bin3->Fill(gen_msd,weight);

      }
      
    }else if((pass_reco_parton_matched_cat2 || pass_gen_parton_matched_cat2) && isTTbarSel){
      int reco_bin_matched(0),gen_bin_matched(0);
      int reco_bin_flow_matched(0),gen_bin_flow_matched(0);
      if(pass_reco_parton_matched_cat2){
        reco_bin_matched = reco_bin;
        reco_bin_flow_matched = reco_bin_flow;
      }
      if(pass_gen_parton_matched_cat2){
        gen_bin_matched = gen_bin;
        gen_bin_flow_matched = gen_bin_flow;
      }

      h_pt_msd_response_matched_cat2->Fill(gen_bin_matched, reco_bin_matched, weight);
      h_pt_msd_response_matched_cat2_flow->Fill(gen_bin_flow_matched, reco_bin_flow_matched, weight);

      if(pass_reco_parton_matched_cat2 && pass_gen_parton_matched_cat2){  
        h_msd_detector_matched_cat2->Fill(reco_bin, weight);        
        h_msd_generator_matched_cat2->Fill(gen_bin_matched, weight);        

        h_msd_detector_matched_cat2_flow->Fill(reco_bin_flow, weight);
        h_msd_generator_matched_cat2_flow->Fill(gen_bin_flow, weight);
      }
      
    }else{
      h_pt_msd_response_unmatched->Fill(gen_bin, reco_bin, weight);
      h_pt_msd_response_unmatched_flow->Fill(gen_bin_flow, reco_bin_flow, weight);

      h_msd_detector_unmatched->Fill(reco_bin, weight);
      h_msd_detector_unmatched_flow->Fill(reco_bin_flow, weight);

      h_msd_generator_unmatched->Fill(gen_bin, weight);
      h_msd_generator_unmatched_flow->Fill(gen_bin_flow, weight);
    }      
  }
  
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
