#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

PFHists::PFHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  TString gfilename = ctx.get("GridFile");
  TFile* gfile = new TFile(gfilename);
  TH2F* grid = (TH2F*) gfile->Get("grid");
  h_count_all = (TH2F*) grid->Clone("h_count_all");
  h_count_chargedH = (TH2F*) grid->Clone("h_count_chargedH");
  h_count_neutralH = (TH2F*) grid->Clone("h_count_neutralH");
  h_count_gamma = (TH2F*) grid->Clone("h_count_gamma");
  h_count_other = (TH2F*) grid->Clone("h_count_other");
  h_count_all->Reset();
  h_count_chargedH->Reset();
  h_count_neutralH->Reset();
  h_count_gamma->Reset();
  h_count_other->Reset();
  h_pt_all = book<TH1F>("pt_all", "p_T", 50, 0, 50);
  h_pt_chargedH = book<TH1F>("pt_chargedH", "p_T", 50, 0, 50);
  h_pt_neutralH = book<TH1F>("pt_neutralH", "p_T", 50, 0, 50);
  h_pt_gamma = book<TH1F>("pt_gamma", "p_T", 50, 0, 50);
  h_pt_other = book<TH1F>("pt_other", "p_T", 50, 0, 50);
  h_eta_all = book<TH1F>("eta_all", "eta", 50, -5, 5);
  h_eta_chargedH = book<TH1F>("eta_chargedH", "eta", 50, -5, 5);
  h_eta_neutralH = book<TH1F>("eta_neutralH", "eta", 50, -5, 5);
  h_eta_gamma = book<TH1F>("eta_gamma", "eta", 50, -5, 5);
  h_eta_other = book<TH1F>("eta_other", "eta", 50, -5, 5);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void PFHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  vector<PFParticle>* allparticles = event.pfparticles;


  // get the PFParticles inside the first topjet
  // if mode is soft drop, get particles from subjets
  vector<PFParticle> particles;
  vector<Jet> subjets = topjets->at(0).subjets();
  for(auto subjet: subjets){
    for(const auto candInd : subjet.pfcand_indexs()){
      particles.push_back(allparticles->at(candInd));
    }
  }

  // fill some particle histograms
  for(auto p: particles){
    h_pt_all->Fill(p.pt(), weight);
    h_eta_all->Fill(p.eta(), weight);
    h_count_all->Fill(p.pt(), p.eta(), weight);
    /*
    Particle IDs
    0 - undefined
    1 - charged hadron
    2 - elec
    3 - muon
    4 - gamma
    5 - neutral hadron
    6 - HF tower identified as a hadron
    7 - HF tower identified as an EM particle
    */
    int id = p.particleID();
    if(id == 1){
      h_pt_chargedH->Fill(p.pt(), weight);
      h_eta_chargedH->Fill(p.eta(), weight);
      h_count_chargedH->Fill(p.pt(), p.eta(), weight);
    }
    else if(id == 4){
      h_pt_gamma->Fill(p.pt(), weight);
      h_eta_gamma->Fill(p.eta(), weight);
      h_count_gamma->Fill(p.pt(), p.eta(), weight);
    }
    else if(id == 5){
      h_pt_neutralH->Fill(p.pt(), weight);
      h_eta_neutralH->Fill(p.eta(), weight);
      h_count_neutralH->Fill(p.pt(), p.eta(), weight);
    }
    else{
      h_pt_other->Fill(p.pt(), weight);
      h_eta_other->Fill(p.eta(), weight);
      h_count_other->Fill(p.pt(), p.eta(), weight);
    }
  }

}


PFHists::~PFHists(){}
