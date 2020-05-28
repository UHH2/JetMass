#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

PFHists::PFHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  TString gfilename = ctx.get("GridFile");
  TFile* gfile = new TFile(locate_file((std::string)gfilename).c_str());
  TH2F* grid = (TH2F*) gfile->Get("grid");
  vector<double> binsX = ExtractBinning(grid, "x");
  vector<double> binsY = ExtractBinning(grid, "y");

  h_count_all = book<TH2F>("count_all", "x=pt y=eta", binsX.size()-1, &binsX[0], binsY.size()-1, &binsY[0]);
  h_count_chargedH = book<TH2F>("count_chargedH", "x=pt y=eta", binsX.size()-1, &binsX[0], binsY.size()-1, &binsY[0]);
  h_count_neutralH = book<TH2F>("count_neutralH", "x=pt y=eta", binsX.size()-1, &binsX[0], binsY.size()-1, &binsY[0]);
  h_count_gamma = book<TH2F>("count_gamma", "x=pt y=eta", binsX.size()-1, &binsX[0], binsY.size()-1, &binsY[0]);
  h_count_other = book<TH2F>("count_other", "x=pt y=eta", binsX.size()-1, &binsX[0], binsY.size()-1, &binsY[0]);
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
  h_eta_Eweight_all = book<TH1F>("eta_Eweight_all", "eta", 50, -5, 5);
  h_eta_Eweight_chargedH = book<TH1F>("eta_Eweight_chargedH", "eta", 50, -5, 5);
  h_eta_Eweight_neutralH = book<TH1F>("eta_Eweight_neutralH", "eta", 50, -5, 5);
  h_eta_Eweight_gamma = book<TH1F>("eta_Eweight_gamma", "eta", 50, -5, 5);
  h_eta_Eweight_other = book<TH1F>("eta_Eweight_other", "eta", 50, -5, 5);
  h_particle_injet = book<TH1F>("particle_injet", "number particles", 30, 0, 300);
  h_density_energy = book<TH1F>("h_density_energy", "energy density", 50, 0, 500);
  h_density_number = book<TH1F>("h_density_number", "number density", 50, 0, 50);
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

  // Calculate densities
  // double R = 0.2;
  // vector<int> number_density;
  // vector<double> energy_density;
  // for(unsigned int i=0; i<particles.size(); i++){
  //   int Nparts = 0;
  //   double Eparts = 0.0;
  //   for(unsigned int j=0; j<particles.size(); j++){
  //     if(i==j) continue;
  //     if(deltaR(particles[i], particles[j]) < R){
  //       Nparts ++;
  //       Eparts += particles[j].energy();
  //     }
  //   }
  // h_density_number->Fill(Nparts, weight);
  // h_density_energy->Fill(Eparts, weight);
  // }

  // N particles in Jet
  h_particle_injet->Fill(particles.size(), weight);

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
      h_eta_Eweight_chargedH->Fill(p.eta(), p.energy()*weight);
      h_count_chargedH->Fill(p.pt(), p.eta(), weight);
    }
    else if(id == 4){
      h_pt_gamma->Fill(p.pt(), weight);
      h_eta_gamma->Fill(p.eta(), weight);
      h_eta_Eweight_gamma->Fill(p.eta(), p.energy()*weight);
      h_count_gamma->Fill(p.pt(), p.eta(), weight);
    }
    else if(id == 5){
      h_pt_neutralH->Fill(p.pt(), weight);
      h_eta_neutralH->Fill(p.eta(), weight);
      h_eta_Eweight_neutralH->Fill(p.eta(), p.energy()*weight);
      h_count_neutralH->Fill(p.pt(), p.eta(), weight);
    }
    else{
      h_pt_other->Fill(p.pt(), weight);
      h_eta_other->Fill(p.eta(), weight);
      h_eta_Eweight_other->Fill(p.eta(), p.energy()*weight);
      h_count_other->Fill(p.pt(), p.eta(), weight);
    }
  }

}

vector<double> PFHists::ExtractBinning(TH2F* hist, TString xy){
  vector<double> bins;
  TAxis *axis;
  if(xy == "x") axis = hist->GetXaxis();
  else          axis = hist->GetYaxis();

  // first bin lower edge:
  bins.push_back(axis->GetBinLowEdge(1));
  // now get all other upper edges
  int Nbins = axis->GetNbins();
  for(int i=1; i<=Nbins; i++) bins.push_back(axis->GetBinUpEdge(i));
  return bins;
}

PFHists::~PFHists(){}
