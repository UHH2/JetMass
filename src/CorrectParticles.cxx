#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;

CorrectParticles::CorrectParticles(TString filename){
    use_grid = true;
    // read configuration from root file
    TFile* file = new TFile(filename);
    TH1F* h_cat = (TH1F*) file->Get("categories");
    vector<TString> categories;
    for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));
    // now search for categories and read th2 hists
    for(unsigned int i=0; i<all_cat.size(); i++){
      bool catFound = false;
      for(auto cat: categories){
        if(all_cat[i] == cat && cat != "other"){
          grid.push_back((TH2F*) file->Get("grid_fit_"+cat));
          catFound = true;
        }
      }
      if(!catFound){
        grid.push_back((TH2F*) file->Get("grid_fit_other"));
      }
    }
}


CorrectParticles::CorrectParticles(){
  use_JECfactor = true;
}

bool CorrectParticles::process(uhh2::Event & event){
  // Correction for soft drop only
  vector<PFParticle>* particles = event.pfparticles;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return false;
  vector<Jet> subjets = topjets->at(0).subjets();
  for(auto subjet: subjets){
    for(const auto candInd : subjet.pfcand_indexs()){
      double factor = 1.0;
      if(use_JECfactor){
        factor = 1/topjets->at(0).JEC_factor_raw();
      }
      else if(use_grid){
        factor = GetGridFactor(particles->at(candInd));
      }
			LorentzVectorXYZE old_v4XYZE = toXYZ(particles->at(candInd).v4());
      particles->at(candInd).set_v4(toPtEtaPhi(old_v4XYZE * factor));
    }
  }
  return true;
}

double CorrectParticles::GetGridFactor(PFParticle p){
  int id = p.particleID();
  int ptbin = grid[id]->GetXaxis()->FindBin(p.pt());
  int etabin = grid[id]->GetYaxis()->FindBin(fabs(p.eta()));
  // if bin is overflow, set manually to last bin
  int Nbins_pt =  grid[id]->GetXaxis()->GetNbins();
  int Nbins_eta = grid[id]->GetYaxis()->GetNbins();
  if(ptbin > Nbins_pt) ptbin = Nbins_pt;
  if(etabin > Nbins_eta) etabin = Nbins_eta;
  double factor = grid[id]->GetBinContent(ptbin, etabin);
  return factor;
}
