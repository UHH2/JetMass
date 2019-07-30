#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;

CorrectParticles::CorrectParticles(TH2* h){
  use_grid = true;
  grid = (TH2F*) h->Clone();
}

CorrectParticles::CorrectParticles(){
  use_JECfactor = true;
}



bool CorrectParticles::process(uhh2::Event & event){
  // Correction Constructed for soft drop only
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
      LorentzVector old_v4 = particles->at(candInd).v4();
      particles->at(candInd).set_v4(old_v4 * factor);
    }
  }
  return true;
}

double CorrectParticles::GetGridFactor(PFParticle p){
  int ptbin = grid->GetXaxis()->FindBin(p.pt());
  int etabin = grid->GetYaxis()->FindBin(p.eta());
  double factor = grid->GetBinContent(ptbin, etabin);
  return factor;
}
