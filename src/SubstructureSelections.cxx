#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/core/include/Event.h"

#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"

#include <stdexcept>

using namespace uhh2;

double uhh2::computeDDTValue(TopJet ak8jet,TH2F* ddtMap,bool useRho, std::vector<PFParticle> particles){
  int pt_bin = ddtMap->GetYaxis()->FindFixBin(ak8jet.pt());
  if(pt_bin > ddtMap->GetYaxis()->GetNbins()){
    pt_bin = ddtMap->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }
  double x=0;
  double mSD=0;
  if(particles.size()>0){
    LorentzVector jet_sd_v4;
    for(auto p: particles){
      jet_sd_v4 += p.v4();
    }
    mSD = jet_sd_v4.M();
  }else{
    mSD = ak8jet.softdropmass();
  }

  if(useRho){
    x = 2 * TMath::Log(ak8jet.softdropmass()/ ak8jet.pt());
  }else{
    x = mSD;
  }
  int x_bin = ddtMap->GetXaxis()->FindFixBin(x);
  if(x_bin > ddtMap->GetXaxis()->GetNbins()){
    x_bin = ddtMap->GetXaxis()->GetNbins();
  }else if(x_bin <= 0){
    x_bin = 1;
  }
  return ak8jet.ecfN2_beta1() - ddtMap->GetBinContent(x_bin,pt_bin);
}

double uhh2::computeTau21DDT(TopJet ak8jet, double slope){
  double rhoddt = 2 * TMath::Log( ak8jet.softdropmass()/ ak8jet.pt() );
  double tau21 = ak8jet.tau2()/ak8jet.tau1();
  return tau21 + slope*rhoddt;
}

N2Selection::N2Selection(double n2_cutValue_,double n2_minValue_):n2_cutValue(n2_cutValue_),n2_minValue(n2_minValue_){}

bool  N2Selection::passes(const uhh2::Event &event){
  assert(event.topjets);
  if(event.topjets->size() < 1) return false;
  auto ak8jet = event.topjets->at(0);
  double n2= ak8jet.ecfN2_beta1();
  return n2 < n2_cutValue && n2 > n2_minValue;
}

N2ddtSelection::N2ddtSelection(TH2F* ddtMap_,bool useRho_,bool usePFMass_):ddtMap(ddtMap_),useRho(useRho_),usePFMass(usePFMass_){
  std::cout << "N2ddt Selection uses " << (usePFMass ? "PFMass":"regular JetMass") << std::endl;
}

bool N2ddtSelection::passes(const uhh2::Event &event){
  assert(event.topjets);
  if(event.topjets->size() < 1) return false;
  auto ak8jet = event.topjets->at(0);
  std::vector<PFParticle> particles=  {};
  if(usePFMass){
    std::vector<PFParticle> * pfparticles = event.pfparticles;
    for(auto subjet: ak8jet.subjets()){
      for(const auto candInd: subjet.pfcand_indexs()){
        particles.push_back(pfparticles->at(candInd));
      }
    } 
  }
  double n2_ddt = computeDDTValue(ak8jet,ddtMap,useRho,particles);
  return n2_ddt < 0;
}

Tau21DDTSelection::Tau21DDTSelection(double cutValue_,double slope_):cutValue(cutValue_),slope(slope_){}

bool Tau21DDTSelection::passes(const uhh2::Event &event){
  assert(event.topjets);
  if(event.topjets->size()<1) return false;
  auto ak8jet = event.topjets->at(0);
  double tau21DDT = computeTau21DDT(ak8jet,slope);
  return tau21DDT < cutValue;
}

RhoSelection::RhoSelection(double rho_min_, double rho_max_):rho_min(rho_min_),rho_max(rho_max_){}

bool RhoSelection::passes(const uhh2::Event &event){
  assert(event.topjets);
  if(event.topjets->size()<1) return false;
  double rho=2*TMath::Log(event.topjets->at(0).softdropmass()/event.topjets->at(0).pt());
  return ( rho < -2.1 ) && ( rho > -6.0);
}
