#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/core/include/Event.h"

#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"

#include <stdexcept>

using namespace uhh2;

double uhh2::computeDDTValue(TopJet ak8jet,TH2F* ddtMap,bool useRho){
  int pt_bin = ddtMap->GetYaxis()->FindFixBin(ak8jet.pt());
  if(pt_bin > ddtMap->GetYaxis()->GetNbins()){
    pt_bin = ddtMap->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }
  double x=0;
  if(useRho){
    x = 2 * TMath::Log(ak8jet.softdropmass()/ ak8jet.pt());
  }else{
    x = ak8jet.softdropmass();
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
  double rhoddt = TMath::Log( (ak8jet.softdropmass()*ak8jet.softdropmass()) / ak8jet.pt() );
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

N2ddtSelection::N2ddtSelection(TH2F* ddtMap_,bool useRho_):ddtMap(ddtMap_),useRho(useRho_){}

bool N2ddtSelection::passes(const uhh2::Event &event){
  assert(event.topjets);
  if(event.topjets->size() < 1) return false;
  auto ak8jet = event.topjets->at(0);
  double n2_ddt = computeDDTValue(ak8jet,ddtMap,useRho);
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
