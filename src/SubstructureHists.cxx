#include "UHH2/JetMass/include/SubstructureHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

SubstructureHists::SubstructureHists(Context & ctx, const string & dirname,TH2F* ddtMap_,bool useRho_): Hists(ctx, dirname),ddtMap(ddtMap_),useRho(useRho_){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");

  book<TH1D>("N2","N_{2}",40,-1.5,1);
  book<TH1D>("N2ddt","N_{2}^{DDT}",40,-1.5,1);
  book<TH1D>("tau21","#tau_{21}",40,-1.5,1);
  book<TH1D>("tau21ddt","tau_{21}^{DDT}",40,-1.5,1);
  rho_v_pt = book<TH2F>("rho_v_pt","x=#rho=#log{m^2_{SD}/p^2_{T}},y=p_{T}",12,-7.0,-1.0,100,200,1200);
}

void SubstructureHists::fill(const Event & event){
  auto weight = event.weight;
  assert(event.topjets);
  if(event.topjets->size() < 1) return;
  hist("N2")->Fill(event.topjets->at(0).ecfN2_beta1(),weight);
  hist("N2ddt")->Fill(computeDDTValue(event.topjets->at(0),ddtMap,useRho),weight);
  hist("tau21")->Fill(event.topjets->at(0).tau2()/event.topjets->at(0).tau1(),weight);
  hist("tau21ddt")->Fill(computeTau21DDT(event.topjets->at(0),0.063),weight);
  double rho = 2 * TMath::Log(event.topjets->at(0).softdropmass()/ event.topjets->at(0).pt());
  rho_v_pt->Fill(rho,event.topjets->at(0).pt(),weight);
}

SubstructureHists::~SubstructureHists(){}

TH2F* SubstructureHists::computeDDTMap(double working_point, TH3D* n2_msd_pt ){
  TH2F * ddtMap= (TH2F*)n2_msd_pt->Project3D("yx");
  char ddtMapName[100];
  sprintf(ddtMapName,"%s_ddt_wp%f",n2_msd_pt->GetName(),working_point);
  ddtMap->SetName(ddtMapName);
  ddtMap->Reset();
  ddtMap->SetStats(0);
  ddtMap->SetDirectory(0);
  int nbins_x = ddtMap->GetNbinsX();
  int nbins_y = ddtMap->GetNbinsY();
  for (int xbin = 1; xbin <= nbins_x; ++xbin) {
    for (int ybin = 1; ybin <= nbins_y; ++ybin) {
      sprintf(ddtMapName, "%s_x%d_y%d", n2_msd_pt->GetName(), xbin, ybin);
      TH1D* n2proj = n2_msd_pt->ProjectionZ(ddtMapName, xbin, xbin, ybin, ybin);
      if (n2proj->Integral() == 0) {
        double xval = n2_msd_pt->GetXaxis()->GetBinCenter(xbin);
        double yval = n2_msd_pt->GetYaxis()->GetBinCenter(ybin);
        std::cout << "N2 integral = 0 for xbin=" << xbin << " (" << xval << ") / ybin=" << ybin << " (" << yval << ") for hist " << n2_msd_pt->GetName() << ". Setting DDT to -1000." << std::endl;
        ddtMap->SetBinContent(xbin, ybin, -1000);
        continue;
      }
      double wp_array[1] = {working_point};
      double quantiles[1] = {0.};
      n2proj->GetQuantiles(1, quantiles, wp_array);
      ddtMap->SetBinContent(xbin, ybin, quantiles[0]);
      delete n2proj;
    }
  }
  return ddtMap;
}
