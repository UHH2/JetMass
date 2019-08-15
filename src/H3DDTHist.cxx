#include "UHH2/JetMass/include/H3DDTHist.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

H3DDTHist::H3DDTHist(Context & ctx, const string & dirname): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");
  N2_v_mSD_v_pT=book<TH3D>("N2_v_mSD_v_pT","x=m_{SD},y=p_{T},z=N_2",100,0,500,100,0,2000,100,-1.5,1);
  N2_v_rho_v_pT=book<TH3D>("N2_v_rho_v_pT","x=#rho,y=p_{T},z=N_2",100,-10,0,100,0,2000,100,-1.5,1);
  N2_v_rho=book<TH2D>("N2_v_rho","x=#rho=#log{m^2_{SD}/p^2_{T}},y=N_2^{#beta=1}",125,0.0,0.5,12,-7.0,-1.0);

}

void H3DDTHist::fill(const Event & event){
  auto weight = event.weight;
  if(event.topjets->size() <1) return;
  auto jet = event.topjets->at(0);
  N2_v_mSD_v_pT->Fill(jet.softdropmass(),jet.pt(),jet.ecfN2_beta1(),weight);
  double rho= 2* TMath::Log(jet.softdropmass()/jet.pt());
  N2_v_rho->Fill(rho,jet.ecfN2_beta1(),weight);
  N2_v_rho_v_pT->Fill(rho,jet.pt(),jet.ecfN2_beta1(),weight);
}

H3DDTHist::~H3DDTHist(){}
