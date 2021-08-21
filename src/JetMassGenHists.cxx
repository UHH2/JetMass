#include "UHH2/JetMass/include/JetMassGenHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassGenHists::JetMassGenHists(Context & ctx, const string & dirname, const std::string & ttbargen_handlename): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");

  h_ttbargen = ctx.get_handle<TTbarGen>(ttbargen_handlename);
  
  book<TH1D>("N2","N_{2}",40,-1.5,1);
}

void JetMassGenHists::fill(const Event & event){
  if(!isMC)return;
  auto weight = event.weight;
  hist("N2")->Fill(0.6,weight);
  TTbarGen ttgen = event.get(h_ttbargen);
  if(ttgen.IsSemiLeptonicDecay()){
    std::cout << "lepton pt: " << ttgen.ChargedLepton().pt() << std::endl;
  }
}

JetMassGenHists::~JetMassGenHists(){}
