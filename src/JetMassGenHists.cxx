#include "UHH2/JetMass/include/JetMassGenHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassGenHists::JetMassGenHists(Context & ctx, const string & dirname, const std::string & ttbargen_handlename, const std::string & genHT_handlename): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");

  h_ttbargen = ctx.get_handle<TTbarGen>(ttbargen_handlename);
  h_genHT = ctx.get_handle<double>(genHT_handlename);
  
  book<TH1D>("genjet_msd","m_{SD,gen AK8}",300,0,300);
  book<TH1D>("genjet_mass","m_{gen AK8}",300,0,300);
  book<TH1D>("genjet_pt","p_{T,gen AK8}",300,0,3000);
  book<TH1D>("genjet_eta","#eta_{gen AK8}",100,-6,6);
  book<TH1D>("ngenjets","N_{gen AK8}",20,0,20);
  
  book<TH1D>("genHT","H_{T,gen}",300,0,3000);
  book<TH1D>("genMET","MET_{gen}",100,0,1000);
  
  book<TH1D>("genlepton_pt","p_{T,gen lepton}",300,0,3000);
  book<TH1D>("genlepton_eta","#eta_{gen lepton}",100,-6,6);
  book<TH1D>("genlepton_pdgId","pdgId_{gen lepton}",40,-20,20);

  
}

void JetMassGenHists::fill(const Event & event){
  if(!isMC)return;
  auto weight = event.weight;

  TTbarGen ttgen = event.get(h_ttbargen);
  if(ttgen.IsSemiLeptonicDecay()){
    GenParticle lepton = ttgen.ChargedLepton();
    hist("genlepton_pt")->Fill(lepton.pt(),weight);
    hist("genlepton_eta")->Fill(lepton.eta(),weight);
    hist("genlepton_pdgId")->Fill(lepton.pdgId(),weight);
  }

  int ngenjets = event.gentopjets->size();
  hist("ngenjets")->Fill(ngenjets,weight);
  if(ngenjets>0){
    GenTopJet genjet = event.gentopjets->at(0);
    hist("genjet_msd")->Fill(genjet.softdropmass(), weight);
    hist("genjet_mass")->Fill(genjet.v4().M(), weight);
    hist("genjet_pt")->Fill(genjet.pt(), weight);
    hist("genjet_eta")->Fill(genjet.eta(), weight);
  }

  hist("genMET")->Fill(event.genmet->pt(),weight);

  float genHT = event.get(h_genHT);
  hist("genHT")->Fill(genHT,weight);
  
}

JetMassGenHists::~JetMassGenHists(){}
