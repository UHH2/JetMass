#include "UHH2/JetMass/include/JetMassGenHists.h"
#include "UHH2/JetMass/include/JetMassUtils.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassGenHists::JetMassGenHists(Context & ctx, const string & dirname, const std::string & ttbargen_handlename, const std::string & genHT_handlename, const std::string & matched_gentopjet_handlename): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  auto version = ctx.get("dataset_version");

  h_ttbargen = ctx.get_handle<TTbarGen>(ttbargen_handlename);
  h_genHT = ctx.get_handle<double>(genHT_handlename);
  h_matched_gentopjet = ctx.get_handle<const GenTopJet*>(matched_gentopjet_handlename);
  
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
  
  h_m_vs_pt = book<TH2D>("m_vs_pt","m_vs_pt",300,0,300,300,0,3000);
  h_m_vs_eta = book<TH2D>("m_vs_eta","m_vs_eta",300,0,300,100,-6,6);
  h_m_vs_dRreco = book<TH2D>("m_vs_dRreco","m_vs_dRreco",300,0,300,200,0,10);  
  h_m_vs_dRV = book<TH2D>("m_vs_dRV","m_vs_dRV",300,0,300,200,0,10);  
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
  const GenTopJet *matched_gentopjet(NULL);
  if(event.is_valid(h_matched_gentopjet))matched_gentopjet = event.get(h_matched_gentopjet);
  else{
    if(ngenjets>0)matched_gentopjet = &event.gentopjets->at(0);
  }
  hist("ngenjets")->Fill(ngenjets,weight);
  if(matched_gentopjet){
    // GenTopJet genjet = event.gentopjets->at(0);
    hist("genjet_msd")->Fill(matched_gentopjet->softdropmass(), weight);
    hist("genjet_mass")->Fill(matched_gentopjet->v4().M(), weight);
    hist("genjet_pt")->Fill(matched_gentopjet->pt(), weight);
    hist("genjet_eta")->Fill(matched_gentopjet->eta(), weight);
    h_m_vs_pt->Fill(matched_gentopjet->softdropmass(),matched_gentopjet->pt(),weight);
    h_m_vs_eta->Fill(matched_gentopjet->softdropmass(),matched_gentopjet->eta(),weight);
    if(event.topjets->size()>0)h_m_vs_dRreco->Fill(matched_gentopjet->softdropmass(),deltaR(event.topjets->at(0),*matched_gentopjet),weight);
    int v_index = get_V_index(*event.genparticles);
    if(v_index>=0){
      double dRV=deltaR(event.genparticles->at(v_index),*matched_gentopjet);
      h_m_vs_dRV->Fill(matched_gentopjet->softdropmass(),dRV,weight);
    }
  }


  hist("genMET")->Fill(event.genmet->pt(),weight);

  float genHT = event.get(h_genHT);
  hist("genHT")->Fill(genHT,weight);

}

JetMassGenHists::~JetMassGenHists(){}
