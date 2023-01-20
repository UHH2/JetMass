#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarReconstruction.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;

template<typename GenericJet>
bool JetIdSelection<GenericJet>::passes(const uhh2::Event & event){
  if(!event.is_valid(jet_handle)) throw std::runtime_error("JetIdSelection: jet handle is invalid!");
  const GenericJet * jet = event.get(jet_handle);
  if(!jet)return false;
  return jet_id(*jet,event);
}


template<typename GenericTopJet>
bool RhoCut<GenericTopJet>::passes(const uhh2::Event &event){
  if(!event.is_valid(jet_handle))return false;
  const GenericTopJet *jet = event.get(jet_handle);
  if(jet == nullptr) return false;
  double rho=2*TMath::Log(jet->softdropmass()/jet->pt());
  return ( rho < rho_max ) && ( rho > rho_min);
}

TopJetPtCut::TopJetPtCut(uhh2::Context & ctx, float pt_min_, float pt_max_):
  pt_min(pt_min_), pt_max(pt_max_){
  h_jecfactor_nominal = ctx.get_handle<double>("jecfactor_nominal");
  h_jecfactor_up = ctx.get_handle<double>("jecfactor_up");
  h_jecfactor_down = ctx.get_handle<double>("jecfactor_down");
}

bool TopJetPtCut::passes(const Event & event){
  if(event.topjets->size()==0) return false;
  float jecfactor_nominal = 1./event.topjets->at(0).JEC_factor_raw();
  LorentzVector leadingjet = event.topjets->at(0).v4();
  LorentzVector leadingjet_raw = event.topjets->at(0).v4()/jecfactor_nominal;

  if(leadingjet.pt() < pt_min) return false;
  if(leadingjet.pt() > pt_max) return false;

  if(event.is_valid(h_jecfactor_nominal)){
    float jecfactor_nominal_standalone = event.get(h_jecfactor_nominal);
    if( fabs(jecfactor_nominal-jecfactor_nominal_standalone)/jecfactor_nominal > 0.001 ){
      std::cout << "ATTENTION: jecfactor and jecfactor_nominal from standalone jet corrector differ when performing pt cut!!!"<< std::endl;
      std::cout << "jec from generic corrector: " << jecfactor_nominal << std::endl;
      std::cout << "jec from standalone corrector: " << jecfactor_nominal_standalone << std::endl;
      return false;
    }
  }
  if(event.is_valid(h_jecfactor_up)){
    float jecfactor_up = event.get(h_jecfactor_up);
    float pt = (leadingjet_raw * jecfactor_up).pt();
    if( pt < pt_min) return false;
    if( pt > pt_max) return false;
  }

  if(event.is_valid(h_jecfactor_down)){
    float jecfactor_down = event.get(h_jecfactor_down);
    float pt = (leadingjet_raw * jecfactor_down).pt();
    if( pt < pt_min) return false;
    if( pt > pt_max) return false;
  }
}

MassCut::MassCut(){}

bool MassCut::passes(const Event & event){
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() != 2) return false;
  vector<Muon>* muons = event.muons;
  vector<Electron>* electrons = event.electrons;
  if(muons->size() + electrons->size() != 1) return false;
  Particle lepton;
  if(muons->size() != 0) lepton = muons->at(0);
  else lepton = electrons->at(0);

  LorentzVector jet_lepton = lepton.v4() + topjets->at(1).v4();
  LorentzVector leading = topjets->at(0).v4();
  if(leading.M() > jet_lepton.M()) return true;
  else return false;
}

bool TwoDCut::passes(const Event & event){
  if(!event.muons->size() && !event.electrons->size()) return false;
  if(!event.jets->size()) return false;
  float drmin, ptrel;
  if(event.muons->size()) std::tie(drmin, ptrel) = drmin_pTrel(event.muons->at(0), *event.jets);
  else std::tie(drmin, ptrel) = drmin_pTrel(event.electrons->at(0), *event.jets);
  return (drmin > min_deltaR) || (ptrel > min_pTrel);
}

METCut::METCut(float minMet_, float maxMet_, bool genMet_):
minMet(minMet_), maxMet(maxMet_) , genMet(genMet_){}

bool METCut::passes(const Event & event){
  float MET(-1.0);
  if(genMet){
    assert(event.genmet);
    MET = event.genmet->pt();
  }else{
    assert(event.met);
    MET = event.met->pt();
  }
  return (MET > minMet) && (MET < maxMet);
}

NMuonBTagSelection::NMuonBTagSelection(int min_nbtag, int max_nbtag, JetId btag, double ptmin, double etamax )
{
  m_min_nbtag=min_nbtag;
  m_max_nbtag=max_nbtag;
  m_btag=btag;
  m_ptmin=ptmin;
  m_etamax=etamax;
}

bool NMuonBTagSelection::passes(const Event & event)
{
  int nbtag=0;

  //Assumes to have only one muon
  std::vector<Jet>* jets = event.jets;
  std::vector<Muon>* muons = event.muons;
  for(unsigned int i=0; i<event.jets->size(); ++i) {
    int jettagged=0;
    Jet jet=jets->at(i);
    if (m_btag(jet, event)) jettagged=1;

    if(muons->size() != 1){
      // std::cout << "ATTENTION!!! muon size " << muons->size() << std::endl;
      return false;
    }

    double deltaphi=deltaPhi(jet,muons->at(0));
    double pi = 3.14159265359;
    if(jettagged&&(deltaphi<(2*pi/3))&&(jet.pt()>m_ptmin)&&(fabs(jet.eta())<m_etamax)){

      nbtag++;

    }
  }

  if(nbtag<m_min_nbtag) return false;
  if(nbtag>m_max_nbtag) return false;
  return true;
}

HTCut::HTCut(uhh2::Context & ctx, float min_ht_,float max_ht_, std::string ht_handlename):min_ht(min_ht_),max_ht(max_ht_){
  h_ht = ctx.get_handle<double>(ht_handlename);
}

bool HTCut::passes(const Event & event){
  float ht = -2.;
  if(event.is_valid(h_ht)){
     ht = event.get(h_ht);     
  }
  return (ht > min_ht) && (ht < max_ht);
}

WToMuNuSelection::WToMuNuSelection(float min_pt_,float max_pt_):min_pt(min_pt_),max_pt(max_pt_){}

bool WToMuNuSelection::passes(const Event & event){
  assert(event.muons);
  assert(event.met);
  if(event.muons->size() < 1) return false;
  LorentzVector neutrino = NeutrinoReconstruction(event.muons->at(0).v4(), event.met->v4())[0];
  float pt = (event.muons->at(0).v4() + neutrino).pt();
  return (pt > min_pt) && (pt < max_pt);
}


TTbarGenSemilepSelection::TTbarGenSemilepSelection(uhh2::Context & ctx, std::string ttgen_handlename, float mu_pt_min_, float mu_pt_max_):mu_pt_min(mu_pt_min_),mu_pt_max(mu_pt_max_){
  h_ttbargen = ctx.get_handle<TTbarGen>(ttgen_handlename);
}

bool TTbarGenSemilepSelection::passes(const Event & event){
  TTbarGen ttgen  = event.get(h_ttbargen);
  if(! ttgen.IsSemiLeptonicDecay()) return false;
  GenParticle muon = ttgen.ChargedLepton();
  if(abs(muon.pdgId()) != 13) return false;
  if(muon.pt() < mu_pt_min || muon.pt() > mu_pt_max) return false;
  
  return true;
}
