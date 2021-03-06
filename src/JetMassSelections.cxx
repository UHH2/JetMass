#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarReconstruction.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;


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

METCut::METCut(float minMet_, float maxMet_):
minMet(minMet_), maxMet(maxMet_) {}

bool METCut::passes(const Event & event){

  assert(event.met);

  float MET = event.met->pt();
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

HTCut::HTCut(uhh2::Context & ctx, float min_ht_,float max_ht_):min_ht(min_ht_),max_ht(max_ht_){
  h_ht = ctx.get_handle<double>("HT");
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
