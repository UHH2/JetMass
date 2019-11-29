#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/core/include/Event.h"

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
