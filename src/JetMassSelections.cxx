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
