#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;

ApplyPuppiToPF::ApplyPuppiToPF(){}

bool ApplyPuppiToPF::process(uhh2::Event & event){
  vector<PFParticle>* particles = event.pfparticles;
  for(unsigned int i=0; i<particles->size();i++){
    LorentzVector old_v4 = particles->at(i).v4();
    double factor = particles->at(i).puppiWeight();
    particles->at(i).set_v4(old_v4 * factor);
  }
  return true;
}
