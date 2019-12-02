#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;

ApplyPuppiToPF::ApplyPuppiToPF(){}

bool ApplyPuppiToPF::process(uhh2::Event & event){
  vector<PFParticle>* particles = event.pfparticles;
	for (unsigned int i=0; i<particles->size();i++) {
		//applying Puppi-Weight to LorentzVectorXYZE since the normal LorentzVector(PtEtaPhi) is behaving
		//strange when it is multiplied by 0.
		LorentzVectorXYZE v4XYZ = toXYZ(particles->at(i).v4());
    double puppiWeight = particles->at(i).puppiWeight();
		particles->at(i).set_v4(toPtEtaPhi(v4XYZ * puppiWeight));
	}
  return true;
}
