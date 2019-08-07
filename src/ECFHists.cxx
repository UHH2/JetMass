#include "UHH2/JetMass/include/ECFHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

ECFHists::ECFHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  h_comparison = book<TH1F>("comparison", "(N2_jet - N2) / N2", 400, -2, 2);
  h_N2_MIT = book<TH1F>("N2_MIT", "N2", 50, -1.5, 1);
  h_N2_UHH = book<TH1F>("N2_UHH", "N2", 50, -1.5, 1);
  h_LessThan3 = book<TH1F>("LessThan3", "Less than 3 particles in jet", 1, 0, 1);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void ECFHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
  vector<PFParticle>* allparticles = event.pfparticles;

  // get the PFParticles inside the first topjet
  // if mode is soft drop, get particles from subjets
  vector<PFParticle> particles;
  vector<Jet> subjets = topjets->at(0).subjets();
  for(auto subjet: subjets){
    for(const auto candInd : subjet.pfcand_indexs()){
      particles.push_back(allparticles->at(candInd));
    }
  }

  if(particles.size()<3) h_LessThan3->Fill(0.5, weight);

  double N2jet = topjets->at(0).ecfN2_beta1();

  double beta = 1.0;
  bool clear = true;
  vector<fastjet::PseudoJet> particles_pseudojet = ConvertParticles(particles);
  EnergyCorrelations ecf;
  ecf.calcECFN(beta, particles_pseudojet, clear);

  double e2_b1 = ecf.manager->ecfns["2_2"];
  double e3_v2_b1 = ecf.manager->ecfns["3_2"];
  double N2;
  if(e2_b1 != 0) N2 = e3_v2_b1/(e2_b1*e2_b1);
  else           N2 = -0.5;
  h_comparison->Fill((N2jet-N2)/N2, weight);
  h_N2_MIT->Fill(N2, weight);
  h_N2_UHH->Fill(N2jet, weight);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// this applies the factor to PF particles according to the grid
vector<fastjet::PseudoJet> ECFHists::ConvertParticles(vector<PFParticle> oldParticles){
  vector<fastjet::PseudoJet> newParticles;
  for(auto p:oldParticles){
    double px = p.v4().Px();
    double py = p.v4().Py();
    double pz = p.v4().Pz();
    double E = p.v4().E();
    fastjet::PseudoJet newP(px, py, pz, E);
    newParticles.push_back(newP);
  }
  return newParticles;
}


ECFHists::~ECFHists(){}
