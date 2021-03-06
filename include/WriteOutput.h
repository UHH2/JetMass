#pragma once
#include "UHH2/core/include/fwd.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/JetMass/include/MatchingSelections.h"
#include "TopJetCorrections.h"
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <algorithm>


using namespace std;
using namespace uhh2;


class WriteOutput: public uhh2::AnalysisModule {
public:
    WriteOutput(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & ) override;

private:
  double CalculateMJet(vector<PFParticle> Particles);
  vector<double> CalculateMJetVariation(vector<PFParticle>particles, int i, int j, TString cat);
  void ConstructOtherIDs();
  bool inCategory(PFParticle p, TString cat);

  uhh2::Event::Handle<double>h_msubjets; // mass of summed v4 of subjets
  uhh2::Event::Handle<double>h_mgensubjets; // mass of summed v4 of subjets
  uhh2::Event::Handle<double>h_mgenparticles; // mass of summed v4 of genparticle Candidate from subjets
  uhh2::Event::Handle<double>h_genpt;
  uhh2::Event::Handle<double>h_mjet;
  uhh2::Event::Handle<double>h_mjet_SD;
  uhh2::Event::Handle<double>h_DeepBoost;
  uhh2::Event::Handle<double>h_pt;
  uhh2::Event::Handle<double>h_N2;
  uhh2::Event::Handle<double>h_tau32;
  uhh2::Event::Handle<double>h_tau21;
  uhh2::Event::Handle<double>h_weight;
  uhh2::Event::Handle<double>h_genjetpt;
  uhh2::Event::Handle<double>h_jecfactor;
  uhh2::Event::Handle<double>h_jecfactor_SD;

  uhh2::Event::Handle<bool>h_IsMergedTop;
  uhh2::Event::Handle<bool>h_IsMergedQB;
  uhh2::Event::Handle<bool>h_IsMergedWZ;
  uhh2::Event::Handle<bool>h_IsNotMerged;

  // discriminant variables for old WfromTop selection
  uhh2::Event::Handle<double>h_ht;
  uhh2::Event::Handle<double>h_lepW_pt;
  uhh2::Event::Handle<int>h_nak4;
  uhh2::Event::Handle<int>h_nbtag;
  uhh2::Event::Handle<double>h_deltaPhiAk8Mu;
  
  
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_jetmass_variations;
  TH2F* grid;
  int Nbins_pt, Nbins_eta, Nbins_cat;
  vector<TString> categories;
  vector<int> otherIDs;
  double variation = 0.1;
  std::unique_ptr<MatchingSelection> matching_selection;
  std::unique_ptr<StandaloneTopJetCorrector> softdrop_jec;
  bool isMC, do_genStudies, is_WSample, is_ZSample, isTopSel, isWSel;

  /*
  Particle IDs
  0 - undefined
  1 - charged hadron
  2 - elec
  3 - muon
  4 - gamma
  5 - neutral hadron
  6 - HF tower identified as a hadron
  7 - HF tower identified as an EM particle
  */
  vector<TString> all_cat = {"undef", "chargedH", "elec", "muon", "gamma", "neutralH", "had", "em"};

};
