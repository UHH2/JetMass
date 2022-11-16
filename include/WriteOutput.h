#pragma once
#include "UHH2/core/include/fwd.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/JetMass/include/JetMassUtils.h"

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
  WriteOutput(uhh2::Context & ctx, const std::string & matching_selection_handlename, const std::string & reco_topjet_handlename, const std::string & reco_topjet_chs_handlename, const std::string & gen_topjet_handlename, int output_level = 0);
    virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<MatchingSelection> h_matching_selection;
  double CalculateMJet(vector<PFParticle> Particles);
  vector<double> CalculateMJetVariation(vector<PFParticle>particles, int i, int j, TString cat);
  // vector<double> CalculateMJetResolutionVariation(vector<PFParticle>particles, int i, int j, TString cat);
  std::vector<LorentzVector> variedJets(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi=true);
  std::vector<double> CalculateMJetVariation1(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi=true);
  std::vector<double> CalculatePtVariation(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi=true);
  LorentzVector nominalJet(std::vector<PFParticle> particles, bool apply_puppi=true);

  void ConstructOtherIDs();
  bool inCategory(PFParticle p, TString cat);

  int countBJetsAroundJet(const uhh2::Event & event, TopJet AK8, std::vector<Jet> AK4_jets, JetId btag = DeepJetBTag(DeepJetBTag::WP_MEDIUM), float deltaR_min = 0.8 );
  
  uhh2::Event::Handle<int>h_n_pv;
  uhh2::Event::Handle<int>h_n_trueint_intime;
  uhh2::Event::Handle<int>h_n_trueint_ootimebefore;
  uhh2::Event::Handle<int>h_n_trueint_ootimeafter;
  uhh2::Event::Handle<float>h_n_trueint;

  uhh2::Event::Handle<std::vector<int>> h_trigger_bits;
  std::vector<std::string> trigger_names;

  uhh2::Event::Handle<double>h_msubjets; // mass of summed v4 of subjets
  uhh2::Event::Handle<double>h_mgensubjets; // mass of summed v4 of subjets
  uhh2::Event::Handle<double>h_mgenparticles; // mass of summed v4 of genparticle Candidate from subjets
  uhh2::Event::Handle<double>h_genpt;
  uhh2::Event::Handle<double>h_mjet;
  uhh2::Event::Handle<double>h_mjet_SD;
  uhh2::Event::Handle<double>h_DeepBoostWQCD;
  uhh2::Event::Handle<double>h_DeepBoostZQCD;
  uhh2::Event::Handle<double>h_DeepBoostZbbQCD;
  uhh2::Event::Handle<double>h_MDDeepBoostWQCD;
  uhh2::Event::Handle<double>h_MDDeepBoostZQCD;
  uhh2::Event::Handle<double>h_MDDeepBoostZbbQCD;

  uhh2::Event::Handle<double>h_DeepDoubleBHbbprob;
  uhh2::Event::Handle<double>h_DeepDoubleBQCDprob;
  uhh2::Event::Handle<double>h_MIDeepDoubleBHbbprob;
  uhh2::Event::Handle<double>h_MIDeepDoubleBQCDprob;


  //raw scores of ParticleNet tagger
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbcq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbqq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbc;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbel;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbmu;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probTbta;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probWcq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probWqq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probZbb;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probZcc;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probZqq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probHbb;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probHcc;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probHqqqq;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCDbb;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCDcc;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCDb;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCDc;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCDothers;
  uhh2::Event::Handle<double> h_btag_ParticleNetJetTags_probQCD;

  //binary scores of ParticleNet, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfParticleNetDiscriminatorsJetTags_cfi.py
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_TvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_WvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_ZvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_ZbbvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_HbbvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_HccvsQCD;
  uhh2::Event::Handle<double> h_btag_ParticleNetDiscriminatorsJetTags_H4qvsQCD;

  //raw scores of mass decorrelated ParticleNet tagger
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probXbb;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probXcc;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probXqq;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCDbb;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCDcc;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCDb;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCDc;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCDothers;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetJetTags_probQCD;

  //binary scores of mass decorrelated ParticleNet tagger, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfMassDecorrelatedParticleNetDiscriminatorsJetTags_cfi.py
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XbbvsQCD;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XccvsQCD;
  uhh2::Event::Handle<double> h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XqqvsQCD;

  //mass from ParticleNet mass regression
  uhh2::Event::Handle<double> h_ParticleNetMassRegressionJetTags_mass;



  

  uhh2::Event::Handle<int>h_NextraMBtagDR0p8;
  uhh2::Event::Handle<int>h_NextraTBtagDR0p8;
  uhh2::Event::Handle<int>h_NextraMBtagDR1p0;
  uhh2::Event::Handle<int>h_NextraTBtagDR1p0;


  uhh2::Event::Handle<double>h_pt;
  uhh2::Event::Handle<double>h_eta;
  uhh2::Event::Handle<double>h_phi;
  uhh2::Event::Handle<double>h_mass;

  uhh2::Event::Handle<double>h_pt_1;
  uhh2::Event::Handle<double>h_eta_1;
  uhh2::Event::Handle<double>h_phi_1;
  uhh2::Event::Handle<double>h_mass_1;
  uhh2::Event::Handle<double>h_jecfactor_1;
  
  uhh2::Event::Handle<double>h_pt_AK4,h_genpt_AK4;
  uhh2::Event::Handle<double>h_N2;
  uhh2::Event::Handle<double>h_tau32;
  uhh2::Event::Handle<double>h_tau21;
  uhh2::Event::Handle<double>h_weight;
  uhh2::Event::Handle<double>h_genjetpt;
  uhh2::Event::Handle<double>h_jecfactor;
  uhh2::Event::Handle<double>h_jecfactor_SD;
  uhh2::Event::Handle<double>h_jerfactor_SD_nominal;
  uhh2::Event::Handle<double>h_jerfactor_SD_up;
  uhh2::Event::Handle<double>h_jerfactor_SD_down;
  uhh2::Event::Handle<double>h_jerfactor_SD_JEC_nominal;
  uhh2::Event::Handle<double>h_jerfactor_SD_JEC_up;
  uhh2::Event::Handle<double>h_jerfactor_SD_JEC_down;

  
  uhh2::Event::Handle<double>h_CHF;
  uhh2::Event::Handle<double>h_NHF;
  
  uhh2::Event::Handle<bool>h_IsMergedHiggs;
  uhh2::Event::Handle<bool>h_IsMergedTop;
  uhh2::Event::Handle<bool>h_IsMergedQB;
  uhh2::Event::Handle<bool>h_IsMergedWZ;
  uhh2::Event::Handle<bool>h_IsNotMerged;

  uhh2::Event::Handle<int>h_pdgId_Q1;
  uhh2::Event::Handle<int>h_pdgId_Q2;  

  uhh2::Event::Handle<int>h_n_ak8_reco;
  uhh2::Event::Handle<int>h_n_ak8_gen;  
  

  
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_jetmass_variations;
  TH2F* grid;
  int Nbins_pt, Nbins_eta, Nbins_cat;
  vector<TString> categories;
  vector<int> otherIDs;
  double variation = 0.005;
  std::unique_ptr<MatchingSelection> matching_selection;
  std::unique_ptr<StandaloneTopJetCorrector> softdrop_jec,softdrop_jec_chs;
  bool isMC, do_genStudies, is_WSample, is_ZSample, isTTbarSel, isVJetsSel;

  bool save_jms_jes_study_specific,save_variations,save_selection_variables;
  
  uhh2::Event::Handle<const GenTopJet*> handle_gentopjet;
  uhh2::Event::Handle<const TopJet*> handle_recotopjet,handle_recotopjet_chs;

  uhh2::Event::Handle<double> h_dR_chs_puppi, h_dR_gen_puppi, h_dR_gen_chs;

  uhh2::Event::Handle<double> h_jecfactor_puppi, h_jecfactor_chs;
  uhh2::Event::Handle<double> h_jecfactor_puppi_sd, h_jecfactor_chs_sd;
  
  uhh2::Event::Handle<double> h_pt_gen_ak8, h_mass_gen_ak8, h_msd_gen_ak8;

  uhh2::Event::Handle<double> h_pt_reco_ak8_puppi, h_mass_reco_ak8_puppi, h_msd_reco_ak8_puppi;
  uhh2::Event::Handle<double> h_pt_pf_reco_ak8_puppi,h_pt_pf_sd_reco_ak8_puppi;

  uhh2::Event::Handle<double> h_pt_reco_ak8_chs, h_mass_reco_ak8_chs, h_msd_reco_ak8_chs;
  uhh2::Event::Handle<double> h_pt_pf_reco_ak8_chs,h_pt_pf_sd_reco_ak8_chs;

  uhh2::Event::Handle<int> h_n_particles_puppi, h_n_particles_sd_puppi;
  uhh2::Event::Handle<int> h_n_particles_chs, h_n_particles_sd_chs;

  uhh2::Event::Handle<std::vector<int>> h_PF_multiplicities_reco_ak8_puppi,h_PF_multiplicities_reco_ak8_chs;
  uhh2::Event::Handle<std::vector<int>> h_PF_multiplicities_reco_ak8_puppi_sd,h_PF_multiplicities_reco_ak8_chs_sd;

  uhh2::Event::Handle<std::vector<float>> h_PF_energy_reco_ak8_puppi,h_PF_energy_reco_ak8_chs;
  uhh2::Event::Handle<std::vector<float>> h_PF_energy_reco_ak8_puppi_sd,h_PF_energy_reco_ak8_chs_sd;

  uhh2::Event::Handle<std::vector<float>> h_PF_pt_reco_ak8_puppi,h_PF_pt_reco_ak8_chs;
  uhh2::Event::Handle<std::vector<float>> h_PF_pt_reco_ak8_puppi_sd,h_PF_pt_reco_ak8_chs_sd;
 
  uhh2::Event::Handle<double> h_V_pt;
  

  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_pt_variations_puppi;
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_mass_variations_puppi;
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_msd_variations_puppi;
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_pt_variations_chs;   
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_mass_variations_chs;
  std::vector<std::vector<std::vector< uhh2::Event::Handle<std::vector<double>> >>> h_reco_msd_variations_chs;
  

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
  vector<TString> all_cat = {"undef", "chargedH", "elec", "muon", "gamma", "neutralH", "had", "em", "all"};

};
