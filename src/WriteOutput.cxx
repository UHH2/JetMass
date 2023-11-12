#include "UHH2/core/include/Event.h"
#include "UHH2/JetMass/include/WriteOutput.h"
#include "UHH2/common/include/TTbarReconstruction.h"


using namespace uhh2;
using namespace std;
WriteOutput::WriteOutput(uhh2::Context & ctx, const std::string & matching_selection_handlename, const std::string & reco_topjet_handlename, const std::string & reco_topjet_chs_handlename, const std::string & gen_topjet_handlename, int output_level){

  if(output_level == 0){
    //save everthing
    save_jms_jes_study_specific = true;
    save_variations = true;
    save_selection_variables = true;
  }else if(output_level == -1){
    //do not save variations and things needed for JMS-JES study
    save_jms_jes_study_specific = false;
    save_variations = false;
    save_selection_variables = true;
  }else if(output_level == -2){
    //just save kinematics and variables needed for selections -> use this to create Trees for DDTMapProcessor
    save_jms_jes_study_specific = false;
    save_variations = false;
    save_selection_variables = true;
  }
  
  softdrop_jec.reset(new StandaloneTopJetCorrector(ctx,"AK8PFPuppi"));
  softdrop_jec_chs.reset(new StandaloneTopJetCorrector(ctx,"AK8PFchs"));


  if(save_jms_jes_study_specific){
    h_n_pv = ctx.declare_event_output<int>("n_pv");
    h_n_trueint_intime = ctx.declare_event_output<int>("n_trueint_inttime");
    h_n_trueint_ootimebefore = ctx.declare_event_output<int>("n_trueint_ootimebefore");
    h_n_trueint_ootimeafter = ctx.declare_event_output<int>("n_trueint_ootimeafter");
    h_n_trueint  = ctx.declare_event_output<float>("n_trueint");
  }

  h_pt = ctx.declare_event_output<double>("pt");
  h_eta = ctx.declare_event_output<double>("eta");
  h_phi = ctx.declare_event_output<double>("phi");
  h_mass = ctx.declare_event_output<double>("mass");

  h_pt_1 = ctx.declare_event_output<double>("pt_1");
  h_eta_1 = ctx.declare_event_output<double>("eta_1");
  h_phi_1 = ctx.declare_event_output<double>("phi_1");
  h_mass_1 = ctx.declare_event_output<double>("mass_1");
  h_jecfactor_1 = ctx.declare_event_output<double>("jecfactor_1");
  
  h_N2 = ctx.declare_event_output<double>("N2");
  h_tau32 = ctx.declare_event_output<double>("tau32");
  h_tau21 = ctx.declare_event_output<double>("tau21");

  h_DeepBoostWQCD = ctx.declare_event_output<double>("DeepBoostWQCD");
  h_DeepBoostZQCD = ctx.declare_event_output<double>("DeepBoostZQCD");
  h_DeepBoostZbbQCD = ctx.declare_event_output<double>("DeepBoostZbbQCD");
  h_MDDeepBoostWQCD = ctx.declare_event_output<double>("MDDeepBoostWQCD");
  h_MDDeepBoostZQCD = ctx.declare_event_output<double>("MDDeepBoostZQCD");
  h_MDDeepBoostZbbQCD = ctx.declare_event_output<double>("MDDeepBoostZbbQCD");

  h_DeepDoubleBHbbprob = ctx.declare_event_output<double>("DeepDoubleBHbbprob");
  h_DeepDoubleBQCDprob = ctx.declare_event_output<double>("DeepDoubleBQCDprob");

  h_MIDeepDoubleBHbbprob = ctx.declare_event_output<double>("MIDeepDoubleBHbbprob");
  h_MIDeepDoubleBQCDprob = ctx.declare_event_output<double>("MIDeepDoubleBQCDprob");
  
  //ParticleNet Handles
  
  
  //raw scores of ParticleNet tagger
  h_btag_ParticleNetJetTags_probTbcq = ctx.declare_event_output<double>("ParticleNet_probTbcq");
  h_btag_ParticleNetJetTags_probTbqq = ctx.declare_event_output<double>("ParticleNet_probTbqq");
  h_btag_ParticleNetJetTags_probTbc = ctx.declare_event_output<double>("ParticleNet_probTbc");
  h_btag_ParticleNetJetTags_probTbq = ctx.declare_event_output<double>("ParticleNet_probTbq");
  h_btag_ParticleNetJetTags_probTbel = ctx.declare_event_output<double>("ParticleNet_probTbel");
  h_btag_ParticleNetJetTags_probTbmu = ctx.declare_event_output<double>("ParticleNet_probTbmu");
  h_btag_ParticleNetJetTags_probTbta = ctx.declare_event_output<double>("ParticleNet_probTbta");
  h_btag_ParticleNetJetTags_probWcq = ctx.declare_event_output<double>("ParticleNet_probWcq");
  h_btag_ParticleNetJetTags_probWqq = ctx.declare_event_output<double>("ParticleNet_probWqq");
  h_btag_ParticleNetJetTags_probZbb = ctx.declare_event_output<double>("ParticleNet_probZbb");
  h_btag_ParticleNetJetTags_probZcc = ctx.declare_event_output<double>("ParticleNet_probZcc");
  h_btag_ParticleNetJetTags_probZqq = ctx.declare_event_output<double>("ParticleNet_probZqq");
  h_btag_ParticleNetJetTags_probHbb = ctx.declare_event_output<double>("ParticleNet_probHbb");
  h_btag_ParticleNetJetTags_probHcc = ctx.declare_event_output<double>("ParticleNet_probHcc");
  h_btag_ParticleNetJetTags_probHqqqq = ctx.declare_event_output<double>("ParticleNet_probHqqqq");
  h_btag_ParticleNetJetTags_probQCDbb = ctx.declare_event_output<double>("ParticleNet_probQCDbb");
  h_btag_ParticleNetJetTags_probQCDcc = ctx.declare_event_output<double>("ParticleNet_probQCDcc");
  h_btag_ParticleNetJetTags_probQCDb = ctx.declare_event_output<double>("ParticleNet_probQCDb");
  h_btag_ParticleNetJetTags_probQCDc = ctx.declare_event_output<double>("ParticleNet_probQCDc");
  h_btag_ParticleNetJetTags_probQCDothers = ctx.declare_event_output<double>("ParticleNet_probQCDothers");
  h_btag_ParticleNetJetTags_probQCD = ctx.declare_event_output<double>("ParticleNet_probQCD");

  //binary scores of ParticleNet, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfParticleNetDiscriminatorsJetTags_cfi.py
  h_btag_ParticleNetDiscriminatorsJetTags_TvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_TvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_WvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_WvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_ZvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_ZvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_ZbbvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_ZbbvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_HbbvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_HbbvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_HccvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_HccvsQCD");
  h_btag_ParticleNetDiscriminatorsJetTags_H4qvsQCD = ctx.declare_event_output<double>("ParticleNetDiscriminators_H4qvsQCD");
  h_btag_MassDecorrelatedParticleNetJetTags_probXbb = ctx.declare_event_output<double>("ParticleNetMD_probXbb");
  h_btag_MassDecorrelatedParticleNetJetTags_probXcc = ctx.declare_event_output<double>("ParticleNetMD_probXcc");
  h_btag_MassDecorrelatedParticleNetJetTags_probXqq = ctx.declare_event_output<double>("ParticleNetMD_probXqq");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCDbb = ctx.declare_event_output<double>("ParticleNetMD_probQCDbb");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCDcc = ctx.declare_event_output<double>("ParticleNetMD_probQCDcc");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCDb = ctx.declare_event_output<double>("ParticleNetMD_probQCDb");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCDc = ctx.declare_event_output<double>("ParticleNetMD_probQCDc");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCDothers = ctx.declare_event_output<double>("ParticleNetMD_probQCDothers");
  h_btag_MassDecorrelatedParticleNetJetTags_probQCD = ctx.declare_event_output<double>("ParticleNetMD_probQCD");

  //binary scores of mass decorrelated ParticleNet tagger, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfMassDecorrelatedParticleNetDiscriminatorsJetTags_cfi.py
  h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XbbvsQCD = ctx.declare_event_output<double>("ParticleNetMDDiscriminators_XbbvsQCD");
  h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XccvsQCD = ctx.declare_event_output<double>("ParticleNetMDDiscriminators_XccvsQCD");
  h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XqqvsQCD = ctx.declare_event_output<double>("ParticleNetMDDiscriminators_XqqvsQCD");

  //mass from ParticleNet mass regression
  h_ParticleNetMassRegressionJetTags_mass = ctx.declare_event_output<double>("ParticleNetMassRegression_mass");

  
  h_NextraMBtagDR0p8 = ctx.declare_event_output<int>("NextraMBtagDR0p8");
  h_NextraTBtagDR0p8 = ctx.declare_event_output<int>("NextraTBtagDR0p8");
  h_NextraMBtagDR1p0 = ctx.declare_event_output<int>("NextraMBtagDR1p0");
  h_NextraTBtagDR1p0 = ctx.declare_event_output<int>("NextraTBtagDR1p0");

  
  h_mjet = ctx.declare_event_output<double>("mjet");
  h_mjet_SD = ctx.declare_event_output<double>("mjet_SD");
  h_msubjets = ctx.declare_event_output<double>("msubjets");
  h_mgensubjets = ctx.declare_event_output<double>("mgensubjets");
  h_mgenparticles = ctx.declare_event_output<double>("mgenparticles");
  h_genpt = ctx.declare_event_output<double>("genpt");  
  h_weight = ctx.declare_event_output<double>("weight");
  // h_genjetpt = ctx.declare_event_output<double>("genjetpt");
  h_jecfactor = ctx.declare_event_output<double>("jecfactor");
  h_jecfactor_SD = ctx.declare_event_output<double>("jecfactor_SD");

  h_jerfactor_SD_nominal = ctx.declare_event_output<double>("jerfactor_SD_nominal");
  h_jerfactor_SD_up = ctx.declare_event_output<double>("jerfactor_SD_up");
  h_jerfactor_SD_down = ctx.declare_event_output<double>("jerfactor_SD_down");

  h_jerfactor_SD_JEC_nominal = ctx.declare_event_output<double>("jerfactor_SD_JEC_nominal");
  h_jerfactor_SD_JEC_up = ctx.declare_event_output<double>("jerfactor_SD_JEC_up");
  h_jerfactor_SD_JEC_down = ctx.declare_event_output<double>("jerfactor_SD_JEC_down");


  h_pt_AK4 = ctx.declare_event_output<double>("AK4pt");
  h_genpt_AK4 = ctx.declare_event_output<double>("AK4genpt");
  
  h_CHF = ctx.declare_event_output<double>("CHF");
  h_NHF = ctx.declare_event_output<double>("NHF");
  
  h_IsMergedHiggs = ctx.declare_event_output<bool>("IsMergedHiggs");
  h_IsMergedTop = ctx.declare_event_output<bool>("IsMergedTop");
  h_IsMergedQB = ctx.declare_event_output<bool>("IsMergedQB");
  h_IsMergedWZ = ctx.declare_event_output<bool>("IsMergedWZ");
  h_GenIsMatchedWZ = ctx.declare_event_output<bool>("GenIsMatchedWZ");
  h_IsNotMerged = ctx.declare_event_output<bool>("IsNotMerged");

  h_pdgId_Q1 = ctx.declare_event_output<int>("pdgIdQ1");
  h_pdgId_Q2 = ctx.declare_event_output<int>("pdgIdQ2");
  h_V_pt = ctx.declare_event_output<double>("V_pt");

  h_ps_weights = ctx.declare_event_output<std::vector<float>>("ps_weights");

  h_prefiringweight = ctx.declare_event_output<float>("prefiringweight");
  h_prefiringweight_up = ctx.declare_event_output<float>("prefiringweight_up");
  h_prefiringweight_down = ctx.declare_event_output<float>("prefiringweight_down");
  
  h_trigger_bits = ctx.declare_event_output<std::vector<int>>("trigger_bits");
  trigger_names = {
    "HLT_PFJet320_v*",//0
    "HLT_PFJet400_v*",//1
    "HLT_PFJet450_v*",//2
    "HLT_PFJet500_v*",//3
    "HLT_PFJet550_v*",//4
    "HLT_AK8PFJet320_v*",//5
    "HLT_AK8PFJet400_v*",//6
    "HLT_AK8PFJet450_v*",//7
    "HLT_AK8PFJet500_v*",//8
    "HLT_AK8PFJet550_v*",//9
  };
  
  // h_n_ak8_reco = ctx.declare_event_output<int>("N_ak8_reco");
  // h_n_ak8_gen = ctx.declare_event_output<int>("N_ak8_gen");
  
  // read from xml file
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";

  auto version = ctx.get("dataset_version");
  is_WSample = version.find("WJets") != std::string::npos;
  is_ZSample = version.find("ZJets") != std::string::npos;

  auto selection_ = ctx.get("selection", "");
  isTTbarSel = selection_ == "ttbar";
  isVJetsSel = selection_ == "vjets";

  do_genStudies = string2bool(ctx.get("doGenStudies", "true"));
  
  if(save_variations){
    // read configuration from root file
    TString gfilename = ctx.get("GridFile");
    TFile* gfile = new TFile(locate_file((std::string)gfilename).c_str());
    grid = (TH2F*) gfile->Get("grid");
    TH1F* h_cat = (TH1F*) gfile->Get("categories");
    for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));
    ConstructOtherIDs();
    if(categories.size() > all_cat.size()+1) throw runtime_error("you gave too many categories");
    if(categories.size() + otherIDs.size() -1 != 9) throw runtime_error("categories and size of 'other' does not match");
    
    // get binnings and size of variation
    Nbins_pt =  grid->GetXaxis()->GetNbins();
    Nbins_eta = grid->GetYaxis()->GetNbins();
    Nbins_cat = categories.size();

    // name and create all handles for mjet variations
    if(isMC){
      h_jetmass_variations.resize(Nbins_pt);

      h_reco_pt_variations_puppi.resize(Nbins_pt);
      h_reco_mass_variations_puppi.resize(Nbins_pt);
      h_reco_msd_variations_puppi.resize(Nbins_pt);
      h_reco_pt_variations_chs.resize(Nbins_pt);
      h_reco_mass_variations_chs.resize(Nbins_pt);
      h_reco_msd_variations_chs.resize(Nbins_pt);

      for(int i=0; i<Nbins_pt; i++){
        h_jetmass_variations[i].resize(Nbins_eta);

        h_reco_pt_variations_puppi[i].resize(Nbins_eta);
        h_reco_mass_variations_puppi[i].resize(Nbins_eta);
        h_reco_msd_variations_puppi[i].resize(Nbins_eta);
        h_reco_pt_variations_chs[i].resize(Nbins_eta);
        h_reco_mass_variations_chs[i].resize(Nbins_eta);
        h_reco_msd_variations_chs[i].resize(Nbins_eta);

        for(int j=0; j<Nbins_eta; j++){
          h_jetmass_variations[i][j].resize(Nbins_cat);

          h_reco_pt_variations_puppi[i][j].resize(Nbins_cat);
          h_reco_mass_variations_puppi[i][j].resize(Nbins_cat);
          h_reco_msd_variations_puppi[i][j].resize(Nbins_cat);
          h_reco_pt_variations_chs[i][j].resize(Nbins_cat);
          h_reco_mass_variations_chs[i][j].resize(Nbins_cat);
          h_reco_msd_variations_chs[i][j].resize(Nbins_cat);

          for(int k=0; k<Nbins_cat; k++){
            TString bin_name = "_" + to_string(i) + "_" + to_string(j) + "_" + categories[k];
            TString handlename = "mjet"+bin_name;
            h_jetmass_variations[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)handlename);

            // pt_reco_ak8_chs
            h_reco_pt_variations_puppi[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("pt_reco_ak8_puppi"+bin_name));
            h_reco_mass_variations_puppi[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("mass_reco_ak8_puppi"+bin_name));
            h_reco_msd_variations_puppi[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("msd_reco_ak8_puppi"+bin_name));
            h_reco_pt_variations_chs[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("pt_reco_ak8_chs"+bin_name));
            h_reco_mass_variations_chs[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("mass_reco_ak8_chs"+bin_name));
            h_reco_msd_variations_chs[i][j][k] = ctx.declare_event_output<std::vector<double>>((string)("msd_reco_ak8_chs"+bin_name));

          }
        }
      }
    }
  }
  h_matching_selection = ctx.get_handle<MatchingSelection>(matching_selection_handlename);

  
  if(save_jms_jes_study_specific){
    // Stuff for pt, msd response studies
    handle_gentopjet = ctx.get_handle<const GenTopJet*>(gen_topjet_handlename);
    handle_recotopjet = ctx.get_handle<const TopJet*>(reco_topjet_handlename);
    handle_recotopjet_chs = ctx.get_handle<const TopJet*>(reco_topjet_chs_handlename);

    h_mass_reco_ak8_puppi = ctx.declare_event_output<double>("mass_reco_ak8_puppi");
    h_msd_reco_ak8_puppi = ctx.declare_event_output<double>("msd_reco_ak8_puppi");
    h_pt_reco_ak8_puppi = ctx.declare_event_output<double>("pt_reco_ak8_puppi");
    h_pt_pf_reco_ak8_puppi = ctx.declare_event_output<double>("pt_pf_reco_ak8_puppi");
    h_pt_pf_sd_reco_ak8_puppi = ctx.declare_event_output<double>("pt_pf_sd_reco_ak8_puppi");

    h_mass_reco_ak8_chs = ctx.declare_event_output<double>("mass_reco_ak8_chs");
    h_msd_reco_ak8_chs = ctx.declare_event_output<double>("msd_reco_ak8_chs");
    h_pt_reco_ak8_chs = ctx.declare_event_output<double>("pt_reco_ak8_chs");
    h_pt_pf_reco_ak8_chs = ctx.declare_event_output<double>("pt_pf_reco_ak8_chs");
    h_pt_pf_sd_reco_ak8_chs = ctx.declare_event_output<double>("pt_pf_sd_reco_ak8_chs");

    h_pt_gen_ak8 = ctx.declare_event_output<double>("pt_gen_ak8");
    h_mass_gen_ak8 = ctx.declare_event_output<double>("mass_gen_ak8");
    h_msd_gen_ak8 = ctx.declare_event_output<double>("msd_gen_ak8"); 
    h_n2_beta1_gen = ctx.declare_event_output<double>("n2_beta1_gen");
    h_n2_beta2_gen = ctx.declare_event_output<double>("n2_beta2_gen");

    // h_gen_pt_0 = ctx.declare_event_output<double>("gen_pt_0");
    // h_gen_eta_0 = ctx.declare_event_output<double>("gen_eta_0");
    // h_gen_phi_0 = ctx.declare_event_output<double>("gen_phi_0");
    // h_gen_pt_0 = ctx.declare_event_output<double>("gen_pt_0");
    
    
    h_dR_chs_puppi = ctx.declare_event_output<double>("dr_chs_puppi");
    h_dR_gen_puppi = ctx.declare_event_output<double>("dr_gen_puppi");
    h_dR_gen_chs = ctx.declare_event_output<double>("dr_gen_chs");

    h_jecfactor_puppi = ctx.declare_event_output<double>("jecfactor_puppi");
    h_jecfactor_puppi_sd = ctx.declare_event_output<double>("jecfactor_puppi_sd");
    h_jecfactor_chs = ctx.declare_event_output<double>("jecfactor_chs");
    h_jecfactor_chs_sd = ctx.declare_event_output<double>("jecfactor_chs_sd");

    h_n_particles_puppi = ctx.declare_event_output<int>("n_particles_puppi");
    h_n_particles_sd_puppi = ctx.declare_event_output<int>("n_particles_sd_puppi");

    h_n_particles_chs = ctx.declare_event_output<int>("n_particles_chs");
    h_n_particles_sd_chs = ctx.declare_event_output<int>("n_particles_sd_chs");


    h_PF_multiplicities_reco_ak8_puppi = ctx.declare_event_output<std::vector<int>>("PF_multiplicity_reco_ak8_puppi");
    h_PF_multiplicities_reco_ak8_chs = ctx.declare_event_output<std::vector<int>>("PF_multiplicity_reco_ak8_chs");
    h_PF_multiplicities_reco_ak8_puppi_sd = ctx.declare_event_output<std::vector<int>>("PF_multiplicity_reco_ak8_puppi_sd");
    h_PF_multiplicities_reco_ak8_chs_sd = ctx.declare_event_output<std::vector<int>>("PF_multiplicity_reco_ak8_chs_sd");
  
    h_PF_energy_reco_ak8_puppi = ctx.declare_event_output<std::vector<float>>("PF_energy_reco_ak8_puppi");
    h_PF_energy_reco_ak8_chs = ctx.declare_event_output<std::vector<float>>("PF_energy_reco_ak8_chs");
    h_PF_energy_reco_ak8_puppi_sd = ctx.declare_event_output<std::vector<float>>("PF_energy_reco_ak8_puppi_sd");
    h_PF_energy_reco_ak8_chs_sd = ctx.declare_event_output<std::vector<float>>("PF_energy_reco_ak8_chs_sd");
  
    h_PF_pt_reco_ak8_puppi = ctx.declare_event_output<std::vector<float>>("PF_pt_reco_ak8_puppi");
    h_PF_pt_reco_ak8_chs = ctx.declare_event_output<std::vector<float>>("PF_pt_reco_ak8_chs");
    h_PF_pt_reco_ak8_puppi_sd = ctx.declare_event_output<std::vector<float>>("PF_pt_reco_ak8_puppi_sd");
    h_PF_pt_reco_ak8_chs_sd = ctx.declare_event_output<std::vector<float>>("PF_pt_reco_ak8_chs_sd");
  }
  



  
}

bool WriteOutput::process(uhh2::Event & event){
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1 && event.gentopjets->size() < 1) return false;

  //PartonShower weights
  std::vector<float> ps_weights = {};
  if(isMC){
    ps_weights = event.genInfo->weights();
  }
  event.set(h_ps_weights, ps_weights);

  event.set(h_prefiringweight, event.prefiringWeight);
  event.set(h_prefiringweight_up, event.prefiringWeightUp);
  event.set(h_prefiringweight_down, event.prefiringWeightDown);
  
  std::vector<int> trigger_results = {};
  for(auto trigger_name: trigger_names){
    auto trigger_index = event.get_trigger_index(trigger_name);
    // check if trigger actually exists for this run (write -1 if not)
    if ( event.lookup_trigger_index( trigger_index ) ) {
      trigger_results.push_back( event.passes_trigger( trigger_index ) );
    }else{
      trigger_results.push_back(-1);
    }
  }
  event.set(h_trigger_bits,trigger_results);
  
  if(save_jms_jes_study_specific){
    float n_pv = event.pvs->size();
    float n_trueint_intime = -999.0;
    float n_trueint_ootimebefore = -999.0;
    float n_trueint_ootimeafter = -999.0;
    float n_trueint =  -999.0;
    
    if(isMC){
      n_trueint_intime = event.genInfo->pileup_NumInteractions_intime();
      n_trueint_ootimebefore = event.genInfo->pileup_NumInteractions_ootbefore();
      n_trueint_ootimeafter = event.genInfo->pileup_NumInteractions_ootafter();
      n_trueint =  event.genInfo->pileup_TrueNumInteractions();
    }

    event.set(h_n_pv, n_pv);
    event.set(h_n_trueint_intime, n_trueint_intime);
    event.set(h_n_trueint_ootimebefore, n_trueint_ootimebefore);
    event.set(h_n_trueint_ootimeafter, n_trueint_ootimeafter);
    event.set(h_n_trueint, n_trueint);
  }
  
  MatchingSelection matching_selection = event.get(h_matching_selection);

  event.set(h_pdgId_Q1, matching_selection.FlavourQ1() );
  event.set(h_pdgId_Q2, matching_selection.FlavourQ2() );
  const GenParticle *genV = matching_selection.get_genV();
  if(genV)event.set(h_V_pt, genV->pt());
  
  LorentzVector softdrop_jet_raw;
  vector<PFParticle>* allparticles = event.pfparticles;
  vector<PFParticle> particles;
  if(event.topjets->size()>0){
    vector<Jet> subjets = topjets->at(0).subjets();  
    // find all pf particles inside the subjets
    for(auto subjet: subjets){
      softdrop_jet_raw += subjet.v4()*subjet.JEC_factor_raw();
      for(const auto candInd : subjet.pfcand_indexs()){
        particles.push_back(allparticles->at(candInd));
      }
    }
  }  
  TopJet candidateJet;
  int n_topjets = event.topjets->size();
  if(event.topjets->size()>0) candidateJet = topjets->at(0);
  double pt = n_topjets>0 ? candidateJet.v4().Pt()*candidateJet.JEC_factor_raw() : -1.0;
  double eta = n_topjets>0 ? candidateJet.v4().Eta() : -1.0;
  double phi = n_topjets>0 ? candidateJet.v4().Phi() : -1.0;
  double mass = n_topjets>0 ? candidateJet.v4().M()*candidateJet.JEC_factor_raw() : -1.0;
  
  double N2 = n_topjets>0 ? candidateJet.ecfN2_beta1() : -1.0;
  double mjet = n_topjets>0 ? CalculateMJet(particles) : -1.0;
  double jecfactor = n_topjets>0 ? 1.0/candidateJet.JEC_factor_raw() : -1.0;
  double jecfactor_SD = n_topjets>0 ? softdrop_jec->getJecFactor(event,softdrop_jet_raw) : -1.0;
  double mjet_SD = n_topjets>0 ? candidateJet.softdropmass() : -1.0;


  int n_gentopjets = event.gentopjets->size();
  double tau32 = 0;
  if(candidateJet.tau2() > 0) tau32 = n_topjets>0 ? candidateJet.tau3()/candidateJet.tau2() : -1.0;
  double tau21 = 0;
  if(candidateJet.tau1() > 0) tau21 = n_topjets>0 ? candidateJet.tau2()/candidateJet.tau1() : -1.0;

  double AK8CHF = n_topjets>0 ? candidateJet.chargedHadronEnergyFraction() : -1.0;
  double AK8NHF = n_topjets>0 ? candidateJet.neutralHadronEnergyFraction() : -1.0;
  
  // set variations for MC
  if(isMC & save_variations){
    for(int i=0; i<Nbins_pt; i++){
      for(int j=0; j<Nbins_eta; j++){
        for(int k=0; k<Nbins_cat; k++){
          int ptbin = i+1;
          int etabin = j+1;
          vector<double> mjet_var = CalculateMJetVariation(particles, ptbin, etabin, categories[k]);
          event.set(h_jetmass_variations[i][j][k], mjet_var);
        }
      }
    }
  }
  bool IsMergedHiggs = n_topjets>0 ? matching_selection.passes_matching(candidateJet,MatchingSelection::oIsMergedHiggs) : false;
  bool IsMergedTop = n_topjets>0 ? matching_selection.passes_matching(candidateJet,MatchingSelection::oIsMergedTop) : false;
  bool IsMergedQB  = n_topjets>0 ? matching_selection.passes_matching(candidateJet,MatchingSelection::oIsMergedQB) : false;
  bool IsMergedWZ  = n_topjets>0 ? matching_selection.passes_matching(candidateJet,MatchingSelection::oIsMergedV) : false;
  bool IsNotMerged = n_topjets>0 ? matching_selection.passes_matching(candidateJet,MatchingSelection::oIsNotMerged) : false;


  float genpt,m_genparticles,m_gensubjets;

  if(isMC && do_genStudies){    
    if(event.gentopjets->size()<1)return false;
    const GenTopJet * matched_gentopjet = closestParticle(candidateJet, *event.gentopjets);    
    vector<GenJet> gensubjets = matched_gentopjet->subjets();
    LorentzVector softdrop_genjet;
    for(auto gensubjet: gensubjets){
      softdrop_genjet += gensubjet.v4();
    }
    m_gensubjets = softdrop_genjet.M();    
    genpt =  matched_gentopjet->pt();
  }else{
    genpt = -1.;
    m_gensubjets = -1.;
    m_genparticles = -1.;
  }
  
  event.set(h_pt, pt);//v4().pt() * jec_factor_raw
  event.set(h_eta, eta);
  event.set(h_phi, phi);
  event.set(h_mass, mass); //v4().M() * jec_factor_raw
  event.set(h_mjet, mjet); //sum(particles).v4().M() * jec_factor_raw -> should be identical to softdropmass
  
  event.set(h_msubjets, softdrop_jet_raw.M());
  event.set(h_mjet_SD,mjet_SD); // softdropmass from userfloat taken from MINIAOD, not 100% identical to mjet, probably due to precisions in MINIAOD userfloats? 

  float pt_1(-999.9),eta_1(-999.9),phi_1(-999.9),mass_1(-999.9),jecfactor_1(-999.9);
  if(event.topjets->size()>1){
    pt_1 = event.topjets->at(1).v4().Pt()*event.topjets->at(1).JEC_factor_raw();
    eta_1 = event.topjets->at(1).v4().Eta();
    phi_1 = event.topjets->at(1).v4().Phi();
    mass_1 = event.topjets->at(1).v4().M()*event.topjets->at(1).JEC_factor_raw();
    jecfactor_1 = 1.0/event.topjets->at(1).JEC_factor_raw();  
  }
  event.set(h_pt_1, pt_1);//v4().pt() * jec_factor_raw
  event.set(h_eta_1, eta_1);
  event.set(h_phi_1, phi_1);
  event.set(h_mass_1, mass_1); //v4().M() * jec_factor_raw
  event.set(h_jecfactor_1, jecfactor_1);
  
  event.set(h_mgensubjets, m_gensubjets);
  event.set(h_mgenparticles, m_genparticles);
  event.set(h_genpt, genpt);

  
  event.set(h_N2, N2);
  event.set(h_tau32, tau32);
  event.set(h_tau21, tau21);
  event.set(h_DeepBoostWQCD, n_topjets>0 ? candidateJet.btag_DeepBoosted_WvsQCD() : -1.0);
  event.set(h_DeepBoostZQCD, n_topjets>0 ? candidateJet.btag_DeepBoosted_ZvsQCD() : -1.0);
  event.set(h_DeepBoostZbbQCD, n_topjets>0 ? candidateJet.btag_DeepBoosted_ZbbvsQCD() : -1.0);
  event.set(h_MDDeepBoostWQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedDeepBoosted_WvsQCD() : -1.0);
  event.set(h_MDDeepBoostZQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedDeepBoosted_ZvsQCD() : -1.0);
  event.set(h_MDDeepBoostZbbQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedDeepBoosted_ZbbvsQCD() : -1.0);

  event.set(h_DeepDoubleBHbbprob, n_topjets>0 ? candidateJet.btag_DeepDoubleBvLJet_probHbb() : -1.0);
  event.set(h_DeepDoubleBQCDprob, n_topjets>0 ? candidateJet.btag_DeepDoubleBvLJet_probQCD() : -1.0);
  event.set(h_MIDeepDoubleBHbbprob, n_topjets>0 ? candidateJet.btag_MassIndependentDeepDoubleBvLJet_probHbb() : -1.0);
  event.set(h_MIDeepDoubleBQCDprob, n_topjets>0 ? candidateJet.btag_MassIndependentDeepDoubleBvLJet_probQCD() : -1.0);

    //raw scores of ParticleNet tagger
  event.set(h_btag_ParticleNetJetTags_probTbcq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbcq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbqq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbqq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbc, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbc() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbel, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbel() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbmu, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbmu() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probTbta, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probTbta() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probWcq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probWcq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probWqq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probWqq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probZbb, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probZbb() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probZcc, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probZcc() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probZqq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probZqq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probHbb, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probHbb() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probHcc, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probHcc() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probHqqqq, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probHqqqq() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCDbb, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCDbb() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCDcc, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCDcc() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCDb, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCDb() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCDc, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCDc() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCDothers, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCDothers() : -1.0);
  event.set(h_btag_ParticleNetJetTags_probQCD, n_topjets>0 ? candidateJet.btag_ParticleNetJetTags_probQCD() : -1.0);
  
  //binary scores of ParticleNet, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfParticleNetDiscriminatorsJetTags_cfi.py
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_TvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_TvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_WvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_WvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_ZvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_ZvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_ZbbvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_ZbbvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_HbbvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_HbbvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_HccvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_HccvsQCD() : -1.0);
  event.set(h_btag_ParticleNetDiscriminatorsJetTags_H4qvsQCD, n_topjets>0 ? candidateJet.btag_ParticleNetDiscriminatorsJetTags_H4qvsQCD() : -1.0);
  
  //raw scores of mass decorrelated ParticleNet tagger
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probXbb, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probXbb() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probXcc, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probXcc() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probXqq, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probXqq() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCDbb, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCDbb() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCDcc, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCDcc() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCDb, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCDb() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCDc, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCDc() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCDothers, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCDothers() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetJetTags_probQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetJetTags_probQCD() : -1.0);

  //binary scores of mass decorrelated ParticleNet tagger, see https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/python/pfMassDecorrelatedParticleNetDiscriminatorsJetTags_cfi.py
  event.set(h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XbbvsQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XbbvsQCD() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XccvsQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XccvsQCD() : -1.0);
  event.set(h_btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XqqvsQCD, n_topjets>0 ? candidateJet.btag_MassDecorrelatedParticleNetDiscriminatorsJetTags_XqqvsQCD() : -1.0);

  //mass from ParticleNet mass regression
  event.set(h_ParticleNetMassRegressionJetTags_mass, n_topjets>0 ? candidateJet.ParticleNetMassRegressionJetTags_mass() : -1.0);

  


  
  event.set(h_NextraMBtagDR0p8, n_topjets>0 ? countBJetsAroundJet(event, candidateJet, *event.jets, DeepJetBTag(DeepJetBTag::WP_MEDIUM),0.8):-1.0);
  event.set(h_NextraTBtagDR0p8, n_topjets>0 ? countBJetsAroundJet(event, candidateJet, *event.jets, DeepJetBTag(DeepJetBTag::WP_TIGHT),0.8):-1.0);
  event.set(h_NextraMBtagDR1p0, n_topjets>0 ? countBJetsAroundJet(event, candidateJet, *event.jets, DeepJetBTag(DeepJetBTag::WP_MEDIUM),1.0):-1.0);
  event.set(h_NextraTBtagDR1p0, n_topjets>0 ? countBJetsAroundJet(event, candidateJet, *event.jets, DeepJetBTag(DeepJetBTag::WP_TIGHT),1.0):-1.0);

  //std::cout << "weight stored in tree: " << event.weight << std::endl;
  event.set(h_weight, event.weight);
  // event.set(h_genjetpt, genjetpt);
  event.set(h_jecfactor, jecfactor);
  event.set(h_jecfactor_SD, jecfactor_SD);

  event.set(h_jerfactor_SD_nominal, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,0,1.0): -1.0);
  event.set(h_jerfactor_SD_up, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,1,1.0):-1.0);
  event.set(h_jerfactor_SD_down, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,-1,1.0):-1.0);

  event.set(h_jerfactor_SD_JEC_nominal, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,0,jecfactor_SD):-1.0);
  event.set(h_jerfactor_SD_JEC_up, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,1,jecfactor_SD):-1.0);
  event.set(h_jerfactor_SD_JEC_down, n_topjets>0 ? softdrop_jec->getJERSmearingFactor(event,softdrop_jet_raw,-1,jecfactor_SD):-1.0);
  
    
  double ak4_pt = -1.0;
  double ak4_genpt = -1.0;
  if(event.jets->size()>0){
    ak4_pt = event.jets->at(0).pt();
    if(isMC && event.genjets && event.genjets->size() >0){
      const GenJet * closest_genjet = closestParticle(event.jets->at(0), *event.genjets);
      ak4_genpt = closest_genjet->pt();
    }
      // cout << "test"<<endl; ak4_genpt = event.genjets->at(event.jets->at(0).genjet_index()).pt();
  }
  event.set(h_pt_AK4,ak4_pt);
  event.set(h_genpt_AK4,ak4_genpt);
  
  event.set(h_CHF,AK8CHF);
  event.set(h_NHF,AK8NHF);
  
  event.set(h_IsMergedHiggs, IsMergedHiggs);
  event.set(h_IsMergedTop, IsMergedTop);
  event.set(h_IsMergedQB, IsMergedQB); 
  event.set(h_IsMergedWZ, IsMergedWZ);   
  event.set(h_IsNotMerged, IsNotMerged);

  // CHS Puppi comparison stuff
  if(save_jms_jes_study_specific){
    const GenTopJet* gen_jet(NULL);
    if(event.is_valid(handle_gentopjet)) gen_jet = event.get(handle_gentopjet);
    const TopJet* reco_puppi_jet(NULL);
    if(event.is_valid(handle_recotopjet)) reco_puppi_jet = event.get(handle_recotopjet);
    const TopJet* reco_chs_jet(NULL);
    if(event.is_valid(handle_recotopjet_chs)) reco_chs_jet = event.get(handle_recotopjet_chs);

  
    std::vector<PFParticle> particles_puppi_sd,particles_puppi;
    std::vector<PFParticle> particles_chs_sd,particles_chs;
    int n_particles_puppi(-1.0),n_particles_sd_puppi(-1.0);
    int n_particles_chs(-1.0),n_particles_sd_chs(-1.0);

  
    float pt_gen_ak8(-1.0),mass_gen_ak8(-1.0),msd_gen_ak8(-1.0);
    float n2_beta1_gen(-1.0),n2_beta2_gen(-1.0);
    bool GenIsMatchedWZ(false);
    if(gen_jet){
      pt_gen_ak8 = gen_jet->pt();
      mass_gen_ak8 = gen_jet->v4().M();
      msd_gen_ak8 = gen_jet->softdropmass();
      n2_beta1_gen = gen_jet->ecfN2_beta1();
      n2_beta2_gen = gen_jet->ecfN2_beta2();
      GenIsMatchedWZ = matching_selection.passes_matching(*gen_jet,MatchingSelection::oIsMergedV);
    }
    event.set(h_GenIsMatchedWZ,GenIsMatchedWZ);
    // std::cout << "mass_gen: \t" << mass_gen_ak8 << std::endl;
    // std::cout << "msd_gen: \t" << msd_gen_ak8 << std::endl;
    event.set(h_pt_gen_ak8,pt_gen_ak8);
    event.set(h_mass_gen_ak8,mass_gen_ak8);
    event.set(h_msd_gen_ak8,msd_gen_ak8);
    event.set(h_n2_beta1_gen, n2_beta1_gen);
    event.set(h_n2_beta2_gen, n2_beta2_gen);

    float pt_reco_ak8_puppi(-1.0),mass_reco_ak8_puppi(-1.0),msd_reco_ak8_puppi(-1.0);
    float pt_reco_ak8_puppi_pf(-1.0),pt_reco_ak8_puppi_pf_sd(-1.0);
    float jecfactor_puppi(0.),jecfactor_puppi_sd(0.);
    if(reco_puppi_jet){
      pt_reco_ak8_puppi = reco_puppi_jet->pt();
      // mass_reco_ak8_puppi = reco_puppi_jet->v4().M();
      // msd_reco_ak8_puppi = reco_puppi_jet->softdropmass();

      for(const auto candInd : reco_puppi_jet->pfcand_indexs()){
        particles_puppi.push_back(allparticles->at(candInd));
      }
      n_particles_puppi = particles_puppi.size();
    
      for(auto subjet: reco_puppi_jet->subjets()){
        for(const auto candInd : subjet.pfcand_indexs()){
          particles_puppi_sd.push_back(allparticles->at(candInd));
        }
      }
      n_particles_sd_puppi = particles_puppi_sd.size();

      LorentzVector nominal_puppi_jet = nominalJet(particles_puppi);
      LorentzVector nominal_puppi_jet_sd = nominalJet(particles_puppi_sd);

      jecfactor_puppi = softdrop_jec->getJecFactor(event,nominal_puppi_jet);
      jecfactor_puppi_sd = softdrop_jec->getJecFactor(event,nominal_puppi_jet_sd);

      mass_reco_ak8_puppi = nominal_puppi_jet.M();
      msd_reco_ak8_puppi = nominal_puppi_jet_sd.M();  
      pt_reco_ak8_puppi_pf = nominal_puppi_jet.pt();
      pt_reco_ak8_puppi_pf_sd = nominal_puppi_jet_sd.pt();
  
    }
    event.set(h_pt_reco_ak8_puppi,pt_reco_ak8_puppi);

    event.set(h_pt_pf_reco_ak8_puppi,pt_reco_ak8_puppi_pf);
    event.set(h_pt_pf_sd_reco_ak8_puppi,pt_reco_ak8_puppi_pf_sd);
  
    event.set(h_mass_reco_ak8_puppi,mass_reco_ak8_puppi);
    event.set(h_msd_reco_ak8_puppi,msd_reco_ak8_puppi);

    event.set(h_n_particles_puppi, n_particles_puppi);
    event.set(h_n_particles_sd_puppi, n_particles_sd_puppi);

    event.set(h_jecfactor_puppi, jecfactor_puppi);
    event.set(h_jecfactor_puppi_sd, jecfactor_puppi_sd);

  
    float pt_reco_ak8_chs(-1.0),mass_reco_ak8_chs(-1.0),msd_reco_ak8_chs(-1.0);
    float pt_reco_ak8_chs_pf(-1.0),pt_reco_ak8_chs_pf_sd(-1.0);
    float jecfactor_chs(0.),jecfactor_chs_sd(0.);

  
    if(reco_chs_jet){
      pt_reco_ak8_chs = reco_chs_jet->pt();
      // mass_reco_ak8_chs = reco_chs_jet->v4().M();
      // msd_reco_ak8_chs = reco_chs_jet->softdropmass();    

      for(const auto candInd : reco_chs_jet->pfcand_indexs()){
        particles_chs.push_back(allparticles->at(candInd));
      }
      n_particles_chs = particles_chs.size();
    
      for(auto subjet: reco_chs_jet->subjets()){
        for(const auto candInd : subjet.pfcand_indexs()){
          particles_chs_sd.push_back(allparticles->at(candInd));
        }
      }
      n_particles_sd_chs = particles_chs_sd.size();

      LorentzVector nominal_chs_jet = nominalJet(particles_chs,false);
      LorentzVector nominal_chs_jet_sd = nominalJet(particles_chs_sd,false);
    
      jecfactor_chs = softdrop_jec_chs->getJecFactor(event,nominal_chs_jet);
      jecfactor_chs_sd = softdrop_jec_chs->getJecFactor(event,nominal_chs_jet_sd);

      mass_reco_ak8_chs = nominal_chs_jet.M();
      msd_reco_ak8_chs = nominal_chs_jet_sd.M();  
      pt_reco_ak8_chs_pf = nominal_chs_jet.pt();
      pt_reco_ak8_chs_pf_sd = nominal_chs_jet_sd.pt();

    }
    event.set(h_pt_reco_ak8_chs,pt_reco_ak8_chs);

    event.set(h_pt_pf_reco_ak8_chs,pt_reco_ak8_chs_pf);
    event.set(h_pt_pf_sd_reco_ak8_chs,pt_reco_ak8_chs_pf_sd);
  
    event.set(h_mass_reco_ak8_chs,mass_reco_ak8_chs);
    event.set(h_msd_reco_ak8_chs,msd_reco_ak8_chs);

    event.set(h_n_particles_chs, n_particles_chs);
    event.set(h_n_particles_sd_chs, n_particles_sd_chs);
  
    event.set(h_jecfactor_chs, jecfactor_chs);
    event.set(h_jecfactor_chs_sd, jecfactor_chs_sd);

  
    // std::cout << "puppi particles summary:" << std::endl;
    // std::cout << "chargedHadrons: " << PFMultiplicity(particles_puppi, PFParticle::EParticleID::eH) << std::endl;
    // std::cout << "jet chargedHadrons multi: " << reco_puppi_jet->chargedMultiplicity() << std::endl;
    std::vector<int> PF_multiplicities_reco_ak8_puppi,PF_multiplicities_reco_ak8_chs,PF_multiplicities_reco_ak8_puppi_sd,PF_multiplicities_reco_ak8_chs_sd;
    std::vector<float> PF_energy_reco_ak8_puppi,PF_energy_reco_ak8_chs,PF_energy_reco_ak8_puppi_sd,PF_energy_reco_ak8_chs_sd;
    std::vector<float> PF_pt_reco_ak8_puppi,PF_pt_reco_ak8_chs,PF_pt_reco_ak8_puppi_sd,PF_pt_reco_ak8_chs_sd;
    
    for(unsigned int i=0;i<all_cat.size();i++){
      PF_multiplicities_reco_ak8_puppi.push_back(PFMultiplicity(particles_puppi, i));
      PF_multiplicities_reco_ak8_chs.push_back(PFMultiplicity(particles_chs, i));
      PF_multiplicities_reco_ak8_puppi_sd.push_back(PFMultiplicity(particles_puppi_sd, i));
      PF_multiplicities_reco_ak8_chs_sd.push_back(PFMultiplicity(particles_chs_sd, i));
    
      PF_energy_reco_ak8_puppi.push_back(PFEnergy(particles_puppi, i));
      PF_energy_reco_ak8_chs.push_back(PFEnergy(particles_chs, i));
      PF_energy_reco_ak8_puppi_sd.push_back(PFEnergy(particles_puppi_sd, i));
      PF_energy_reco_ak8_chs_sd.push_back(PFEnergy(particles_chs_sd, i));
    
      PF_pt_reco_ak8_puppi.push_back(PFTransverseMomentum(particles_puppi, i));
      PF_pt_reco_ak8_chs.push_back(PFTransverseMomentum(particles_chs, i));
      PF_pt_reco_ak8_puppi_sd.push_back(PFTransverseMomentum(particles_puppi_sd, i));
      PF_pt_reco_ak8_chs_sd.push_back(PFTransverseMomentum(particles_chs_sd, i));
    }
    event.set(h_PF_multiplicities_reco_ak8_puppi,PF_multiplicities_reco_ak8_puppi);
    event.set(h_PF_multiplicities_reco_ak8_chs,PF_multiplicities_reco_ak8_chs);
    event.set(h_PF_multiplicities_reco_ak8_puppi_sd,PF_multiplicities_reco_ak8_puppi_sd);
    event.set(h_PF_multiplicities_reco_ak8_chs_sd,PF_multiplicities_reco_ak8_chs_sd);
  
    event.set(h_PF_energy_reco_ak8_puppi,PF_energy_reco_ak8_puppi);
    event.set(h_PF_energy_reco_ak8_chs,PF_energy_reco_ak8_chs);
    event.set(h_PF_energy_reco_ak8_puppi_sd,PF_energy_reco_ak8_puppi_sd);
    event.set(h_PF_energy_reco_ak8_chs_sd,PF_energy_reco_ak8_chs_sd);
  
    event.set(h_PF_pt_reco_ak8_puppi,PF_pt_reco_ak8_puppi);
    event.set(h_PF_pt_reco_ak8_chs,PF_pt_reco_ak8_chs);
    event.set(h_PF_pt_reco_ak8_puppi_sd,PF_pt_reco_ak8_puppi_sd);
    event.set(h_PF_pt_reco_ak8_chs_sd,PF_pt_reco_ak8_chs_sd);

  

    // set variations for MC
    if(isMC){
      for(int i=0; i<Nbins_pt; i++){
        for(int j=0; j<Nbins_eta; j++){
          for(int k=0; k<Nbins_cat; k++){
            int ptbin = i+1;
            int etabin = j+1;
            vector<double> puppi_mass_var = CalculateMJetVariation1(particles_puppi, ptbin, etabin, categories[k]);
            event.set(h_reco_mass_variations_puppi[i][j][k], puppi_mass_var);
            vector<double> puppi_msd_var = CalculateMJetVariation1(particles_puppi_sd, ptbin, etabin, categories[k]);
            event.set(h_reco_msd_variations_puppi[i][j][k], puppi_msd_var);
            vector<double> puppi_pt_var = CalculatePtVariation(particles_puppi, ptbin, etabin, categories[k]);
            event.set(h_reco_pt_variations_puppi[i][j][k], puppi_pt_var);

            vector<double> chs_mass_var = CalculateMJetVariation1(particles_chs, ptbin, etabin, categories[k],false);
            event.set(h_reco_mass_variations_chs[i][j][k], chs_mass_var);
            vector<double> chs_msd_var = CalculateMJetVariation1(particles_chs_sd, ptbin, etabin, categories[k],false);
            event.set(h_reco_msd_variations_chs[i][j][k], chs_msd_var);
            vector<double> chs_pt_var = CalculatePtVariation(particles_chs, ptbin, etabin, categories[k],false);
            event.set(h_reco_pt_variations_chs[i][j][k], chs_pt_var);
          }
        }
      }
    }

    double dr_chs_puppi(-999.);
    if(reco_chs_jet && reco_puppi_jet) dr_chs_puppi = deltaR(*reco_chs_jet,*reco_puppi_jet);
    event.set(h_dR_chs_puppi, dr_chs_puppi);

    double dr_gen_puppi(-999.);
    if(gen_jet && reco_puppi_jet) dr_gen_puppi = deltaR(*gen_jet,*reco_puppi_jet);
    event.set(h_dR_gen_puppi, dr_gen_puppi);

    double dr_gen_chs(-999.);
    if(gen_jet && reco_chs_jet) dr_gen_chs = deltaR(*gen_jet,*reco_chs_jet);
    event.set(h_dR_gen_chs, dr_gen_chs);
  }

  return true;
}


std::vector<LorentzVector> WriteOutput::variedJets(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi){
  LorentzVector up, down;
  for(auto part:particles){
    LorentzVectorXYZE old_v4XYZ = toXYZ(part.v4());
    double puppiWeight = apply_puppi ? part.puppiWeight() : 1.0;
    PFParticle p = part;
    p.set_v4(toPtEtaPhi(old_v4XYZ * puppiWeight));
    double pt = p.pt();
    double eta = fabs(p.eta());
    int ptbin = grid->GetXaxis()->FindBin(pt);
    int etabin = grid->GetYaxis()->FindBin(eta);
    // if bin is overflow, set manually to last bin
    if(ptbin > Nbins_pt) ptbin = Nbins_pt;
    if(etabin > Nbins_eta) etabin = Nbins_eta;
    double factor = 1.0;
    if(isMC && ptbin == i && etabin == j && inCategory(p,cat)){
      up += (factor+variation) * p.v4();
      down += (factor-variation) * p.v4();
    }
    else{
      up += p.v4();
      down += p.v4();
    }
  }  
  return (std::vector<LorentzVector>) {up, down};  
}


std::vector<double> WriteOutput::CalculateMJetVariation1(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi){
  std::vector<LorentzVector> vP = variedJets(particles, i, j, cat, apply_puppi);
  return ( std::vector<double> {vP[0].M(),vP[1].M()});
}

    
std::vector<double> WriteOutput::CalculatePtVariation(std::vector<PFParticle> particles, int i, int j, TString cat, bool apply_puppi){
  std::vector<LorentzVector> vP = variedJets(particles, i, j, cat, apply_puppi);
  return ( std::vector<double> {vP[0].pt(),vP[1].pt()});
}

LorentzVector WriteOutput::nominalJet(std::vector<PFParticle> particles, bool apply_puppi){
  LorentzVector jet_v4;
  for(auto part:particles){
    LorentzVectorXYZE old_v4XYZ = toXYZ(part.v4());
    double puppiWeight = apply_puppi ? part.puppiWeight() : 1.0;
    PFParticle p = part;
    p.set_v4(toPtEtaPhi(old_v4XYZ * puppiWeight));
    jet_v4 += p.v4();
  }
  return jet_v4;
}

double WriteOutput::CalculateMJet(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto part:Particles){
    LorentzVectorXYZE old_v4XYZ = toXYZ(part.v4());
    double puppiWeight = part.puppiWeight();
    PFParticle p = part;
    p.set_v4(toPtEtaPhi(old_v4XYZ * puppiWeight));
    jet_v4 += p.v4();
  }
  double mjet = jet_v4.M();
  return mjet;
}



vector<double> WriteOutput::CalculateMJetVariation(vector<PFParticle>particles, int i, int j, TString cat){
  LorentzVector up, down;
  for(auto part:particles){
    LorentzVectorXYZE old_v4XYZ = toXYZ(part.v4());
    double puppiWeight = part.puppiWeight();
    PFParticle p = part;
    p.set_v4(toPtEtaPhi(old_v4XYZ * puppiWeight));
    double pt = p.pt();
    double eta = fabs(p.eta());
    int ptbin = grid->GetXaxis()->FindBin(pt);
    int etabin = grid->GetYaxis()->FindBin(eta);
    // if bin is overflow, set manually to last bin
    if(ptbin > Nbins_pt) ptbin = Nbins_pt;
    if(etabin > Nbins_eta) etabin = Nbins_eta;
    double factor = 1.0;
    if(isMC && ptbin == i && etabin == j && inCategory(p,cat)){
      up += (factor+variation) * p.v4();
      down += (factor-variation) * p.v4();
    }
    else{
      up += p.v4();
      down += p.v4();
    }
  }
  return (vector<double>) {up.M(), down.M()};
}

bool WriteOutput::inCategory(PFParticle p, TString cat){
  if(cat == "all") return true;
  bool inCat = false;
  int id = p.particleID();
  for(unsigned int i=0; i<all_cat.size(); i++){
    if(cat == all_cat[i] && abs(id) == i) inCat = true;
  }
  if(cat == "other"){
    for(auto oID: otherIDs){
      if(id == oID) inCat = true;
    }
  }
  return inCat;
}

void WriteOutput::ConstructOtherIDs(){
  vector<bool> AddIDToOther(9, true); //going to 8+1 since i added "all" to all_cat vector
  for(unsigned int i=0; i<all_cat.size(); i++){
    for(auto cat: categories){
      if(cat == all_cat[i]) AddIDToOther[i] = false;
    }
  }
  for(unsigned int id=0; id<AddIDToOther.size(); id++){
    if(AddIDToOther[id]) otherIDs.push_back(id);
  }
  return;
}


int WriteOutput::countBJetsAroundJet(const uhh2::Event & event, TopJet AK8, std::vector<Jet> AK4_jets, JetId btag, float deltaR_min ){
  int nbtag = 0;
  for(unsigned int ijet = 0; ijet<AK4_jets.size(); ijet++){
    Jet thisjet = AK4_jets.at(ijet);
    bool jettagged = btag(thisjet,event);
    bool not_overlapping = deltaR(AK8,thisjet)>=deltaR_min;
    if(jettagged && not_overlapping) nbtag++;
  }
  return nbtag;
}
