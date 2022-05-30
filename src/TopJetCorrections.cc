#include "../include/TopJetCorrections.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TriggerSelection.h"

using namespace uhh2;
using namespace std;

void TopJetCorrections::fail_if_init() const{
  if(init_done){
    throw invalid_argument("TopJetCorrections::init already called!");
  }
}


TopJetCorrections::TopJetCorrections(const std::string jet_collection_name):jet_collection_name_(jet_collection_name){
  tjec_tjet_coll = "dummy";
}


void TopJetCorrections::init(Context & ctx){
  if(init_done){
    throw invalid_argument("TopJetCorrections::init called twice!");
  }
  init_done = true;

  jer_smearing_ = true;

  is_mc = ctx.get("dataset_type") == "MC";
  year = extract_year(ctx);

  // setup correction tjet type for JECs
  std::string userTopJetColl = string2lowercase(ctx.get("TopJetCollection"));
  recojet_collection = "topjets";
  genjet_collection = "gentopjets";

  if(jet_collection_name_ != ""){
    recojet_collection = jet_collection_name_;
    userTopJetColl = string2lowercase(recojet_collection);
  }
  
  std::string algo = "";
  // algo size
  if (userTopJetColl.find("ak4") != std::string::npos) {
    algo = "AK4";
  }
  else if(userTopJetColl.find("ak8") != std::string::npos){
    algo = "AK8";
  }
  else if (userTopJetColl.find("ak8") == std::string::npos) {
    std::cout << "TopJetCorrections.cxx: Cannot determine tjet cone + radius (neither AK4 nor AK8) - going to assume it is AK8 for JECs" << '\n';
    algo = "AK8";
  }

  std::string pus = "PFchs";
  // Pileup subtraction
  if (userTopJetColl.find("puppi") != std::string::npos) {
    pus = "PFPuppi";
  } else if (userTopJetColl.find("chs") == std::string::npos) {
    std::cout << "Cannot determine pileup subtraction (neither CHS nor PUPPI) - going to assume it is CHS for JECs" << std::endl;
  }
  tjec_tjet_coll = algo + pus;

  if(is_mc){
    tjet_corrector_MC.reset(new YearSwitcher(ctx));

    tjet_corrector_MC->setup2016(       std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_2016,        tjec_ver_2016,        tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setup2017(       std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_2017,        tjec_ver_2017,        tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setup2018(       std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_2018,        tjec_ver_2018,        tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setupUL16preVFP( std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_UL16preVFP,  tjec_ver_UL16preVFP,  tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setupUL16postVFP(std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_UL16postVFP, tjec_ver_UL16postVFP, tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setupUL17(       std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_UL17,        tjec_ver_UL17,        tjec_tjet_coll),recojet_collection));
    tjet_corrector_MC->setupUL18(       std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(tjec_tag_UL18,        tjec_ver_UL18,        tjec_tjet_coll),recojet_collection));

    std::string jer_tag = "";
    if (year == Year::is2016v2 || year == Year::is2016v3) {
      jer_tag = tjer_tag_2016;
    } else if (year == Year::is2017v1 || year == Year::is2017v2) {
      jer_tag = tjer_tag_2017;
    } else if (year == Year::is2018) {
      jer_tag = tjer_tag_2018;
    } else if (year == Year::isUL16preVFP) {
      jer_tag = tjer_tag_UL16preVFP;
    } else if (year == Year::isUL16postVFP) {
      jer_tag = tjer_tag_UL16postVFP;
    } else if (year == Year::isUL17) {
      jer_tag = tjer_tag_UL17;
    } else if (year == Year::isUL18) {
      jer_tag = tjer_tag_UL18;
    } else {
      throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
    }

    // tjet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "gentopjets", JERFiles::JERPathStringMC(jer_tag,tjec_tjet_coll,"SF"), JERFiles::JERPathStringMC(jer_tag,tjec_tjet_coll,"PtResolution")));
    tjet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, recojet_collection, genjet_collection, JERFiles::JERPathStringMC(jer_tag,tjec_tjet_coll,"SF"), JERFiles::JERPathStringMC(jer_tag,tjec_tjet_coll,"PtResolution")));
  }
  else{
    tjec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      // tjec_switcher_16->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2016, tjec_ver_2016, tjec_tjet_coll, runItr)));
      tjec_switcher_16->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2016, tjec_ver_2016, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    for (const auto & runItr : runPeriods2017) {
      // tjec_switcher_17->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2017, tjec_ver_2017, tjec_tjet_coll, runItr)));
      tjec_switcher_17->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2017, tjec_ver_2017, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    for (const auto & runItr : runPeriods2018) {
      // tjec_switcher_18->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2018, tjec_ver_2018, tjec_tjet_coll, runItr)));
      tjec_switcher_18->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_2018, tjec_ver_2018, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_UL16preVFP.reset(new RunSwitcher(ctx, "UL16preVFP"));
    for (const auto & runItr : runPeriodsUL16preVFP) {
      tjec_switcher_UL16preVFP->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_UL16preVFP, tjec_ver_UL16preVFP, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_UL16postVFP.reset(new RunSwitcher(ctx, "UL16postVFP"));
    for (const auto & runItr : runPeriodsUL16postVFP) {
      tjec_switcher_UL16postVFP->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_UL16postVFP, tjec_ver_UL16postVFP, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_UL17.reset(new RunSwitcher(ctx, "UL17"));
    for (const auto & runItr : runPeriods2017) {
      tjec_switcher_UL17->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_UL17, tjec_ver_UL17, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjec_switcher_UL18.reset(new RunSwitcher(ctx, "UL18"));
    for (const auto & runItr : runPeriods2018) {
      tjec_switcher_UL18->setupRun(runItr, std::make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(tjec_tag_UL18, tjec_ver_UL18, tjec_tjet_coll, runItr),recojet_collection));
    }

    tjet_corrector_data.reset(new YearSwitcher(ctx));
    tjet_corrector_data->setup2016(tjec_switcher_16);
    tjet_corrector_data->setup2017(tjec_switcher_17);
    tjet_corrector_data->setup2018(tjec_switcher_18);

    tjet_corrector_data->setupUL16preVFP(tjec_switcher_UL16preVFP);
    tjet_corrector_data->setupUL16postVFP(tjec_switcher_UL16postVFP);
    tjet_corrector_data->setupUL17(tjec_switcher_UL17);
    tjet_corrector_data->setupUL18(tjec_switcher_UL18);
  }
}

bool TopJetCorrections::process(uhh2::Event & event){
  if(!init_done){
    throw runtime_error("TopJetCorrections::init not called (has to be called in AnalysisModule constructor)");
  }
  if (is_mc) {
    tjet_corrector_MC->process(event);
  } else {
    tjet_corrector_data->process(event);
  }

  if(is_mc && jer_smearing_) tjet_resolution_smearer->process(event);
  return true;
}


TopJetLeptonDeltaRCleaner::TopJetLeptonDeltaRCleaner(uhh2::Context & ctx, float mindr, const std::string jet_collection_name):minDR_(mindr),jet_collection_name_(jet_collection_name){
  h_topjets = ctx.get_handle<std::vector<TopJet> >(jet_collection_name_);  
}
bool TopJetLeptonDeltaRCleaner::process(uhh2::Event & event) {
  const auto topjets = &event.get(h_topjets);
  assert(topjets);
  
  std::vector<TopJet> cleaned_topjets;

  for(const auto & tjet : *topjets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
        if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
        if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }

  topjets->clear();
  topjets->reserve(cleaned_topjets.size());
  for(auto & j : cleaned_topjets) topjets->push_back(j);

  return true;
}


StandaloneTopJetCorrector::StandaloneTopJetCorrector(uhh2::Context& ctx, std::string tjec_tjet_coll){

  is_mc = ctx.get("dataset_type") == "MC";
  Year year = extract_year(ctx);

  const std::string year_str = year_str_map.at(year);
  if(year_str.find("16") != std::string::npos) short_year = "2016";
  if(year_str.find("17") != std::string::npos) short_year = "2017";
  if(year_str.find("18") != std::string::npos) short_year = "2018";
  
  std::cout << "Standalone TopJetCorrector for " << tjec_tjet_coll << std::endl;
  std::unordered_map<std::string,std::vector<std::string>> jec_mc_files = {
    {"2016",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_2016, TopJetCorrections::tjec_ver_2016, tjec_tjet_coll)},
    {"2017",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_2017, TopJetCorrections::tjec_ver_2017, tjec_tjet_coll)},
    {"2018",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_2018, TopJetCorrections::tjec_ver_2018, tjec_tjet_coll)},
    {"UL16preVFP",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_UL16preVFP, TopJetCorrections::tjec_ver_UL16preVFP, tjec_tjet_coll)},
    {"UL16postVFP",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_UL16postVFP, TopJetCorrections::tjec_ver_UL16postVFP, tjec_tjet_coll)},
    {"UL17",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_UL17, TopJetCorrections::tjec_ver_UL17, tjec_tjet_coll)},
    {"UL18",JERFiles::JECFilesMC(TopJetCorrections::tjec_tag_UL18, TopJetCorrections::tjec_ver_UL18, tjec_tjet_coll)}
  };
    
  std::unordered_map<std::string,std::map<std::string,std::vector<std::string>>> jec_data_files;
  std::map<std::string,std::vector<std::string>> jec_data_files_16,jec_data_files_17,jec_data_files_18;
  std::map<std::string,std::vector<std::string>> jec_data_files_UL16preVFP,jec_data_files_UL16postVFP,jec_data_files_UL17,jec_data_files_UL18;
  for (const auto runPeriod : year2runPeriods("2016")) {
    jec_data_files_16[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_2016, TopJetCorrections::tjec_ver_2016, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("2017")) {  
    jec_data_files_17[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_2017, TopJetCorrections::tjec_ver_2017, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("2018")) {
    jec_data_files_18[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_2018, TopJetCorrections::tjec_ver_2018, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("UL16preVFP")) {
    jec_data_files_UL16preVFP[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_UL16preVFP, TopJetCorrections::tjec_ver_UL16preVFP, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("UL16postVFP")) {
    jec_data_files_UL16postVFP[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_UL16postVFP, TopJetCorrections::tjec_ver_UL16postVFP, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("UL17")) {
    jec_data_files_UL17[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_UL17, TopJetCorrections::tjec_ver_UL17, tjec_tjet_coll,runPeriod);
  }
  for (const auto runPeriod : year2runPeriods("UL18")) {
    jec_data_files_UL18[runPeriod] = JERFiles::JECFilesDATA(TopJetCorrections::tjec_tag_UL18, TopJetCorrections::tjec_ver_UL18, tjec_tjet_coll,runPeriod);
  }
  jec_data_files["2016"] = jec_data_files_16;
  jec_data_files["2017"] = jec_data_files_17;
  jec_data_files["2018"] = jec_data_files_18;
  jec_data_files["UL16preVFP"] = jec_data_files_UL16preVFP;
  jec_data_files["UL16postVFP"] = jec_data_files_UL16postVFP;
  jec_data_files["UL17"] = jec_data_files_UL17;
  jec_data_files["UL18"] = jec_data_files_UL18;
    
  std::vector<JetCorrectorParameters> jec_pars_mc;
  for(const auto & filename: jec_mc_files.find(short_year)->second){
    jec_pars_mc.emplace_back(locate_file(filename));
  }      
  jet_corrector["MC"] = std::make_unique<FactorizedJetCorrector>(jec_pars_mc);

  for (const auto runPeriod : year2runPeriods(short_year)) {
    std::vector<JetCorrectorParameters> jec_pars_data;
    auto run_files_map = jec_data_files[short_year];
    for(const auto & filename: run_files_map.find(runPeriod)->second){
      jec_pars_data.emplace_back(locate_file(filename));
    }            
    jet_corrector[runPeriod] = std::make_unique<FactorizedJetCorrector>(jec_pars_data);    
  }

  std::unordered_map<std::string,std::string> jer_tag = {
    {"2016",TopJetCorrections::tjer_tag_2016},
    {"2017",TopJetCorrections::tjer_tag_2017},
    {"2018",TopJetCorrections::tjer_tag_2018},
    {"UL16preVFP",TopJetCorrections::tjer_tag_UL16preVFP},
    {"UL16postVFP",TopJetCorrections::tjer_tag_UL16postVFP},
    {"UL17",TopJetCorrections::tjer_tag_UL17},
    {"UL18",TopJetCorrections::tjer_tag_UL18}
  };
  res_sf_= JME::JetResolutionScaleFactor(locate_file(JERFiles::JERPathStringMC(jer_tag[short_year],tjec_tjet_coll,"SF")));
  resolution_ = JME::JetResolution(locate_file(JERFiles::JERPathStringMC(jer_tag[short_year],tjec_tjet_coll,"PtResolution")));
  
}

float StandaloneTopJetCorrector::getJecFactor(const uhh2::Event & event, LorentzVector topjet){
  std::string period;  
  if(is_mc){
    period="MC";
  }else{
    //This following bit is taken from the YearRunSwitcher
    auto foundYear = run_number_map.find(short_year);
    if(foundYear == run_number_map.end()){
      std::string valid = "";
      for (const auto & itr : run_number_map) {
        valid += itr.first;
        valid += ", ";      
      }
      throw std::runtime_error("year for CustomRunSwitcher (for SoftdropJetCorrector) must be one of: " + valid);
    }
    std::map<std::string, std::pair<int, int>> runNumberMap_ = foundYear->second;

    for (const auto & [key, val] : runNumberMap_) {
      if (event.run >= val.first && event.run <= val.second) {
        period = key;
      }
    }
  }
  
  jet_corrector[period]->setJetPt(topjet.pt());
  jet_corrector[period]->setJetEta(topjet.eta());
  jet_corrector[period]->setJetE(topjet.E());
  jet_corrector[period]->setJetA(0);
  jet_corrector[period]->setJetPhi(topjet.phi());
  jet_corrector[period]->setRho(event.rho);

  auto correction_factors = jet_corrector[period]->getSubCorrections();
  
  
  return correction_factors.back();
}

float StandaloneTopJetCorrector::getJERSmearingFactor(const uhh2::Event &event, Particle topjet, int direction,float jec_factor){
  if(! is_mc) return -1.0;
  float radius = 0.8;
  float recopt = topjet.pt()*jec_factor;
  float recoeta = topjet.eta();
  float abseta = fabs(recoeta);
  float rho = event.rho;

  assert(event.gentopjets);
  auto closest_genjet = closestParticle(topjet, *event.gentopjets);
  float genpt = -1.0;
  float resolution = resolution_.getResolution({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}, {JME::Binning::Rho, rho}});

  if (isnan(resolution)) {
    if (recopt < 35) { // leniency in this problematic region, hopefully fixed in future version of JER
      cout << "WARNING: getResolution() evaluated to nan. Since this jet is in problematic region, it will instead be set to 0." << endl;
      cout << "Input eta : rho : pt = " << recoeta << " : " << rho << ": " << recopt << endl;
      resolution = 0.;
    } else {
      throw std::runtime_error("getResolution() evaluated to nan. Input eta : rho : pt = " + double2string(recoeta) + " : " + double2string(rho) + " : " + double2string(recopt));
      }
  }
    if(!(closest_genjet == nullptr) && uhh2::deltaR(*closest_genjet, topjet) < 0.5*radius){
      genpt = closest_genjet->pt();
    }
    if( fabs(genpt-recopt) > 3*resolution*recopt){
      genpt=-1;
    }
    if(genpt < 15.0f) {
      genpt=-1.;
    }

    // Get the scale factor for this jet
    float c = -1.0;
    if (direction == 0) {
      c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}});
    } else if (direction == 1) {
      c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}}, Variation::UP);
    } else {
      c = res_sf_.getScaleFactor({{JME::Binning::JetPt, recopt}, {JME::Binning::JetEta, recoeta}}, Variation::DOWN);
    }
    
    if (c < 0) {
      std::cout << "WARNING: GenericJetResolutionSmearer: no scale factor found for this jet with pt : eta = " << recopt << " : " << recoeta << std::endl;
      std::cout << "         No JER smearing will be applied." << std::endl;
      return c;
    }

    // Calculate the new pt
    float new_pt = -1.;
    // Use scaling method in case a matching generator jet was found
    if(genpt>0){
      new_pt = std::max(0.0f, genpt + c * (recopt - genpt));
    }
    // Use stochastic method if no generator jet could be matched to the reco jet
    else{
      // Initialize random generator with eta-dependend random seed to be reproducible
      TRandom rand((int)(1000*abseta));
      float random_gauss = rand.Gaus(0, resolution);
      new_pt = recopt * (1 + random_gauss*sqrt(std::max(c*c-1, 0.0f)));
    }
  return new_pt / recopt;
}
