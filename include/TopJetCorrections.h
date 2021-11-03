#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"
#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"


class TopJetCorrections: public uhh2::AnalysisModule {
public:
    TopJetCorrections();

    virtual bool process(uhh2::Event & event) override;
    void init(uhh2::Context & ctx);
     void disable_jersmear(){jer_smearing_ = false;}

  //these are defined as static members here so the StandalineTopJetCorrector can access them as well.
  inline static const std::string tjec_tag_2016 = "Summer16_07Aug2017";
  inline static const std::string tjec_ver_2016 = "11";
  inline static const std::string tjec_tag_2017 = "Fall17_17Nov2017";
  inline static const std::string tjec_ver_2017 = "32";
  inline static const std::string tjec_tag_2018 = "Autumn18";
  inline static const std::string tjec_ver_2018 = "7";

  inline static const std::string tjec_tag_UL16preVFP = "Summer19UL16APV";
  inline static const std::string tjec_ver_UL16preVFP = "7";
  inline static const std::string tjec_tag_UL16postVFP = "Summer19UL16";
  inline static const std::string tjec_ver_UL16postVFP = "7";
  inline static const std::string tjec_tag_UL17 = "Summer19UL17";
  inline static const std::string tjec_ver_UL17 = "5";
  inline static const std::string tjec_tag_UL18 = "Summer19UL18";
  inline static const std::string tjec_ver_UL18 = "5";

  inline static const std::string tjer_tag_2016 = "Summer16_25nsV1";
  inline static const std::string tjer_tag_2017 = "Fall17_V3";
  inline static const std::string tjer_tag_2018 = "Autumn18_V7";

  inline static const std::string tjer_tag_UL16preVFP = "Summer20UL16APV_JRV3";
  inline static const std::string tjer_tag_UL16postVFP = "Summer20UL16_JRV3";
  inline static const std::string tjer_tag_UL17 = "Summer19UL17_JRV3";
  inline static const std::string tjer_tag_UL18 = "Summer19UL18_JRV2";
  

private:
    void fail_if_init() const;

    std::unique_ptr<YearSwitcher> tjet_corrector_MC, tjet_corrector_data;
    std::shared_ptr<RunSwitcher> tjec_switcher_16, tjec_switcher_17, tjec_switcher_18;
    std::shared_ptr<RunSwitcher> tjec_switcher_UL16preVFP, tjec_switcher_UL16postVFP, tjec_switcher_UL17, tjec_switcher_UL18;
    std::unique_ptr<GenericJetResolutionSmearer> tjet_resolution_smearer;

    bool is_mc;
    bool init_done = false;
    bool jer_smearing_;

    Year year;

    // Parameters for JEC & JLC sets
  std::string tjec_tjet_coll;
};

class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

 private:
  float minDR_;
};

class StandaloneTopJetCorrector{
public:
  StandaloneTopJetCorrector(uhh2::Context& ctx);
  float getJecFactor(const uhh2::Event & event, LorentzVector topjet);
  float getJERSmearingFactor(const uhh2::Event &event, Particle topjet, int direction, float jec_factor = 1.0);
  float getJERSmearingFactor(const uhh2::Event &event, LorentzVector topjet, int direction, float jec_factor = 1.0){
    Particle topjet_particle;
    topjet_particle.set_v4(topjet);
    return getJERSmearingFactor(event,topjet_particle,direction,jec_factor);
  };
  
private:
  bool is_mc;
  std::string short_year;
  std::map<std::string,std::unique_ptr<FactorizedJetCorrector>> jet_corrector;
  JME::JetResolution resolution_;
  JME::JetResolutionScaleFactor res_sf_;

};
