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

class TopJetCorrections: public uhh2::AnalysisModule {
public:
    TopJetCorrections();

    virtual bool process(uhh2::Event & event) override;
    void init(uhh2::Context & ctx);

private:
    void fail_if_init() const;

    std::unique_ptr<YearSwitcher> tjet_corrector_MC, tjet_corrector_data;
    std::shared_ptr<RunSwitcher> tjec_switcher_16, tjec_switcher_17, tjec_switcher_18;
    std::unique_ptr<GenericJetResolutionSmearer> tjet_resolution_smearer;

    bool is_mc;
    bool init_done = false;

    Year year;

    // Parameters for JEC & JLC sets
    std::string tjec_tag_2016, tjec_ver_2016;
    std::string tjec_tag_2017, tjec_ver_2017;
    std::string tjec_tag_2018, tjec_ver_2018;
    std::string tjec_tjet_coll;
};
