#pragma once
#include "UHH2/core/include/fwd.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include <vector>
#include <TH2.h>
#include <TString.h>
#include <TFile.h>

using namespace std;

class ApplyPuppiToPF: public uhh2::AnalysisModule {
public:
    ApplyPuppiToPF();
    virtual bool process(uhh2::Event & ) override;
private:
};
