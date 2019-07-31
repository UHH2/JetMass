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

class CorrectParticles: public uhh2::AnalysisModule {
public:
    CorrectParticles(TString filename);
    CorrectParticles();  
    virtual bool process(uhh2::Event & ) override;
private:
    double GetGridFactor(PFParticle p);
    vector<TH2F*> grid;
    bool use_JECfactor = false;
    bool use_grid = false;
    vector<TString> all_cat = {"undef", "chargedH", "elec", "muon", "gamma", "neutralH", "had", "em"};
};
