#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h" //mainly for runPeriods vectors

#include "UHH2/JetMass/include/SubstructureHists.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
#include "UHH2/JetMass/include/H3DDTHist.h"
#include "UHH2/JetMass/include/ECFHists.h"
#include "UHH2/common/include/MCLargeWeightKiller.h"

#include "TH3D.h"

#include <unistd.h>


#define EXTRAOUT false


using namespace std;
using namespace uhh2;

// bool repairPuppi(uhh2::Event & event){
//   vector<PFParticle>* particles = event.pfparticles;
//   for(unsigned int i=0; i<particles->size();i++){
//     LorentzVector old_v4 = particles->at(i).v4();
//     double factor = particles->at(i).puppiWeight();
//     particles->at(i).set_v4(old_v4 * 1/factor);
//   }
//   for (unsigned int i=0; i<particles->size();i++) {
//     LorentzVectorXYZE v4XYZ = toXYZ(particles->at(i).v4());
//     double puppiWeight = particles->at(i).puppiWeight();
//     particles->at(i).set_v4(toPtEtaPhi(v4XYZ * puppiWeight));
//   }
//   return true;
// }



namespace uhh2examples {

  class WMassModule: public AnalysisModule{
  public:
    explicit WMassModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> MC_LumiWeight,MC_PUWeight;
    std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;

    std::unique_ptr<Hists> h_PreSel,h_PreSel_Substr;
    std::unique_ptr<Hists> h_spikeKiller_1p5,h_spikeKiller_1p6,h_spikeKiller_1p7,h_spikeKiller_1p8,h_spikeKiller_2p0,h_spikeKiller_2p5;

    std::unique_ptr<Selection> n2ddtsel,n2sel,dak8ddtsel;

    std::unique_ptr<Selection> NEleSel,NMuonSel;

    std::unique_ptr<Hists> h_N2ddt_Substr_pass,h_N2ddt_Substr_fail;

    std::unique_ptr<Hists> h_N2ddt_pass,h_N2ddt_fail;
    std::unique_ptr<Hists> h_N2ddt_pt450To500_Sel_pass,h_N2ddt_pt450To500_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt500To550_Sel_pass,h_N2ddt_pt500To550_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt550To600_Sel_pass,h_N2ddt_pt550To600_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt600To675_Sel_pass,h_N2ddt_pt600To675_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt675To800_Sel_pass,h_N2ddt_pt675To800_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt800To1200_Sel_pass,h_N2ddt_pt800To1200_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt1200ToInf_Sel_pass,h_N2ddt_pt1200ToInf_Sel_fail;

    std::unique_ptr<Hists> h_DAK8ddt_pass,h_DAK8ddt_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt450To500_Sel_pass,h_DAK8ddt_pt450To500_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt500To550_Sel_pass,h_DAK8ddt_pt500To550_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt550To600_Sel_pass,h_DAK8ddt_pt550To600_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt600To675_Sel_pass,h_DAK8ddt_pt600To675_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt675To800_Sel_pass,h_DAK8ddt_pt675To800_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt800To1200_Sel_pass,h_DAK8ddt_pt800To1200_Sel_fail;
    std::unique_ptr<Hists> h_DAK8ddt_pt1200ToInf_Sel_pass,h_DAK8ddt_pt1200ToInf_Sel_fail;

    TH2F * n2_ddtMap, *dak8_ddtMap;
    

    // TH3D * N2_v_mSD_v_pT;
  };

  WMassModule::WMassModule(Context & ctx){
    is_mc = ctx.get("dataset_type") == "MC";
    is_WSample = ctx.get("dataset_version").find("WJets") != std::string::npos;
    is_ZSample = ctx.get("dataset_version").find("ZJets") != std::string::npos;
    
    if(is_mc){
      MC_LumiWeight.reset(new MCLumiWeight(ctx));
      MC_PUWeight.reset(new MCPileupReweight(ctx, "central"));
      mcSpikeKiller.reset(new MCLargeWeightKiller(ctx,infinity,infinity,infinity,2.5,3,3,"jets","genjets"));
    }

    apply_kfactors = string2bool(ctx.get("apply_kfactor", "true")) && (is_WSample ||is_ZSample);
    std::string ddt_FileName = ctx.get("ddtMapFile");
    TFile * mapFile = new TFile(ddt_FileName.c_str());

    std::string n2_ddt_HistName = ctx.get("N2ddtMapName");
    n2_ddtMap = (TH2F*) mapFile->Get(n2_ddt_HistName.c_str());
    n2ddtsel.reset(new N2ddtSelection(n2_ddtMap,true,TString(n2_ddt_HistName).Contains("PFMass")));

    std::string dak8_ddt_HistName = ctx.get("DAK8ddtMapName");
    dak8_ddtMap = (TH2F*) mapFile->Get(dak8_ddt_HistName.c_str());
    dak8ddtsel.reset(new DeepBoosted_WvsQCD_ddtSelection(dak8_ddtMap,true,TString(dak8_ddt_HistName).Contains("PFMass")));

    n2sel.reset(new N2Selection(0.2));

    if(is_mc)h_PreSel.reset(new JetMassHists(ctx,"JetMass_PreSel","SDVAR"));
    h_PreSel_Substr.reset(new SubstructureHists(ctx,"SubstructureHists_PreSel",n2_ddtMap,true));

    TString plotMode = "SDVAR";
    std::cout << "Creating Hists with mode: " << plotMode << std::endl;
    

    if(is_mc){
      h_N2ddt_Substr_pass.reset(new SubstructureHists(ctx,"N2ddt_SubstructureHists_pass",n2_ddtMap,true));
      h_N2ddt_Substr_fail.reset(new SubstructureHists(ctx,"N2ddt_SubstructureHists_fail",n2_ddtMap,true));
    }
    h_N2ddt_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pass",plotMode));
    h_N2ddt_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_fail",plotMode));
    
    h_N2ddt_pt450To500_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt450To500_pass",plotMode));
    h_N2ddt_pt450To500_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt450To500_fail",plotMode));
    h_N2ddt_pt500To550_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt500To550_pass",plotMode));
    h_N2ddt_pt500To550_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt500To550_fail",plotMode));
    h_N2ddt_pt550To600_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt550To600_pass",plotMode));
    h_N2ddt_pt550To600_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt550To600_fail",plotMode));
    h_N2ddt_pt600To675_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt600To675_pass",plotMode));
    h_N2ddt_pt600To675_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt600To675_fail",plotMode));
    h_N2ddt_pt675To800_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt675To800_pass",plotMode));
    h_N2ddt_pt675To800_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt675To800_fail",plotMode));
    h_N2ddt_pt800To1200_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt800To1200_pass",plotMode));
    h_N2ddt_pt800To1200_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt800To1200_fail",plotMode));
    h_N2ddt_pt1200ToInf_Sel_pass.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt1200ToInf_pass",plotMode));
    h_N2ddt_pt1200ToInf_Sel_fail.reset(new JetMassHists(ctx,"JetMass_N2DDT_pt1200ToInf_fail",plotMode));

    h_DAK8ddt_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pass",plotMode));
    h_DAK8ddt_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_fail",plotMode));
    
    h_DAK8ddt_pt450To500_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt450To500_pass",plotMode));
    h_DAK8ddt_pt450To500_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt450To500_fail",plotMode));
    h_DAK8ddt_pt500To550_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt500To550_pass",plotMode));
    h_DAK8ddt_pt500To550_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt500To550_fail",plotMode));
    h_DAK8ddt_pt550To600_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt550To600_pass",plotMode));
    h_DAK8ddt_pt550To600_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt550To600_fail",plotMode));
    h_DAK8ddt_pt600To675_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt600To675_pass",plotMode));
    h_DAK8ddt_pt600To675_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt600To675_fail",plotMode));
    h_DAK8ddt_pt675To800_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt675To800_pass",plotMode));
    h_DAK8ddt_pt675To800_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt675To800_fail",plotMode));
    h_DAK8ddt_pt800To1200_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt800To1200_pass",plotMode));
    h_DAK8ddt_pt800To1200_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt800To1200_fail",plotMode));
    h_DAK8ddt_pt1200ToInf_Sel_pass.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt1200ToInf_pass",plotMode));
    h_DAK8ddt_pt1200ToInf_Sel_fail.reset(new JetMassHists(ctx,"JetMass_DAK8DDT_pt1200ToInf_fail",plotMode));
	
  }

  bool WMassModule::process(Event & event){
    if(is_mc){
      MC_LumiWeight->process(event);
      MC_PUWeight->process(event);
      if(!mcSpikeKiller->passes(event)) return false;
    }
    if(is_mc)h_PreSel->fill(event);
    h_PreSel_Substr->fill(event);
    if(event.topjets->size() < 1 ) return false;
    double pt = event.topjets->at(0).pt();
    
    bool N2DDT_Selection = n2ddtsel->passes(event);

    if(N2DDT_Selection){
      if(is_mc){
        h_N2ddt_Substr_pass->fill(event);
      }
      h_N2ddt_pass->fill(event);
      if(pt > 450 && pt < 500)h_N2ddt_pt450To500_Sel_pass->fill(event);
      else if(pt >= 500 && pt < 550)h_N2ddt_pt500To550_Sel_pass->fill(event);
      else if(pt >= 550 && pt < 600)h_N2ddt_pt550To600_Sel_pass->fill(event);
      else if (pt >= 600 && pt < 675)h_N2ddt_pt600To675_Sel_pass->fill(event);
      else if (pt >= 675 && pt < 800)h_N2ddt_pt675To800_Sel_pass->fill(event);
      else if (pt >= 800 && pt < 1200)h_N2ddt_pt800To1200_Sel_pass->fill(event);
      else h_N2ddt_pt1200ToInf_Sel_pass->fill(event);
    }else{
      if(is_mc){
        h_N2ddt_Substr_fail->fill(event);
      }
      h_N2ddt_fail->fill(event);
      if(pt > 450 && pt < 500)h_N2ddt_pt450To500_Sel_fail->fill(event);
      else if(pt >= 500 && pt < 550)h_N2ddt_pt500To550_Sel_fail->fill(event);
      else if(pt >= 550 && pt < 600)h_N2ddt_pt550To600_Sel_fail->fill(event);
      else if (pt >= 600 && pt < 675)h_N2ddt_pt600To675_Sel_fail->fill(event);
      else if (pt >= 675 && pt < 800)h_N2ddt_pt675To800_Sel_fail->fill(event);
      else if (pt >= 800 && pt < 1200)h_N2ddt_pt800To1200_Sel_fail->fill(event);
      else h_N2ddt_pt1200ToInf_Sel_fail->fill(event);
    }

    bool DAK8DDT_Selection = dak8ddtsel->passes(event);

    if(DAK8DDT_Selection){
      h_DAK8ddt_pass->fill(event);
      if(pt > 450 && pt < 500)h_DAK8ddt_pt450To500_Sel_pass->fill(event);
      else if(pt >= 500 && pt < 550)h_DAK8ddt_pt500To550_Sel_pass->fill(event);
      else if(pt >= 550 && pt < 600)h_DAK8ddt_pt550To600_Sel_pass->fill(event);
      else if (pt >= 600 && pt < 675)h_DAK8ddt_pt600To675_Sel_pass->fill(event);
      else if (pt >= 675 && pt < 800)h_DAK8ddt_pt675To800_Sel_pass->fill(event);
      else if (pt >= 800 && pt < 1200)h_DAK8ddt_pt800To1200_Sel_pass->fill(event);
      else h_DAK8ddt_pt1200ToInf_Sel_pass->fill(event);
    }else{
      h_DAK8ddt_fail->fill(event);
      if(pt > 450 && pt < 500)h_DAK8ddt_pt450To500_Sel_fail->fill(event);
      else if(pt >= 500 && pt < 550)h_DAK8ddt_pt500To550_Sel_fail->fill(event);
      else if(pt >= 550 && pt < 600)h_DAK8ddt_pt550To600_Sel_fail->fill(event);
      else if (pt >= 600 && pt < 675)h_DAK8ddt_pt600To675_Sel_fail->fill(event);
      else if (pt >= 675 && pt < 800)h_DAK8ddt_pt675To800_Sel_fail->fill(event);
      else if (pt >= 800 && pt < 1200)h_DAK8ddt_pt800To1200_Sel_fail->fill(event);
      else h_DAK8ddt_pt1200ToInf_Sel_fail->fill(event);
    }

    
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(WMassModule)
}
