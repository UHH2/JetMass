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

#include "TH3D.h"

#include <unistd.h>


#define EXTRAOUT false


using namespace std;
using namespace uhh2;

bool repairPuppi(uhh2::Event & event){
  vector<PFParticle>* particles = event.pfparticles;
  for(unsigned int i=0; i<particles->size();i++){
    LorentzVector old_v4 = particles->at(i).v4();
    double factor = particles->at(i).puppiWeight();
    particles->at(i).set_v4(old_v4 * 1/factor);
  }
  for (unsigned int i=0; i<particles->size();i++) {
    LorentzVectorXYZE v4XYZ = toXYZ(particles->at(i).v4());
    double puppiWeight = particles->at(i).puppiWeight();
    particles->at(i).set_v4(toPtEtaPhi(v4XYZ * puppiWeight));
  }
  return true;
}



namespace uhh2examples {

  class WMassModule: public AnalysisModule{
  public:
    explicit WMassModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> MC_LumiWeight,MC_PUWeight;

    std::unique_ptr<Hists> h_PreSel,h_PreSel_Substr,h_PreSel_Substr_MyDDT;
    std::unique_ptr<Hists> h_spikeKiller_1p5,h_spikeKiller_1p6,h_spikeKiller_1p7,h_spikeKiller_1p8,h_spikeKiller_2p0,h_spikeKiller_2p5;

    std::unique_ptr<Selection> n2ddtsel,MYn2ddtsel,n2sel;

    std::unique_ptr<Selection> NEleSel,NMuonSel;

    std::unique_ptr<Hists> h_MyN2ddt_Substr_pass,h_MyN2ddt_Substr_fail;

    std::unique_ptr<Hists> h_MyN2ddt_pass,h_MyN2ddt_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt450To500_Sel_pass,h_MyN2ddt_pt450To500_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt500To550_Sel_pass,h_MyN2ddt_pt500To550_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt550To600_Sel_pass,h_MyN2ddt_pt550To600_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt600To675_Sel_pass,h_MyN2ddt_pt600To675_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt675To800_Sel_pass,h_MyN2ddt_pt675To800_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt800To1200_Sel_pass,h_MyN2ddt_pt800To1200_Sel_fail;
    std::unique_ptr<Hists> h_MyN2ddt_pt1200ToInf_Sel_pass,h_MyN2ddt_pt1200ToInf_Sel_fail;

    std::unique_ptr<Hists> h_N2ddt_Substr_pass,h_N2ddt_Substr_fail;

    std::unique_ptr<Hists> h_N2ddt_pass,h_N2ddt_fail;
    std::unique_ptr<Hists> h_N2ddt_pt450To500_Sel_pass,h_N2ddt_pt450To500_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt500To550_Sel_pass,h_N2ddt_pt500To550_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt550To600_Sel_pass,h_N2ddt_pt550To600_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt600To675_Sel_pass,h_N2ddt_pt600To675_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt675To800_Sel_pass,h_N2ddt_pt675To800_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt800To1200_Sel_pass,h_N2ddt_pt800To1200_Sel_fail;
    std::unique_ptr<Hists> h_N2ddt_pt1200ToInf_Sel_pass,h_N2ddt_pt1200ToInf_Sel_fail;

    TH2F * ddtMap;
    TH2F * MYddtMap;
    bool is_mc;
    // TH3D * N2_v_mSD_v_pT;
  };

  WMassModule::WMassModule(Context & ctx){
    is_mc = ctx.get("dataset_type") == "MC";

    if(is_mc){
      MC_LumiWeight.reset(new MCLumiWeight(ctx));
      MC_PUWeight.reset(new MCPileupReweight(ctx, "central"));
    }

    std::string h2_ddt_FileName = ctx.get("ddtMapFile");
    std::string My_h2_ddt_FileName = ctx.get("MYddtMapFile");
    std::string h2_ddt_HistName = ctx.get("ddtMapName");
    TFile * mapFile = new TFile(h2_ddt_FileName.c_str());
    ddtMap = (TH2F*) mapFile->Get(h2_ddt_HistName.c_str());
    TFile * myMapFile = new TFile(My_h2_ddt_FileName.c_str());
    MYddtMap = (TH2F*) myMapFile->Get(h2_ddt_HistName.c_str());


    // std::string h3_fileName = ctx.get("h3File");
    // std::string h3_histName = ctx.get("h3Name");
    // TFile * mapFile = new TFile(h3_fileName.c_str(),"READ");
    // N2_v_mSD_v_pT = (TH3D*) mapFile->Get(h3_histName.c_str());

    // TH2D* ddtMap_QCD_0p01=SubstructureHists::computeDDTMap(0.01, N2_v_mSD_v_pT);


    MYn2ddtsel.reset(new N2ddtSelection(MYddtMap,true));
    n2ddtsel.reset(new N2ddtSelection(ddtMap,true));
    n2sel.reset(new N2Selection(0.2));
    // if(is_mc)h_ECF_comparison.reset(new ECFHists(ctx,"ECFComparisonHists"));

    if(is_mc)h_PreSel.reset(new JetMassHists(ctx,"JetMass_PreSel","SDVAR"));

    h_PreSel_Substr_MyDDT.reset(new SubstructureHists(ctx,"SubstructureHists_PreSel_MyDDT",MYddtMap,true));
    h_PreSel_Substr.reset(new SubstructureHists(ctx,"SubstructureHists_PreSel",ddtMap,true));

    TString plotMode = "";
    if(is_mc)plotMode = "SDVAR";
    else plotMode = "SDNC";
    std::cout << "Creating Hists with mode: " << plotMode << std::endl;
    
    
    if(is_mc){
      h_MyN2ddt_Substr_pass.reset(new SubstructureHists(ctx,"MyN2ddt_SubstructureHists_pass",ddtMap,true));
      h_MyN2ddt_Substr_fail.reset(new SubstructureHists(ctx,"MyN2ddt_SubstructureHists_fail",ddtMap,true));
    }
    
    h_MyN2ddt_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pass",plotMode));
    h_MyN2ddt_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_fail",plotMode));
    
    h_MyN2ddt_pt450To500_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt450To500_pass",plotMode));
    h_MyN2ddt_pt450To500_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt450To500_fail",plotMode));
    h_MyN2ddt_pt500To550_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt500To550_pass",plotMode));
    h_MyN2ddt_pt500To550_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt500To550_fail",plotMode));
    h_MyN2ddt_pt550To600_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt550To600_pass",plotMode));
    h_MyN2ddt_pt550To600_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt550To600_fail",plotMode));
    h_MyN2ddt_pt600To675_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt600To675_pass",plotMode));
    h_MyN2ddt_pt600To675_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt600To675_fail",plotMode));
    h_MyN2ddt_pt675To800_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt675To800_pass",plotMode));
    h_MyN2ddt_pt675To800_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt675To800_fail",plotMode));
    h_MyN2ddt_pt800To1200_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt800To1200_pass",plotMode));
    h_MyN2ddt_pt800To1200_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt800To1200_fail",plotMode));
    h_MyN2ddt_pt1200ToInf_Sel_pass.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt1200ToInf_pass",plotMode));
    h_MyN2ddt_pt1200ToInf_Sel_fail.reset(new JetMassHists(ctx,"JetMass_MyN2ddt_pt1200ToInf_fail",plotMode));
    

    if(is_mc){
      h_N2ddt_Substr_pass.reset(new SubstructureHists(ctx,"N2ddt_SubstructureHists_pass",ddtMap,true));
      h_N2ddt_Substr_fail.reset(new SubstructureHists(ctx,"N2ddt_SubstructureHists_fail",ddtMap,true));
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
	

  }

  bool WMassModule::process(Event & event){
    if(is_mc){
      MC_LumiWeight->process(event);
      MC_PUWeight->process(event);
    }
    if(is_mc)h_PreSel->fill(event);
    h_PreSel_Substr->fill(event);
    h_PreSel_Substr_MyDDT->fill(event);
    if(event.topjets->size() < 1 ) return false;
    double pt = event.topjets->at(0).pt();

    double pTave=0.0;
    double ratio=0.0;
    if(is_mc){
      if(event.jets->size()>0){
        pTave=event.jets->at(0).pt();
      }else	if(event.jets->size()>1){
        pTave=(event.jets->at(0).pt()+event.jets->at(1).pt())/2.0;
      }
      double pThat=event.genInfo->qScale();
      ratio=pThat/pTave;
    }

    double qScale_max=1.5;
    if(is_mc && (ratio > qScale_max)) return false;

    bool MyN2DDT_Selection = MYn2ddtsel->passes(event);

    if(MyN2DDT_Selection){
      if(is_mc) h_MyN2ddt_Substr_pass->fill(event);
      h_MyN2ddt_pass->fill(event);
      if(pt > 450 && pt < 500)h_MyN2ddt_pt450To500_Sel_pass->fill(event);
      else if(pt >= 500 && pt < 550)h_MyN2ddt_pt500To550_Sel_pass->fill(event);
      else if(pt >= 550 && pt < 600)h_MyN2ddt_pt550To600_Sel_pass->fill(event);
      else if (pt >= 600 && pt < 675)h_MyN2ddt_pt600To675_Sel_pass->fill(event);
      else if (pt >= 675 && pt < 800)h_MyN2ddt_pt675To800_Sel_pass->fill(event);
      else if (pt >= 800 && pt < 1200)h_MyN2ddt_pt800To1200_Sel_pass->fill(event);
      else h_MyN2ddt_pt1200ToInf_Sel_pass->fill(event);
    }else{
      if(is_mc) h_MyN2ddt_Substr_fail->fill(event);
      h_MyN2ddt_fail->fill(event);
      if(pt > 450 && pt < 500)h_MyN2ddt_pt450To500_Sel_fail->fill(event);
      else if(pt >= 500 && pt < 550)h_MyN2ddt_pt500To550_Sel_fail->fill(event);
      else if(pt >= 550 && pt < 600)h_MyN2ddt_pt550To600_Sel_fail->fill(event);
      else if (pt >= 600 && pt < 675)h_MyN2ddt_pt600To675_Sel_fail->fill(event);
      else if (pt >= 675 && pt < 800)h_MyN2ddt_pt675To800_Sel_fail->fill(event);
      else if (pt >= 800 && pt < 1200)h_MyN2ddt_pt800To1200_Sel_fail->fill(event);
      else h_MyN2ddt_pt1200ToInf_Sel_fail->fill(event);
    }
            
    bool N2DDT_Selection = n2ddtsel->passes(event);
    
    if(N2DDT_Selection){
      if(is_mc){
        h_N2ddt_pass->fill(event);
        h_N2ddt_Substr_pass->fill(event);
      }
      if(pt > 450 && pt < 500)h_N2ddt_pt450To500_Sel_pass->fill(event);
      else if(pt >= 500 && pt < 550)h_N2ddt_pt500To550_Sel_pass->fill(event);
      else if(pt >= 550 && pt < 600)h_N2ddt_pt550To600_Sel_pass->fill(event);
      else if (pt >= 600 && pt < 675)h_N2ddt_pt600To675_Sel_pass->fill(event);
      else if (pt >= 675 && pt < 800)h_N2ddt_pt675To800_Sel_pass->fill(event);
      else if (pt >= 800 && pt < 1200)h_N2ddt_pt800To1200_Sel_pass->fill(event);
      else h_N2ddt_pt1200ToInf_Sel_pass->fill(event);
    }else{
      if(is_mc){
        h_N2ddt_fail->fill(event);
        h_N2ddt_Substr_fail->fill(event);
      }
      if(pt > 450 && pt < 500)h_N2ddt_pt450To500_Sel_fail->fill(event);
      else if(pt >= 500 && pt < 550)h_N2ddt_pt500To550_Sel_fail->fill(event);
      else if(pt >= 550 && pt < 600)h_N2ddt_pt550To600_Sel_fail->fill(event);
      else if (pt >= 600 && pt < 675)h_N2ddt_pt600To675_Sel_fail->fill(event);
      else if (pt >= 675 && pt < 800)h_N2ddt_pt675To800_Sel_fail->fill(event);
      else if (pt >= 800 && pt < 1200)h_N2ddt_pt800To1200_Sel_fail->fill(event);
      else h_N2ddt_pt1200ToInf_Sel_fail->fill(event);
    }
    
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(WMassModule)
}
