#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/JetMass/include/JetMassHists_central.h"
#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/JetMass/include/CorrectParticles.h"

using namespace std;
using namespace uhh2;

class TopMassModule: public AnalysisModule {
public:

  explicit TopMassModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  uhh2::Event::Handle<vector<PFParticle>> h_pfparticles;
  uhh2::Event::Handle<bool> h_passed_rec;
  std::unique_ptr<AnalysisModule> PUreweight, lumiweight, sf_btag, muo_tight_noniso_SF, muo_trigger_SF;
  std::unique_ptr<AnalysisModule> pf_jec;
  std::unique_ptr<AnalysisModule> pf_var1, pf_var2, pf_var3, pf_var4;
  std::unique_ptr<Selection> masscut;

  std::unique_ptr<Hists> h_pf;
  std::unique_ptr<Hists> h_mjet, h_mjet_pt400, h_mjet_pt300, h_mjet_pt200;
  std::unique_ptr<Hists> h_mjet_masscut, h_mjet_masscut_pt200, h_mjet_masscut_pt300, h_mjet_masscut_pt400;
  std::unique_ptr<Hists> h_mjet_failmasscut, h_mjet_failmasscut_pt200, h_mjet_failmasscut_pt300, h_mjet_failmasscut_pt400;
  std::unique_ptr<Hists> h_pf_aftervar1;

  std::unique_ptr<Hists> h_pf_var1;
  std::unique_ptr<Hists> h_mjet_var1, h_mjet_pt400_var1, h_mjet_pt300_var1, h_mjet_pt200_var1;
  std::unique_ptr<Hists> h_mjet_masscut_var1, h_mjet_masscut_pt200_var1, h_mjet_masscut_pt300_var1, h_mjet_masscut_pt400_var1;
  std::unique_ptr<Hists> h_mjet_failmasscut_var1, h_mjet_failmasscut_pt200_var1, h_mjet_failmasscut_pt300_var1, h_mjet_failmasscut_pt400_var1;

  std::unique_ptr<Hists> h_pf_var2;
  std::unique_ptr<Hists> h_mjet_var2, h_mjet_pt400_var2, h_mjet_pt300_var2, h_mjet_pt200_var2;
  std::unique_ptr<Hists> h_mjet_masscut_var2, h_mjet_masscut_pt200_var2, h_mjet_masscut_pt300_var2, h_mjet_masscut_pt400_var2;
  std::unique_ptr<Hists> h_mjet_failmasscut_var2, h_mjet_failmasscut_pt200_var2, h_mjet_failmasscut_pt300_var2, h_mjet_failmasscut_pt400_var2;

  std::unique_ptr<Hists> h_pf_var3;
  std::unique_ptr<Hists> h_mjet_var3, h_mjet_pt400_var3, h_mjet_pt300_var3, h_mjet_pt200_var3;
  std::unique_ptr<Hists> h_mjet_masscut_var3, h_mjet_masscut_pt200_var3, h_mjet_masscut_pt300_var3, h_mjet_masscut_pt400_var3;
  std::unique_ptr<Hists> h_mjet_failmasscut_var3, h_mjet_failmasscut_pt200_var3, h_mjet_failmasscut_pt300_var3, h_mjet_failmasscut_pt400_var3;

  std::unique_ptr<Hists> h_pf_var4;
  std::unique_ptr<Hists> h_mjet_var4, h_mjet_pt400_var4, h_mjet_pt300_var4, h_mjet_pt200_var4;
  std::unique_ptr<Hists> h_mjet_masscut_var4, h_mjet_masscut_pt200_var4, h_mjet_masscut_pt300_var4, h_mjet_masscut_pt400_var4;
  std::unique_ptr<Hists> h_mjet_failmasscut_var4, h_mjet_failmasscut_pt200_var4, h_mjet_failmasscut_pt300_var4, h_mjet_failmasscut_pt400_var4;


  bool isMC;

};


TopMassModule::TopMassModule(Context & ctx){

  isMC = (ctx.get("dataset_type") == "MC");
  h_passed_rec = ctx.get_handle<bool>("h_passed_rec");
  h_pfparticles = ctx.get_handle<std::vector<PFParticle>>("pfparticles");

  // apply JEC
  pf_jec.reset(new CorrectParticles());

  // apply some variation
  TString f1 = ctx.get("Grid_highptUP");
  TString f2 = ctx.get("Grid_netrualHDOWN_gammaDOWN");
  TString f3 = ctx.get("Grid_highetaDOWN");
  TString f4 = ctx.get("Grid_chargedHUP_oneBin");
  pf_var1.reset(new CorrectParticles(f1));
  pf_var2.reset(new CorrectParticles(f2));
  pf_var3.reset(new CorrectParticles(f3));
  pf_var4.reset(new CorrectParticles(f4));

  masscut.reset(new MassCut());

  lumiweight.reset(new MCLumiWeight(ctx));
  PUreweight.reset(new MCPileupReweight(ctx, "central"));
  sf_btag.reset(new MCBTagScaleFactor(ctx, BTag::DEEPCSV, BTag::WP_MEDIUM, "jets", "central"));
  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID", false, "central"));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger", false, "central"));

  double variation = 0.1;
  h_pf.reset(new PFHists(ctx, "PFHists"));
  h_mjet.reset(new JetMassHists(ctx, "JetMass", variation, "SD"));
  h_mjet_pt200.reset(new JetMassHists(ctx, "JetMass_pt200", variation, "SD"));
  h_mjet_pt300.reset(new JetMassHists(ctx, "JetMass_pt300", variation, "SD"));
  h_mjet_pt400.reset(new JetMassHists(ctx, "JetMass_pt400", variation, "SD"));
  h_mjet_masscut.reset(new JetMassHists(ctx, "JetMass_masscut", variation, "SD"));
  h_mjet_masscut_pt200.reset(new JetMassHists(ctx, "JetMass_masscut_pt200", variation, "SD"));
  h_mjet_masscut_pt300.reset(new JetMassHists(ctx, "JetMass_masscut_pt300", variation, "SD"));
  h_mjet_masscut_pt400.reset(new JetMassHists(ctx, "JetMass_masscut_pt400", variation, "SD"));
  h_mjet_failmasscut.reset(new JetMassHists(ctx, "JetMass_failmasscut", variation, "SD"));
  h_mjet_failmasscut_pt200.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt200", variation, "SD"));
  h_mjet_failmasscut_pt300.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300", variation, "SD"));
  h_mjet_failmasscut_pt400.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400", variation, "SD"));

  h_pf_aftervar1.reset(new PFHists(ctx, "PFHists_aftervar1"));

  h_pf_var1.reset(new PFHists(ctx, "PFHists_var1"));
  h_mjet_var1.reset(new JetMassHists_central(ctx, "JetMass_var1"));
  h_mjet_pt200_var1.reset(new JetMassHists_central(ctx, "JetMass_pt200_var1"));
  h_mjet_pt300_var1.reset(new JetMassHists_central(ctx, "JetMass_pt300_var1"));
  h_mjet_pt400_var1.reset(new JetMassHists_central(ctx, "JetMass_pt400_var1"));
  h_mjet_masscut_var1.reset(new JetMassHists_central(ctx, "JetMass_masscut_var1"));
  h_mjet_masscut_pt200_var1.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt200_var1"));
  h_mjet_masscut_pt300_var1.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt300_var1"));
  h_mjet_masscut_pt400_var1.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt400_var1"));
  h_mjet_failmasscut_var1.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_var1"));
  h_mjet_failmasscut_pt200_var1.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt200_var1"));
  h_mjet_failmasscut_pt300_var1.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt300_var1"));
  h_mjet_failmasscut_pt400_var1.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt400_var1"));

  h_pf_var2.reset(new PFHists(ctx, "PFHists_var2"));
  h_mjet_var2.reset(new JetMassHists_central(ctx, "JetMass_var2"));
  h_mjet_pt200_var2.reset(new JetMassHists_central(ctx, "JetMass_pt200_var2"));
  h_mjet_pt300_var2.reset(new JetMassHists_central(ctx, "JetMass_pt300_var2"));
  h_mjet_pt400_var2.reset(new JetMassHists_central(ctx, "JetMass_pt400_var2"));
  h_mjet_masscut_var2.reset(new JetMassHists_central(ctx, "JetMass_masscut_var2"));
  h_mjet_masscut_pt200_var2.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt200_var2"));
  h_mjet_masscut_pt300_var2.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt300_var2"));
  h_mjet_masscut_pt400_var2.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt400_var2"));
  h_mjet_failmasscut_var2.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_var2"));
  h_mjet_failmasscut_pt200_var2.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt200_var2"));
  h_mjet_failmasscut_pt300_var2.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt300_var2"));
  h_mjet_failmasscut_pt400_var2.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt400_var2"));

  h_pf_var3.reset(new PFHists(ctx, "PFHists_var3"));
  h_mjet_var3.reset(new JetMassHists_central(ctx, "JetMass_var3"));
  h_mjet_pt200_var3.reset(new JetMassHists_central(ctx, "JetMass_pt200_var3"));
  h_mjet_pt300_var3.reset(new JetMassHists_central(ctx, "JetMass_pt300_var3"));
  h_mjet_pt400_var3.reset(new JetMassHists_central(ctx, "JetMass_pt400_var3"));
  h_mjet_masscut_var3.reset(new JetMassHists_central(ctx, "JetMass_masscut_var3"));
  h_mjet_masscut_pt200_var3.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt200_var3"));
  h_mjet_masscut_pt300_var3.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt300_var3"));
  h_mjet_masscut_pt400_var3.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt400_var3"));
  h_mjet_failmasscut_var3.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_var3"));
  h_mjet_failmasscut_pt200_var3.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt200_var3"));
  h_mjet_failmasscut_pt300_var3.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt300_var3"));
  h_mjet_failmasscut_pt400_var3.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt400_var3"));

  h_pf_var4.reset(new PFHists(ctx, "PFHists_var4"));
  h_mjet_var4.reset(new JetMassHists_central(ctx, "JetMass_var4"));
  h_mjet_pt200_var4.reset(new JetMassHists_central(ctx, "JetMass_pt200_var4"));
  h_mjet_pt300_var4.reset(new JetMassHists_central(ctx, "JetMass_pt300_var4"));
  h_mjet_pt400_var4.reset(new JetMassHists_central(ctx, "JetMass_pt400_var4"));
  h_mjet_masscut_var4.reset(new JetMassHists_central(ctx, "JetMass_masscut_var4"));
  h_mjet_masscut_pt200_var4.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt200_var4"));
  h_mjet_masscut_pt300_var4.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt300_var4"));
  h_mjet_masscut_pt400_var4.reset(new JetMassHists_central(ctx, "JetMass_masscut_pt400_var4"));
  h_mjet_failmasscut_var4.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_var4"));
  h_mjet_failmasscut_pt200_var4.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt200_var4"));
  h_mjet_failmasscut_pt300_var4.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt300_var4"));
  h_mjet_failmasscut_pt400_var4.reset(new JetMassHists_central(ctx, "JetMass_failmasscut_pt400_var4"));
}


bool TopMassModule::process(Event & event) {

  if(isMC){
    lumiweight->process(event);
    PUreweight->process(event);
    sf_btag->process(event);
    muo_tight_noniso_SF->process(event);
    muo_trigger_SF->process(event);
  }

  bool passed_presel = event.get(h_passed_rec);

  // SELECTION
  vector<TopJet>* topjets = event.topjets;
  if(!passed_presel) return false;
  if(topjets->size() < 1) return false;
  if(topjets->at(0).pt() < 200) return false;

  // CORRECT PARTICLES
  pf_jec->process(event);
  h_pf->fill(event);

  // store input particles
  // (do not use pointer here!)
  vector<PFParticle> pf_input = event.get(h_pfparticles);

  // FILL MJET HISTOGRAMS
  double pt = topjets->at(0).pt();
  h_mjet->fill(event);
  if(pt > 200 && pt < 300) h_mjet_pt200->fill(event);
  else if(pt > 300 && pt < 400) h_mjet_pt300->fill(event);
  else if(pt > 400) h_mjet_pt400->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut->fill(event);
    if(pt > 200 && pt < 300) h_mjet_masscut_pt200->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_masscut_pt300->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400->fill(event);
  }
  else{
    h_mjet_failmasscut->fill(event);
    if(pt > 200 && pt < 300) h_mjet_failmasscut_pt200->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400->fill(event);
  }

  // include some variation
  pf_var1->process(event);
  h_pf_var1->fill(event);
  h_mjet_var1->fill(event);
  if(pt > 200 && pt < 300) h_mjet_pt200_var1->fill(event);
  else if(pt > 300 && pt < 400) h_mjet_pt300_var1->fill(event);
  else if(pt > 400) h_mjet_pt400_var1->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var1->fill(event);
    if(pt > 200 && pt < 300) h_mjet_masscut_pt200_var1->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var1->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var1->fill(event);
  }
  else{
    h_mjet_failmasscut_var1->fill(event);
    if(pt > 200 && pt < 300) h_mjet_failmasscut_pt200_var1->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var1->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var1->fill(event);
  }

  event.set(h_pfparticles, pf_input);// set pf back to input values
  h_pf_aftervar1->fill(event);

  // different variation
  pf_var2->process(event);
  h_pf_var2->fill(event);
  h_mjet_var2->fill(event);
  if(pt > 200 && pt < 300) h_mjet_pt200_var2->fill(event);
  else if(pt > 300 && pt < 400) h_mjet_pt300_var2->fill(event);
  else if(pt > 400) h_mjet_pt400_var2->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var2->fill(event);
    if(pt > 200 && pt < 300) h_mjet_masscut_pt200_var2->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var2->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var2->fill(event);
  }
  else{
    h_mjet_failmasscut_var2->fill(event);
    if(pt > 200 && pt < 300) h_mjet_failmasscut_pt200_var2->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var2->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var2->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values

  // different variation
  pf_var3->process(event);
  h_pf_var3->fill(event);
  h_mjet_var3->fill(event);
  if(pt > 200 && pt < 300) h_mjet_pt200_var3->fill(event);
  else if(pt > 300 && pt < 400) h_mjet_pt300_var3->fill(event);
  else if(pt > 400) h_mjet_pt400_var3->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var3->fill(event);
    if(pt > 200 && pt < 300) h_mjet_masscut_pt200_var3->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var3->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var3->fill(event);
  }
  else{
    h_mjet_failmasscut_var3->fill(event);
    if(pt > 200 && pt < 300) h_mjet_failmasscut_pt200_var3->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var3->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var3->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values

  // different variation
  pf_var4->process(event);
  h_pf_var4->fill(event);
  h_mjet_var4->fill(event);
  if(pt > 200 && pt < 300) h_mjet_pt200_var4->fill(event);
  else if(pt > 300 && pt < 400) h_mjet_pt300_var4->fill(event);
  else if(pt > 400) h_mjet_pt400_var4->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var4->fill(event);
    if(pt > 200 && pt < 300) h_mjet_masscut_pt200_var4->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var4->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var4->fill(event);
  }
  else{
    h_mjet_failmasscut_var4->fill(event);
    if(pt > 200 && pt < 300) h_mjet_failmasscut_pt200_var4->fill(event);
    else if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var4->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var4->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values
  // DONE
  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(TopMassModule)
