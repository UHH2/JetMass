#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/JetMass/include/EFJetHists.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"
#include "UHH2/JetMass/include/ECFHists.h"

using namespace std;
using namespace uhh2;

class TopMassModule: public AnalysisModule {
public:

  explicit TopMassModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  // Handles
  uhh2::Event::Handle<vector<PFParticle>> h_pfparticles;

  // Scale Factors
  std::unique_ptr<AnalysisModule> topreweight;
  std::unique_ptr<AnalysisModule> PUreweight, lumiweight, sf_btag, muo_tight_noniso_SF, muo_trigger_SF;

  //Selections
  std::unique_ptr<Selection> masscut;

  // PF variations
  std::unique_ptr<AnalysisModule> pf_var1, pf_var2, pf_var3, pf_var4;

  // Histograms
  std::unique_ptr<Hists> h_ecf;

  std::unique_ptr<Hists> h_EnergyFractions, h_EnergyFractions_W, h_EnergyFractions_top;
  std::unique_ptr<Hists> h_EnergyFractions_pt300, h_EnergyFractions_pt400;
  std::unique_ptr<Hists> h_mjet_loCHF_pt300, h_mjet_hiCHF_pt300;
  std::unique_ptr<Hists> h_mjet_loCHF_pt400, h_mjet_hiCHF_pt400;

  std::unique_ptr<Hists> h_pf, h_pf_W, h_pf_top;
  std::unique_ptr<Hists> h_mjet, h_mjet_pt400, h_mjet_pt300;
  std::unique_ptr<Hists> h_mjet_masscut, h_mjet_masscut_pt300, h_mjet_masscut_pt400;
  std::unique_ptr<Hists> h_mjet_failmasscut, h_mjet_failmasscut_pt300, h_mjet_failmasscut_pt400;
  std::unique_ptr<Hists> h_pf_aftervar1;

  std::unique_ptr<Hists> h_pf_var1;
  std::unique_ptr<Hists> h_mjet_var1, h_mjet_pt400_var1, h_mjet_pt300_var1;
  std::unique_ptr<Hists> h_mjet_masscut_var1, h_mjet_masscut_pt300_var1, h_mjet_masscut_pt400_var1;
  std::unique_ptr<Hists> h_mjet_failmasscut_var1, h_mjet_failmasscut_pt300_var1, h_mjet_failmasscut_pt400_var1;

  std::unique_ptr<Hists> h_pf_var2;
  std::unique_ptr<Hists> h_mjet_var2, h_mjet_pt400_var2, h_mjet_pt300_var2;
  std::unique_ptr<Hists> h_mjet_masscut_var2, h_mjet_masscut_pt300_var2, h_mjet_masscut_pt400_var2;
  std::unique_ptr<Hists> h_mjet_failmasscut_var2, h_mjet_failmasscut_pt300_var2, h_mjet_failmasscut_pt400_var2;

  std::unique_ptr<Hists> h_pf_var3;
  std::unique_ptr<Hists> h_mjet_var3, h_mjet_pt400_var3, h_mjet_pt300_var3;
  std::unique_ptr<Hists> h_mjet_masscut_var3, h_mjet_masscut_pt300_var3, h_mjet_masscut_pt400_var3;
  std::unique_ptr<Hists> h_mjet_failmasscut_var3, h_mjet_failmasscut_pt300_var3, h_mjet_failmasscut_pt400_var3;

  std::unique_ptr<Hists> h_pf_var4;
  std::unique_ptr<Hists> h_mjet_var4, h_mjet_pt400_var4, h_mjet_pt300_var4;
  std::unique_ptr<Hists> h_mjet_masscut_var4, h_mjet_masscut_pt300_var4, h_mjet_masscut_pt400_var4;
  std::unique_ptr<Hists> h_mjet_failmasscut_var4, h_mjet_failmasscut_pt300_var4, h_mjet_failmasscut_pt400_var4;

  // Variables
  bool isMC;

};


TopMassModule::TopMassModule(Context & ctx){

  isMC = (ctx.get("dataset_type") == "MC");
  h_pfparticles = ctx.get_handle<std::vector<PFParticle>>("pfparticles");

  // Top pT reweighting
  topreweight.reset(new TopPtReweight(ctx, 0.0615, 0.0005));

  // apply some variation
  TString f1 = ctx.get("Grid_highptUP");
  TString f2 = ctx.get("Grid_netrualHDOWN_gammaDOWN");
  TString f3 = ctx.get("Grid_highetaDOWN");
  TString f4 = ctx.get("Grid_chargedHUP_oneBin");
  pf_var1.reset(new CorrectParticles(f1));
  pf_var2.reset(new CorrectParticles(f2));
  pf_var3.reset(new CorrectParticles(f3));
  pf_var4.reset(new CorrectParticles(f4));

  // Selection
  masscut.reset(new MassCut());

  // Scale Factors
  lumiweight.reset(new MCLumiWeight(ctx));
  PUreweight.reset(new MCPileupReweight(ctx, "central"));
  sf_btag.reset(new MCBTagScaleFactor(ctx, BTag::DEEPCSV, BTag::WP_MEDIUM, "jets", "central"));
  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root","NUM_TightID_DEN_genTracks_eta_pt",1, "tightID", false, "central"));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger", false, "central"));

  // Histograms
  h_ecf.reset(new ECFHists(ctx, "ECFHists"));

  h_EnergyFractions.reset(new EFJetHists(ctx, "EnergyFractions"));
  h_EnergyFractions_W.reset(new EFJetHists(ctx, "EnergyFractions_w"));
  h_EnergyFractions_top.reset(new EFJetHists(ctx, "EnergyFractions_top"));

  h_EnergyFractions_pt300.reset(new EFJetHists(ctx, "h_EnergyFractions_pt300"));
  h_EnergyFractions_pt400.reset(new EFJetHists(ctx, "h_EnergyFractions_pt400"));

  h_mjet_loCHF_pt300.reset(new JetMassHists(ctx, "JetMass_loCHF_pt300", "SD"));
  h_mjet_hiCHF_pt300.reset(new JetMassHists(ctx, "JetMass_hiCHF_pt300", "SD"));
  h_mjet_loCHF_pt400.reset(new JetMassHists(ctx, "JetMass_loCHF_pt400", "SD"));
  h_mjet_hiCHF_pt400.reset(new JetMassHists(ctx, "JetMass_hiCHF_pt400", "SD"));

  h_pf.reset(new PFHists(ctx, "PFHists"));
  h_pf_W.reset(new PFHists(ctx, "PFHists_W"));
  h_pf_top.reset(new PFHists(ctx, "PFHists_top"));

  h_mjet.reset(new JetMassHists(ctx, "JetMass", "SDVAR"));
  h_mjet_pt300.reset(new JetMassHists(ctx, "JetMass_pt300", "SDVAR"));
  h_mjet_pt400.reset(new JetMassHists(ctx, "JetMass_pt400", "SDVAR"));
  h_mjet_masscut.reset(new JetMassHists(ctx, "JetMass_masscut", "SDVAR"));
  h_mjet_masscut_pt300.reset(new JetMassHists(ctx, "JetMass_masscut_pt300", "SDVAR"));
  h_mjet_masscut_pt400.reset(new JetMassHists(ctx, "JetMass_masscut_pt400", "SDVAR"));
  h_mjet_failmasscut.reset(new JetMassHists(ctx, "JetMass_failmasscut", "SDVAR"));
  h_mjet_failmasscut_pt300.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300", "SDVAR"));
  h_mjet_failmasscut_pt400.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400", "SDVAR"));

  h_pf_aftervar1.reset(new PFHists(ctx, "PFHists_aftervar1"));

  h_pf_var1.reset(new PFHists(ctx, "PFHists_var1"));
  h_mjet_var1.reset(new JetMassHists(ctx, "JetMass_var1", "SD"));
  h_mjet_pt300_var1.reset(new JetMassHists(ctx, "JetMass_pt300_var1", "SD"));
  h_mjet_pt400_var1.reset(new JetMassHists(ctx, "JetMass_pt400_var1", "SD"));
  h_mjet_masscut_var1.reset(new JetMassHists(ctx, "JetMass_masscut_var1", "SD"));
  h_mjet_masscut_pt300_var1.reset(new JetMassHists(ctx, "JetMass_masscut_pt300_var1", "SD"));
  h_mjet_masscut_pt400_var1.reset(new JetMassHists(ctx, "JetMass_masscut_pt400_var1", "SD"));
  h_mjet_failmasscut_var1.reset(new JetMassHists(ctx, "JetMass_failmasscut_var1", "SD"));
  h_mjet_failmasscut_pt300_var1.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300_var1", "SD"));
  h_mjet_failmasscut_pt400_var1.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400_var1", "SD"));

  h_pf_var2.reset(new PFHists(ctx, "PFHists_var2"));
  h_mjet_var2.reset(new JetMassHists(ctx, "JetMass_var2", "SD"));
  h_mjet_pt300_var2.reset(new JetMassHists(ctx, "JetMass_pt300_var2", "SD"));
  h_mjet_pt400_var2.reset(new JetMassHists(ctx, "JetMass_pt400_var2", "SD"));
  h_mjet_masscut_var2.reset(new JetMassHists(ctx, "JetMass_masscut_var2", "SD"));
  h_mjet_masscut_pt300_var2.reset(new JetMassHists(ctx, "JetMass_masscut_pt300_var2", "SD"));
  h_mjet_masscut_pt400_var2.reset(new JetMassHists(ctx, "JetMass_masscut_pt400_var2", "SD"));
  h_mjet_failmasscut_var2.reset(new JetMassHists(ctx, "JetMass_failmasscut_var2", "SD"));
  h_mjet_failmasscut_pt300_var2.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300_var2", "SD"));
  h_mjet_failmasscut_pt400_var2.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400_var2", "SD"));

  h_pf_var3.reset(new PFHists(ctx, "PFHists_var3"));
  h_mjet_var3.reset(new JetMassHists(ctx, "JetMass_var3", "SD"));
  h_mjet_pt300_var3.reset(new JetMassHists(ctx, "JetMass_pt300_var3", "SD"));
  h_mjet_pt400_var3.reset(new JetMassHists(ctx, "JetMass_pt400_var3", "SD"));
  h_mjet_masscut_var3.reset(new JetMassHists(ctx, "JetMass_masscut_var3", "SD"));
  h_mjet_masscut_pt300_var3.reset(new JetMassHists(ctx, "JetMass_masscut_pt300_var3", "SD"));
  h_mjet_masscut_pt400_var3.reset(new JetMassHists(ctx, "JetMass_masscut_pt400_var3", "SD"));
  h_mjet_failmasscut_var3.reset(new JetMassHists(ctx, "JetMass_failmasscut_var3", "SD"));
  h_mjet_failmasscut_pt300_var3.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300_var3", "SD"));
  h_mjet_failmasscut_pt400_var3.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400_var3", "SD"));

  h_pf_var4.reset(new PFHists(ctx, "PFHists_var4"));
  h_mjet_var4.reset(new JetMassHists(ctx, "JetMass_var4", "SD"));
  h_mjet_pt300_var4.reset(new JetMassHists(ctx, "JetMass_pt300_var4", "SD"));
  h_mjet_pt400_var4.reset(new JetMassHists(ctx, "JetMass_pt400_var4", "SD"));
  h_mjet_masscut_var4.reset(new JetMassHists(ctx, "JetMass_masscut_var4", "SD"));
  h_mjet_masscut_pt300_var4.reset(new JetMassHists(ctx, "JetMass_masscut_pt300_var4", "SD"));
  h_mjet_masscut_pt400_var4.reset(new JetMassHists(ctx, "JetMass_masscut_pt400_var4", "SD"));
  h_mjet_failmasscut_var4.reset(new JetMassHists(ctx, "JetMass_failmasscut_var4", "SD"));
  h_mjet_failmasscut_pt300_var4.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt300_var4", "SD"));
  h_mjet_failmasscut_pt400_var4.reset(new JetMassHists(ctx, "JetMass_failmasscut_pt400_var4", "SD"));
}


bool TopMassModule::process(Event & event) {
  // first, get topjets
  vector<TopJet>* topjets = event.topjets;

  // apply SFs
  if(isMC){
    lumiweight->process(event);
    PUreweight->process(event);
    sf_btag->process(event);
    muo_tight_noniso_SF->process(event);
    muo_trigger_SF->process(event);
    topreweight->process(event);
  }

  h_ecf->fill(event);
  h_pf->fill(event);
  h_EnergyFractions->fill(event);


  // store input particles
  // (do not use pointer here!)
  vector<PFParticle> pf_input = event.get(h_pfparticles);

  // FILL MJET HISTOGRAMS
  double pt = topjets->at(0).pt();
  double mass = topjets->at(0).softdropmass();
  double chf = topjets->at(0).chargedHadronEnergyFraction();
  h_mjet->fill(event);
  if(pt > 300 && pt < 400){
    h_mjet_pt300->fill(event);
    h_EnergyFractions_pt300->fill(event);
    if(chf < 0.6) h_mjet_loCHF_pt300->fill(event);
    if(chf > 0.6) h_mjet_hiCHF_pt300->fill(event);
  }
  else if(pt > 400){
    h_EnergyFractions_pt400->fill(event);
    h_mjet_pt400->fill(event);
    if(chf < 0.6) h_mjet_loCHF_pt400->fill(event);
    if(chf > 0.6) h_mjet_hiCHF_pt400->fill(event);
    if(fabs(mass-80) < 10){
      h_pf_W->fill(event);
      h_EnergyFractions_W->fill(event);
    }
    else if(fabs(mass-172) < 10){
      h_pf_top->fill(event);
      h_EnergyFractions_top->fill(event);
    }
  }

  if(masscut->passes(event)){
    h_mjet_masscut->fill(event);
    if(pt > 300 && pt < 400) h_mjet_masscut_pt300->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400->fill(event);
  }
  else{
    h_mjet_failmasscut->fill(event);
    if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400->fill(event);
  }

  // include some variation
  pf_var1->process(event);
  h_pf_var1->fill(event);
  h_mjet_var1->fill(event);
  if(pt > 300 && pt < 400) h_mjet_pt300_var1->fill(event);
  else if(pt > 400) h_mjet_pt400_var1->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var1->fill(event);
    if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var1->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var1->fill(event);
  }
  else{
    h_mjet_failmasscut_var1->fill(event);
    if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var1->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var1->fill(event);
  }

  event.set(h_pfparticles, pf_input);// set pf back to input values
  h_pf_aftervar1->fill(event);

  // different variation
  pf_var2->process(event);
  h_pf_var2->fill(event);
  h_mjet_var2->fill(event);
  if(pt > 300 && pt < 400) h_mjet_pt300_var2->fill(event);
  else if(pt > 400) h_mjet_pt400_var2->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var2->fill(event);
    if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var2->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var2->fill(event);
  }
  else{
    h_mjet_failmasscut_var2->fill(event);
    if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var2->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var2->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values

  // different variation
  pf_var3->process(event);
  h_pf_var3->fill(event);
  h_mjet_var3->fill(event);
  if(pt > 300 && pt < 400) h_mjet_pt300_var3->fill(event);
  else if(pt > 400) h_mjet_pt400_var3->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var3->fill(event);
    if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var3->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var3->fill(event);
  }
  else{
    h_mjet_failmasscut_var3->fill(event);
    if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var3->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var3->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values

  // different variation
  pf_var4->process(event);
  h_pf_var4->fill(event);
  h_mjet_var4->fill(event);
  if(pt > 300 && pt < 400) h_mjet_pt300_var4->fill(event);
  else if(pt > 400) h_mjet_pt400_var4->fill(event);

  if(masscut->passes(event)){
    h_mjet_masscut_var4->fill(event);
    if(pt > 300 && pt < 400) h_mjet_masscut_pt300_var4->fill(event);
    else if(pt > 400) h_mjet_masscut_pt400_var4->fill(event);
  }
  else{
    h_mjet_failmasscut_var4->fill(event);
    if(pt > 300 && pt < 400) h_mjet_failmasscut_pt300_var4->fill(event);
    else if(pt > 400) h_mjet_failmasscut_pt400_var4->fill(event);
  }
  event.set(h_pfparticles, pf_input);// set pf back to input values
  // DONE, don't store Analysis Tree
  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(TopMassModule)
