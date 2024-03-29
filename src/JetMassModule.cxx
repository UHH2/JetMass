#include <iostream>
#include <memory>


#ifdef CMSSW_10_6_28
#include "UHH2/JetMass/include/JetMassUtils106X.h"
#else
#include "UHH2/JetMass/include/JetMassUtils102X.h"
#endif

  
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h" //mainly for runPeriods vectors

#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MuonHists.h"

#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/JetMass/include/JetMassSelections.h"
#include "UHH2/JetMass/include/TopJetCorrections.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/common/include/MCLargeWeightKiller.h"
#include "UHH2/JetMass/include/WriteOutput.h"
#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/JetMass/include/JetMassGenHists.h"
#include "UHH2/JetMass/include/TriggerHists.h"
#include "UHH2/JetMass/include/JetMassUtils.h"

#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/PartonHT.h"
#include "UHH2/common/include/EventVariables.h"

#include <unistd.h>

#include "TFile.h"
#include "TH1F.h"

#define EXTRAOUT false


using namespace std;
using namespace uhh2;


class JetMassModule: public AnalysisModule {
public:

  explicit JetMassModule(Context & ctx);
  virtual bool process(Event & event) override;
  double get_N2DDT_value(double n2, double mass, double pt);

private:
  TString version;
  bool is_mc, is_QCD, matchV, is_WSample, is_ZSample, is_buggyPU, isUL;
  bool isttbarSel = false;
  bool isvjetsSel = false;
  bool do_genStudies = false;
  bool doOnlyTriggerHists = false;
  bool save_all_jets = false;
  int write_output_level = 0;

  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

  std::unique_ptr<N2DDTComputer<TopJet>> reco_n2ddt_computer;
  std::unique_ptr<N2DDTComputer<GenTopJet>> gen_n2ddt_computer;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> topPtReweighting;
  std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;
  std::unique_ptr<TopJetCorrections> topjetCorr,topjetCorr_chs;
  std::unique_ptr<StandaloneTopJetCorrector> topjet_jec_ev;
  std::unique_ptr<AnalysisModule>ttgen_producer;  
  std::unique_ptr<AnalysisModule> matching_selection_producer;
  std::unique_ptr<AnalysisModule> nlo_weights, nlo_weights_UL;
  std::unique_ptr<AnalysisModule> recojet_selector;
  std::unique_ptr<AnalysisModule> genjet_selector;
  std::unique_ptr<AnalysisModule> ht_calculator;
  
  std::unique_ptr<JetCleaner> ak4cleaner, ak4cleaner15;
  std::unique_ptr<TopJetCleaner> ak8cleaner,ak8cleaner_chs;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> ak8cleaner_dRlep,ak8cleaner_chs_dRlep;
  std::unique_ptr<GenTopJetCleaner> genak8cleaner;
  std::unique_ptr<GenJetCleaner> genak4cleaner;

  std::unique_ptr<AndSelection> reco_selection_part1,reco_selection_part2;
  std::unique_ptr<AndSelection> gen_selection_part1,gen_selection_part2;
  
  std::vector<std::unique_ptr<uhh2::Hists>> hists_nosel, hists;
  std::unique_ptr<uhh2::Hists> h_pfhists_200to500, h_pfhists_500to1000, h_pfhists_1000to2000, h_pfhists_2000to3000, h_pfhists_3000to4000, h_pfhists_4000to5000;
  std::unique_ptr<uhh2::Hists> h_pfhists_inclusive, h_pfhists_500to550, h_pfhists_550to600, h_pfhists_600to675, h_pfhists_675to800, h_pfhists_800to1200, h_pfhists_1200toInf;

  std::unique_ptr<uhh2::Hists> h_hlt_eff;
  std::unique_ptr<uhh2::Hists> h_gen_hists_commonmodules,h_gen_hists_gensel;
  
  std::unique_ptr<AnalysisModule> writer;


  uhh2::Event::Handle<bool>handle_reco_selection;
  uhh2::Event::Handle<bool>handle_gen_selection;

  uhh2::Event::Handle<float>handle_dR_reco_gen;
  uhh2::Event::Handle<bool>handle_jetpfid;

  uhh2::Event::Handle<bool>handle_reco_selection_no_merged_partons;
  uhh2::Event::Handle<bool>handle_gen_selection_no_merged_partons;

  uhh2::Event::Handle<int>handle_n_merged_partons_reco_jet;
  uhh2::Event::Handle<int>handle_n_merged_partons_gen_jet;
  
  uhh2::Event::Handle<double>handle_gen_HT,handle_HT;
  
  uhh2::Event::Handle<TTbarGen>handle_ttbar_gen;
  uhh2::Event::Handle<float> handle_toppt_weight;

  uhh2::Event::Handle<const GenTopJet*> handle_gentopjet,handle_gentopjet_nocut;
  uhh2::Event::Handle<const TopJet*> handle_recotopjet,handle_recotopjet_chs;
  uhh2::Event::Handle<double> handle_genjetpt, handle_V_pt;
  uhh2::Event::Handle<MatchingSelection> handle_matching_selection;

  uhh2::Event::Handle<std::vector<TopJet>> handle_chs_jets;
  uhh2::Event::Handle<std::vector<GenTopJet>> handle_sd_gen_jets;

  //systweights

  //JEC factors
  uhh2::Event::Handle<double> handle_jecfactor_nominal, handle_jecfactor_up, handle_jecfactor_down;

  //genjet values handles
  uhh2::Event::Handle<double> handle_gentopjet_pt_0,handle_gentopjet_msd_0,handle_gentopjet_mass_0,handle_gentopjet_dR_reco_0,handle_gentopjet_n2_0;
  uhh2::Event::Handle<double> handle_gentopjet_pt_1,handle_gentopjet_msd_1,handle_gentopjet_mass_1,handle_gentopjet_dR_reco_1,handle_gentopjet_n2_1;

  TopJetId jetpfid;

  
};


JetMassModule::JetMassModule(Context & ctx){

  std::cout << "Setting up JetMassModule with UHH2::" << UHH2_Release() << " release!" << std::endl;
  // Set some boolians
  version = ctx.get("dataset_version");
  is_mc = ctx.get("dataset_type") == "MC";
  isUL = is_UL(extract_year(ctx));
  is_QCD = ctx.get("dataset_version").find("QCD") != std::string::npos;
  is_WSample = ctx.get("dataset_version").find("WJets") != std::string::npos;
  is_ZSample = ctx.get("dataset_version").find("ZJets") != std::string::npos;
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;
  std::cout << "This sample has "<< (is_buggyPU ? "buggyPU" : "normalPU") << std::endl;

  save_all_jets =  string2bool(ctx.get("SaveAllJets","false"));

  write_output_level = (int) string2double(ctx.get("WriteOutputLevel","0"));

    
  const std::string& selection_ = ctx.get("selection", "");
  if     (selection_ == "ttbar") isttbarSel = true;
  else if(selection_ == "vjets")   isvjetsSel = true;
  else if(selection_ == "none") std::cout << "ATTENTION: You are about to run on a sample and will not apply any selection!" << std::endl;
  else throw runtime_error("JetMassModule: Select 'ttbar' or 'vjets' selection");

  std::string matching_selection_handlename("matching_selection");
  
  matching_selection_producer.reset(new MatchingSelectionProducer(ctx, matching_selection_handlename));
  handle_matching_selection = ctx.get_handle<MatchingSelection>(matching_selection_handlename);

  std::string genjetpt_handlename("genjetpt");
  nlo_weights.reset(new NLOWeights(ctx,genjetpt_handlename));
  std::string V_pt_handlename("V_pt");
  nlo_weights_UL.reset(new NLOWeights(ctx,V_pt_handlename,true));
  do_genStudies = string2bool(ctx.get("doGenStudies", "true"));
  doOnlyTriggerHists = string2bool(ctx.get("justDoTriggerHists", "false"));

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));

  jetpfid = JetPFID(JetPFID::WP_TIGHT_PUPPI);

  if(is_buggyPU){
    TString pu_file_path = (TString) ctx.get("dataset_version");
    pu_file_path = pu_file_path.ReplaceAll("_buggyPU","");
    pu_file_path = pu_file_path.ReplaceAll("_test","");
    pu_file_path = "common/data/2017/Pileup_QCD_PtBinned/MyMCPileupHistogram_"+pu_file_path+".root";
    ctx.set("pileup_directory",(std::string) pu_file_path);
  }
  std::cout << "reweighting mc pileup using " << ctx.get("pileup_directory")<<" as mc profile dir" <<std::endl;
  ctx.undeclare_all_event_output();

  common.reset(new CommonModules());
  //muid and eleid are taken from release dependent header file (JetMassUtils10[2,6]X.h)
  if(isttbarSel && !doOnlyTriggerHists){
    common->set_muon_id(muid);
    common->set_electron_id(eleid);
  }else if(isvjetsSel){
    common->set_muon_id(muid_loose);
    common->set_electron_id(eleid_loose);
  }
  common->set_jet_id(jetid);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection(true);
  common->switch_jetPtSorter();
  common->set_HTjetid(jetid);
  
  common->init(ctx);

  if(is_mc && isvjetsSel && is_QCD){
    mcSpikeKiller.reset(new MCLargeWeightKiller(ctx,
                                                2, // maximum allowed ratio of leading reco jet pT / generator HT
                                                2, // maximum allowed ratio of leading gen jet pT / generator HT
                                                2, // maximum allowed ratio of leading reco jet pT / Q scale
                                                2, // maximum allowed ratio of PU maximum pTHat / gen HT (ensures scale of PU < scale of hard interaction)
                                                2, // maximum allowed ratio of leading reco jet pT / pTHat
                                                2, // maximum allowed ratio of leading gen jet pT / pTHat
                                                "jets", // name of jet collection to be used
                                                "genjets" // name of genjet collection to be used
                          ));
  }

  //HandleNames
  std::string reco_selection_handlename("pass_reco_selection");
  std::string gen_selection_handlename("pass_gen_selection");

  std::string dR_reco_gen_handlename("dR_reco_gen");
  std::string jetpfid_handlename("jetpfid");
  
  std::string n_merged_partons_reco_jet_handlename("n_merged_partons_reco_jet");
  std::string n_merged_partons_gen_jet_handlename("n_merged_partons_gen_jet");

  std::string ttbargen_handlename("ttbar_gen_system");

  std::string genHT_handlename("genHT");
  std::string HT_handlename("HT");

  std::string gentopjet_handlename("gentopjet");
  std::string gentopjet_nocut_handlename("gentopjet_nocut");
  std::string recotopjet_handlename("recotopjet");

  std::string recotopjet_chs_handlename("recotopjet_chs");

  //Handles
  std::string chs_jets_collection = ctx.get("chsjets");
  handle_chs_jets = ctx.get_handle<std::vector<TopJet>>(chs_jets_collection);

  std::string sd_gen_jets_collection = ctx.get("sdgenjets");
  // handle_sd_gen_jets = ctx.get_handle<std::vector<GenTopJet>>(sd_gen_jets_collection);
  handle_sd_gen_jets = ctx.get_handle<std::vector<GenTopJet>>(sd_gen_jets_collection);

  
  handle_reco_selection = ctx.declare_event_output<bool>(reco_selection_handlename);
  handle_gen_selection = ctx.declare_event_output<bool>(gen_selection_handlename);

  handle_dR_reco_gen = ctx.declare_event_output<float>(dR_reco_gen_handlename);
  handle_jetpfid = ctx.declare_event_output<bool>(jetpfid_handlename);
    
  handle_reco_selection_no_merged_partons = ctx.declare_event_output<bool>(reco_selection_handlename+"_no_merged_partons");
  handle_gen_selection_no_merged_partons = ctx.declare_event_output<bool>(gen_selection_handlename+"_no_merged_partons");

  handle_n_merged_partons_reco_jet = ctx.declare_event_output<int>(n_merged_partons_reco_jet_handlename);
  handle_n_merged_partons_gen_jet = ctx.declare_event_output<int>(n_merged_partons_gen_jet_handlename);

  
  handle_gentopjet = ctx.get_handle<const GenTopJet*>(gentopjet_handlename);
  handle_gentopjet_nocut = ctx.get_handle<const GenTopJet*>(gentopjet_nocut_handlename);
  handle_recotopjet = ctx.get_handle<const TopJet*>(recotopjet_handlename);
  handle_recotopjet_chs = ctx.get_handle<const TopJet*>(recotopjet_chs_handlename);

  handle_V_pt = ctx.declare_event_output<double>(V_pt_handlename);
  handle_genjetpt = ctx.declare_event_output<double>(genjetpt_handlename);
  handle_gen_HT = ctx.declare_event_output<double>(genHT_handlename);
  handle_HT = ctx.declare_event_output<double>(HT_handlename);
  ht_calculator.reset(new HTCalculator(ctx,jetid,HT_handlename));

  handle_jecfactor_nominal = ctx.declare_event_output<double>("jecfactor_nominal");
  handle_jecfactor_up = ctx.declare_event_output<double>("jecfactor_up");
  handle_jecfactor_down = ctx.declare_event_output<double>("jecfactor_down");

  // save gentopjets
  handle_gentopjet_pt_0 = ctx.declare_event_output<double>("gentopjet_pt_0");
  handle_gentopjet_msd_0 = ctx.declare_event_output<double>("gentopjet_msd_0");
  handle_gentopjet_mass_0 = ctx.declare_event_output<double>("gentopjet_mass_0");
  handle_gentopjet_dR_reco_0 = ctx.declare_event_output<double>("gentopjet_dR_reco_0");
  handle_gentopjet_n2_0 = ctx.declare_event_output<double>("gentopjet_n2_0");

  handle_gentopjet_pt_1 = ctx.declare_event_output<double>("gentopjet_pt_1");
  handle_gentopjet_msd_1 = ctx.declare_event_output<double>("gentopjet_msd_1");
  handle_gentopjet_mass_1 = ctx.declare_event_output<double>("gentopjet_mass_1");
  handle_gentopjet_dR_reco_1 = ctx.declare_event_output<double>("gentopjet_dR_reco_1");
  handle_gentopjet_n2_1 = ctx.declare_event_output<double>("gentopjet_n2_1");

  
  //AnalysisModules
  ttgen_producer.reset(new TTbarGenProducer(ctx,ttbargen_handlename,false));

  handle_toppt_weight = ctx.declare_event_output<float>("toppt_weight");
  topPtReweighting.reset(new TopPtReweight(ctx, 0.0615, -0.0005,ttbargen_handlename,"toppt_weight",false));

  recojet_selector.reset(new JetSelector<TopJet>(ctx,recotopjet_handlename, matching_selection_handlename));
  genjet_selector.reset(new JetSelector<GenTopJet>(ctx,gentopjet_handlename));

  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);
  topjetCorr->disable_jersmear();
  topjet_jec_ev.reset(new StandaloneTopJetCorrector(ctx, "AK8PFPuppi"));

  // AK8 CHS JEC/JER
  topjetCorr_chs.reset(new TopJetCorrections(chs_jets_collection));
  topjetCorr_chs->init(ctx);
  topjetCorr_chs->disable_jersmear();
  
  
  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;

  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  ak8cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(ctx, 0.8));
  ak8cleaner_chs.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta)),chs_jets_collection));
  ak8cleaner_chs_dRlep.reset(new TopJetLeptonDeltaRCleaner(ctx, 0.8,chs_jets_collection));

  genak4cleaner.reset(new GenJetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  genak8cleaner.reset(new GenTopJetCleaner(ctx,GenTopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  
  
  //Selections
  reco_selection_part1.reset(new AndSelection(ctx,"reco_selection_part1"));
  reco_selection_part2.reset(new AndSelection(ctx,"reco_selection_part2"));
  if(isttbarSel){
    OrSelection trigger_sel =  OrSelection();
    trigger_sel.add<TriggerSelection>("HLT_Mu50_v*");
    if(extract_year(ctx) == Year::isUL16preVFP || extract_year(ctx) == Year::isUL16postVFP){
      trigger_sel.add<TriggerSelection>("HLT_TkMu50_v*", true);
    }
    reco_selection_part1->add<OrSelection>("Trigger selection", trigger_sel);
    reco_selection_part1->add<NTopJetSelection>("N_{AK8} #geq 1", 1);
    reco_selection_part1->add<TopJetPtCut>("p_{T, AK8} #geq 200 GeV", ctx, 200);
    reco_selection_part1->add<NMuonSelection>("N_{#mu} #geq 1", 1,1);
    reco_selection_part1->add<NElectronSelection>("N_{e} = 0", 0,0);
    reco_selection_part1->add<METCut>("MET > 50 GeV", 50.,100000.);
    reco_selection_part1->add<TwoDCut>("2D-Cut",0.4,25.);
    reco_selection_part1->add<NMuonBTagSelection>("b-jet in muon hemisphere", 1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM));

  }else if(isvjetsSel){
    reco_selection_part1->add<NElectronSelection>("ele-veto",0,0);
    reco_selection_part1->add<NMuonSelection>("muon-veto",0,0);
    reco_selection_part1->add<NTopJetSelection>("N_{AK8} #geq 1", 1);
    reco_selection_part1->add<TopJetPtCut>("p_{T, AK8} #geq 500 GeV", ctx, 500);
    reco_selection_part2->add<RhoCut<TopJet>>("rhocut",ctx, -6.0, -2.1,recotopjet_handlename);
  }


  gen_selection_part1.reset(new AndSelection(ctx,"gen_selection_part1"));
  gen_selection_part2.reset(new AndSelection(ctx,"gen_selection_part2"));
  if(isttbarSel){
    gen_selection_part1->add<NGenTopJetSelection>("N_{gen,AK8} #geq 1, p_{T} > 200 GeV", 1,-1,GenTopJetId(PtEtaCut(200.,100000.)));
    gen_selection_part1->add<TTbarGenSemilepSelection>("gen semilep-selection",ctx, ttbargen_handlename,55.);
    gen_selection_part1->add<METCut>("MET > 50 GeV", 50.,100000.,true);
  }else if(isvjetsSel){
    gen_selection_part1->add<NGenTopJetSelection>("N_{gen,AK8} #geq 1, p_{T} > 500 GeV", 1,-1,GenTopJetId(PtEtaCut(500.,100000.)));
  }
  
  
  // HISTOGRAMS
  hists.emplace_back(new ElectronHists(ctx, "ElectronHists"));
  hists.emplace_back(new EventHists(ctx, "EventHists"));
  hists.emplace_back(new MuonHists(ctx, "MuonHists"));
  hists.emplace_back(new JetHists(ctx, "JetHists"));
  hists.emplace_back(new TopJetHists(ctx, "TopJetHists"));

  hists_nosel.emplace_back(new ElectronHists(ctx, "ElectronHists_nosel"));
  hists_nosel.emplace_back(new EventHists(ctx, "EventHists_nosel"));
  hists_nosel.emplace_back(new MuonHists(ctx, "MuonHists_nosel"));
  hists_nosel.emplace_back(new JetHists(ctx, "JetHists_nosel"));
  hists_nosel.emplace_back(new TopJetHists(ctx, "TopJetHists_nosel"));

  
  h_hlt_eff.reset(new TriggerHists(ctx, "HLTEffHists"));

  h_gen_hists_commonmodules.reset(new JetMassGenHists(ctx,"GenHistsCommonModules",ttbargen_handlename,genHT_handlename));
  h_gen_hists_gensel.reset(new JetMassGenHists(ctx, "GenHistsGenSel",ttbargen_handlename,genHT_handlename));

  writer.reset(new WriteOutput(ctx,matching_selection_handlename,recotopjet_handlename,recotopjet_chs_handlename,gentopjet_handlename,write_output_level));

  std::string N2DDT_file_path = "JetMass/Histograms/ddtmaps/QCD_2017_PFMass_smooth_gaus4p00sigma.root";
  std::string N2DDT_hist_name = "N2_v_pT_v_rho_0p05_smooth_gaus4p00sigma_maps_cleaner_PFMass";
  reco_n2ddt_computer.reset(new N2DDTComputer<TopJet>(N2DDT_file_path,N2DDT_hist_name));
  gen_n2ddt_computer.reset(new N2DDTComputer<GenTopJet>(N2DDT_file_path,N2DDT_hist_name));
}

bool JetMassModule::process(Event & event) {
  bool extra_out(false); // for rare debug cases where one want to decide for verbose output on Runtime level
  if(EXTRAOUT || extra_out){
    cout << "JetMassModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }
  if(is_mc)ttgen_producer->process(event);

  matching_selection_producer->process(event);
  MatchingSelection matching_selection = event.get(handle_matching_selection);
  
  double genHT(-1.0);
  if(is_mc)genHT = calculateGenJetPtSum(event);
  event.set(handle_gen_HT,genHT);
  ht_calculator->process(event);
  
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: ExtraHandles done!" << std::endl;
  
  // Throw away events with NVTX in buggy area for buggy samples
  if(is_mc && is_buggyPU){
    float n_true = event.genInfo->pileup_TrueNumInteractions();
    if(n_true < 10. || n_true > 72.) return false;
    if(EXTRAOUT || extra_out) std::cout << "JetMassModule: BuggyPU Treatment done!" << std::endl;
  }

  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;

  // apply prefiring weight
  event.weight *= event.prefiringWeight;
  
  h_gen_hists_commonmodules->fill(event);
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: CommonModules done!" << std::endl;

  //remove MC Events with very large unphysical weights
  if(is_mc && isvjetsSel && is_QCD){
    if(!mcSpikeKiller->passes(event)) return false;
    if(EXTRAOUT || extra_out) std::cout << "JetMassModule: SpikeKiller done!" << std::endl;
  }

  event.set(handle_toppt_weight, 1.0);
  if(is_mc)topPtReweighting->process(event);

  std::vector<TopJet> & chs_jets = event.get(handle_chs_jets);
  // std::vector<GenTopJet> & sd_gen_jets = event.get(handle_sd_gen_jets);
  std::vector<GenTopJet> & sd_gen_jets = *(new std::vector<GenTopJet>({}));
  if(is_mc && event.is_valid(handle_sd_gen_jets))sd_gen_jets = event.get(handle_sd_gen_jets);

// AK8 JEC
  topjetCorr->process(event);
  topjetCorr_chs->process(event);

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);
  sort_by_pt<TopJet>(chs_jets);
  if(is_mc)sort_by_pt<GenTopJet>(*event.gentopjets);
  if(is_mc)sort_by_pt<GenTopJet>(sd_gen_jets);  
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: TopJetCorrections done!" << std::endl;
  if(EXTRAOUT || extra_out) std::cout << "N topjets after topjetcorrections: " << event.topjets->size() << std::endl;
  
  // CLEANER
  ak4cleaner->process(event);
  ak8cleaner->process(event);
  ak8cleaner_dRlep->process(event);

  ak8cleaner_chs->process(event);
  ak8cleaner_chs_dRlep->process(event);
  if(is_mc){
    genak8cleaner->process(event);
    genak4cleaner->process(event);
  }
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: Cleaner done!" << std::endl;
  if(EXTRAOUT || extra_out) std::cout << "N topjets after cleaner: " << event.topjets->size() << std::endl;

  double genjetpt(-1.0), bosonpt(-1.0);
  if(is_mc){
    if(isvjetsSel && (is_WSample || is_ZSample)){
      const Particle gen_V_particle = event.genparticles->at(get_V_index(*event.genparticles));
      bosonpt = gen_V_particle.pt();
      if(event.topjets->size()>0){
        //get genjet pt for k factors
        const GenJet * closest_genjet_1 = closestParticle(event.topjets->at(0), *event.genjets);
        const GenJet * closest_genjet_2 = event.topjets->size() > 1 ? closestParticle(event.topjets->at(1), *event.genjets) : closest_genjet_1;
        float gen_pt_1 = closest_genjet_1 ? closest_genjet_1->pt() : -999;
        float gen_pt_2 = closest_genjet_2 ? closest_genjet_2->pt() : -999;
        genjetpt = matching_selection.passes_matching(event.topjets->at(0), MatchingSelection::oIsMergedV)? gen_pt_1 : gen_pt_2;
      }
    }
  }
  event.set(handle_V_pt, bosonpt);
  event.set(handle_genjetpt, genjetpt);

  if(is_mc){
    if(isUL) nlo_weights_UL->process(event);
    else nlo_weights->process(event);
  }
  h_hlt_eff->fill(event);
  if(doOnlyTriggerHists) return false;

  float jecfactor_standalone_nominal(-1.),jecfactor_standalone_up(-1.),jecfactor_standalone_down(-1.);
  if(event.topjets->size() >0){
    jecfactor_standalone_nominal = topjet_jec_ev->getJecFactor(event, event.topjets->at(0), 0);
    jecfactor_standalone_up = topjet_jec_ev->getJecFactor(event, event.topjets->at(0), +1);
    jecfactor_standalone_down = topjet_jec_ev->getJecFactor(event, event.topjets->at(0), -1);
  }
  event.set(handle_jecfactor_nominal, jecfactor_standalone_nominal);
  event.set(handle_jecfactor_up, jecfactor_standalone_up);
  event.set(handle_jecfactor_down, jecfactor_standalone_down);
  
  bool pass_reco_selection_part1 = reco_selection_part1->passes(event);
  bool pass_reco_selection_part2(false);
  //selecting reco puppi jet and setting respective handle // TODO: make sure selection is applies correctly since i now also set genjet handle with jets that do not explicitily pass selection!!!
  recojet_selector->process(event);
  if(pass_reco_selection_part1){
    
    pass_reco_selection_part2 = reco_selection_part2->passes(event);
  }
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: RecoSelection done!" << std::endl;


  bool pass_gen_selection_part1(false),pass_gen_selection_part2(false);
  if(is_mc)pass_gen_selection_part1 = gen_selection_part1->passes(event);
  else pass_gen_selection_part1 = true;
  
  float gentopjet_pt_0(-1.0);
  float gentopjet_msd_0(-1.0);
  float gentopjet_mass_0(-1.0);
  float gentopjet_dR_0(-1.0);
  float gentopjet_n2_0(-1.0);
  float gentopjet_pt_1(-1.0);
  float gentopjet_msd_1(-1.0);
  float gentopjet_mass_1(-1.0);
  float gentopjet_dR_1(-1.0);
  float gentopjet_n2_1(-1.0);
  
  //putting leading genjet into handle used for Unfolding stuff // TODO: make sure selection is applies correctly since i now also set genjet handle with jets that do not explicitily pass selection!!!
  if(is_mc){

    if(event.gentopjets->size()>0){
        gentopjet_pt_0 = event.gentopjets->at(0).pt();
        gentopjet_mass_0 = event.gentopjets->at(0).v4().M();
        const TopJet* recojet(NULL);
        if(event.is_valid(handle_recotopjet))recojet = event.get(handle_recotopjet);
        if(recojet) gentopjet_dR_0 = deltaR(event.gentopjets->at(0),*recojet);
        if(sd_gen_jets.size()>0){
          const GenTopJet * groomed = closestParticle(event.gentopjets->at(0),sd_gen_jets);
          if(deltaR(*groomed,event.gentopjets->at(0)) < 0.4){
            gentopjet_n2_0 = groomed->ecfN2_beta1();
            gentopjet_msd_0 = groomed->softdropmass();            
          }
        }
      }
    if(event.gentopjets->size()>1){
        gentopjet_pt_1 = event.gentopjets->at(1).pt();
        gentopjet_mass_1 = event.gentopjets->at(1).v4().M();
        const TopJet* recojet(NULL);
        if(event.is_valid(handle_recotopjet))recojet = event.get(handle_recotopjet);
        if(recojet) gentopjet_dR_1 = deltaR(event.gentopjets->at(1),*recojet);
        if(sd_gen_jets.size()>0){
          const GenTopJet * groomed = closestParticle(event.gentopjets->at(1),sd_gen_jets);
          if(deltaR(*groomed,event.gentopjets->at(1)) < 0.4){
            gentopjet_n2_1 = groomed->ecfN2_beta1();
            gentopjet_msd_1 = groomed->softdropmass();
          }
        }
      }
    
    genjet_selector->process(event);
        
    const GenTopJet * ungroomed_gentopjet(NULL);
    if(event.is_valid(handle_gentopjet))ungroomed_gentopjet = event.get(handle_gentopjet);
    GenTopJet * ungroomed_gentopjet_new = new GenTopJet();
    if(ungroomed_gentopjet && sd_gen_jets.size()>0){
      ungroomed_gentopjet_new->set_v4(ungroomed_gentopjet->v4());
      ungroomed_gentopjet_new->set_tau1(ungroomed_gentopjet->tau1());
      ungroomed_gentopjet_new->set_tau2(ungroomed_gentopjet->tau2());
      ungroomed_gentopjet_new->set_tau3(ungroomed_gentopjet->tau3());
      ungroomed_gentopjet_new->set_tau4(ungroomed_gentopjet->tau4());
      const GenTopJet * groomed_gentopjet = closestParticle(*ungroomed_gentopjet_new,sd_gen_jets);
      ungroomed_gentopjet_new->set_ecfN2_beta1(groomed_gentopjet->ecfN2_beta1());
      ungroomed_gentopjet_new->set_ecfN2_beta2(groomed_gentopjet->ecfN2_beta2());
      ungroomed_gentopjet_new->set_ecfN3_beta1(groomed_gentopjet->ecfN3_beta1());
      ungroomed_gentopjet_new->set_ecfN3_beta2(groomed_gentopjet->ecfN3_beta2());
      
      for(unsigned int isj=0; isj < groomed_gentopjet->subjets().size(); isj++) {
        GenJet subjet = groomed_gentopjet->subjets().at(isj); 
        ungroomed_gentopjet_new->add_subjet(subjet);
      }
     
      event.set(handle_gentopjet, ungroomed_gentopjet_new);
    }
  }

  event.set(handle_gentopjet_pt_0,gentopjet_pt_0);
  event.set(handle_gentopjet_msd_0,gentopjet_msd_0);
  event.set(handle_gentopjet_mass_0,gentopjet_mass_0);
  event.set(handle_gentopjet_dR_reco_0,gentopjet_dR_0);
  event.set(handle_gentopjet_n2_0,gentopjet_n2_0);
  event.set(handle_gentopjet_pt_1,gentopjet_pt_1);
  event.set(handle_gentopjet_msd_1,gentopjet_msd_1);
  event.set(handle_gentopjet_mass_1,gentopjet_mass_1);
  event.set(handle_gentopjet_dR_reco_1,gentopjet_dR_1);
  event.set(handle_gentopjet_n2_1,gentopjet_n2_1);
  
  if(pass_gen_selection_part1 && is_mc){
    pass_gen_selection_part2 = gen_selection_part2->passes(event);
  }
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: GenSelection done!" << std::endl;

  //putting closest genjet without selection and chs-jet into handle used for CHSPuppiResponseStudy
  const GenTopJet* genjet(NULL);
  const TopJet * chs_jet(NULL);
  const TopJet* recojet(NULL);
  if(is_mc){
    if(event.is_valid(handle_recotopjet)) recojet=event.get(handle_recotopjet);
    if(recojet){
      //matching chs jet to reco puppi jet and setting respective handle
      chs_jet = closestParticle(*recojet, chs_jets);

      // matching ungroomed gen jet to reco puppi jet
      // const_cast should be highly illegal, but unfortunately closestParticle returns a const GenTopJet * so i would not be able to copy subjets to this - so doing it this hacky way
      GenTopJet * genjet_ungroomed = const_cast<GenTopJet*>(closestParticle(*recojet, *event.gentopjets));
      //matching groomed gen jet to ungroomed one if there are any groomed jets and copy subjets to ungroomed jet
      if(genjet_ungroomed && sd_gen_jets.size()>0){
        const GenTopJet* sd_genjet = closestParticle(*genjet_ungroomed, sd_gen_jets);
        for(unsigned int isj=0; isj<sd_genjet->subjets().size(); isj++) {
          GenJet subjet = sd_genjet->subjets().at(isj); 
          genjet_ungroomed->add_subjet(subjet);
        }
      }
      genjet = genjet_ungroomed;
    }
  }

  event.set(handle_gentopjet_nocut,genjet);
  event.set(handle_recotopjet_chs,chs_jet);

  
  //set reco selection to fail if topjet-gentopjet matching within R/2=0.4 fails.  
  bool pass_dR_reco_gen(true),pass_gen_substructure(false);
  float dR_reco_gen(-1.0);
  if( (pass_gen_selection_part1 && pass_reco_selection_part1 )){
    double dR(numeric_limits<double>::infinity());
    if(event.is_valid(handle_gentopjet) && event.is_valid(handle_recotopjet)){
      const GenTopJet* genjet = event.get(handle_gentopjet);
      const TopJet* recojet = event.get(handle_recotopjet);
      dR_reco_gen = deltaR(*genjet,*recojet);
      pass_gen_substructure = (safe_tau21(genjet) < 0.45) & (safe_tau32(genjet)>0.5);
      
    }
    if(dR > 0.4) pass_dR_reco_gen = false;
  }
  event.set(handle_dR_reco_gen,dR_reco_gen);
  
   
  bool passing_jetpfid(true);
 
  int n_merged_partons_reco_jet(-1),n_merged_partons_gen_jet(-1);
  if(event.is_valid(handle_recotopjet)){
    const TopJet* recojet = event.get(handle_recotopjet);
    n_merged_partons_reco_jet = matching_selection.n_merged_partons(*recojet);
    passing_jetpfid = jetpfid(*recojet,event);
  }

  event.set(handle_jetpfid, passing_jetpfid);
  
  if(event.is_valid(handle_gentopjet)){
    const GenTopJet* genjet = event.get(handle_gentopjet);
    n_merged_partons_gen_jet = matching_selection.n_merged_partons(*genjet);
  }
  event.set(handle_n_merged_partons_reco_jet, n_merged_partons_reco_jet);
  event.set(handle_n_merged_partons_gen_jet, n_merged_partons_gen_jet);


  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: Reco-Gen Matching done!" << std::endl;

  bool pass_substructure_cut(false);  
  if(pass_reco_selection_part1){
    const TopJet * recotopjet = NULL;
    if(event.is_valid(handle_recotopjet))recotopjet=event.get(handle_recotopjet);
    if(isvjetsSel){
      double n2ddt = reco_n2ddt_computer->computeDDTValue(recotopjet);
      pass_substructure_cut = n2ddt < 0;
    }
    if(isttbarSel){
      double tau32 = safe_tau32(recotopjet);
      pass_substructure_cut = tau32<0.5;
    }
  }
  if(EXTRAOUT || extra_out) std::cout << "JetMassModule: SubstructureSelection done!" << std::endl;

  
  //fill second round of unfolding hists with part2 of selections
  bool pass_gen_selection = pass_gen_selection_part1; // keep part 2 out since it is empty at the moment anyway! && pass_gen_selection_part2;

  event.set(handle_reco_selection, pass_reco_selection_part1 );

  event.set(handle_gen_selection, pass_gen_selection);// && pass_gen_substructure);

  event.set(handle_reco_selection_no_merged_partons, pass_reco_selection_part1 && pass_substructure_cut && pass_dR_reco_gen && (n_merged_partons_reco_jet==0));
  event.set(handle_gen_selection_no_merged_partons, pass_gen_selection && (n_merged_partons_gen_jet==0));

  event.set(handle_reco_selection_no_merged_partons, pass_reco_selection_part1 && pass_substructure_cut && pass_reco_selection_part2 && pass_dR_reco_gen && (n_merged_partons_reco_jet==0));
  //Write everything used for JetMassCalibration to Tree if first part of reco-selection passes
  for(auto & h: hists_nosel) h->fill(event);
  if(
    (event.topjets->size() >0 && save_all_jets) ||
    pass_reco_selection_part1
    ){  
    // FILL HISTS
    for(auto & h: hists) h->fill(event);
  }
  if(
    (event.topjets->size() >0 && save_all_jets) ||
    ((pass_reco_selection_part1) ||  (pass_gen_selection && is_WSample && isvjetsSel))
    ){  
    // STORE EVENT
    writer->process(event);
    
    if(EXTRAOUT || extra_out) std::cout << "JetMassModule: OutputTreeWriter done!" << std::endl;
    return true;
  }else return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the JetMassModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(JetMassModule)
