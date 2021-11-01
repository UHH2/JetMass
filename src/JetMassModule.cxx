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
#include "UHH2/JetMass/include/UnfoldingHists.h"
#include "UHH2/JetMass/include/JetMassGenHists.h"
#include "UHH2/JetMass/include/JetMassUtils.h"

#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/PartonHT.h"

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
  bool is_mc, is_QCD, matchV, is_WSample, is_ZSample, is_buggyPU;
  bool isTTbarSel = false;
  bool isVJetsSel = false;
  bool do_genStudies = false;

  Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;

  std::unique_ptr<N2DDTComputer<TopJet>> reco_n2ddt_computer;
  std::unique_ptr<N2DDTComputer<GenTopJet>> gen_n2ddt_computer;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> topPtReweighting;
  std::unique_ptr<MCLargeWeightKiller> mcSpikeKiller;
  std::unique_ptr<TopJetCorrections> topjetCorr;
  std::unique_ptr<AnalysisModule>ttgen_producer;  
  std::unique_ptr<AnalysisModule> matching_selection_producer;
  std::unique_ptr<AnalysisModule> nlo_weights;
  std::unique_ptr<AnalysisModule> recojet_selector;
  std::unique_ptr<AnalysisModule> genjet_selector;

  std::unique_ptr<JetCleaner> ak4cleaner, ak4cleaner15;
  std::unique_ptr<TopJetCleaner> ak8cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> ak8cleaner_dRlep;
  std::unique_ptr<GenTopJetCleaner> genak8cleaner;
  std::unique_ptr<GenJetCleaner> genak4cleaner;

  std::unique_ptr<AndSelection> reco_selection_part1,reco_selection_part2;
  std::unique_ptr<AndSelection> gen_selection_part1,gen_selection_part2;
  
  std::vector<std::unique_ptr<uhh2::Hists>> hists;
  std::unique_ptr<uhh2::Hists> h_pfhists_200to500, h_pfhists_500to1000, h_pfhists_1000to2000, h_pfhists_2000to3000, h_pfhists_3000to4000, h_pfhists_4000to5000;
  std::unique_ptr<uhh2::Hists> h_pfhists_inclusive, h_pfhists_500to550, h_pfhists_550to600, h_pfhists_600to675, h_pfhists_675to800, h_pfhists_800to1200, h_pfhists_1200toInf;

  std::unique_ptr<uhh2::Hists> h_gen_hists_commonmodules,h_gen_hists_gensel;
  std::unique_ptr<uhh2::Hists> h_unfolding_hists, h_unfolding_hists_fine;
  std::unique_ptr<uhh2::Hists> h_unfolding_hists_sel_part1, h_unfolding_hists_sel_part2;
  
  std::unique_ptr<AnalysisModule> writer;


  uhh2::Event::Handle<bool>handle_reco_selection;
  uhh2::Event::Handle<bool>handle_gen_selection;

  uhh2::Event::Handle<double>handle_gen_HT;
  uhh2::Event::Handle<TTbarGen>handle_ttbar_gen;

  uhh2::Event::Handle<const GenTopJet*> handle_gentopjet;
  uhh2::Event::Handle<const TopJet*> handle_recotopjet;
  uhh2::Event::Handle<const GenTopJet*> handle_gentopjet_matched_V;
  uhh2::Event::Handle<MatchingSelection> handle_matching_selection;  


  
};


JetMassModule::JetMassModule(Context & ctx){

  std::cout << "Setting up JetMassModule with UHH2::" << UHH2_Release() << " release!" << std::endl;
  
  // Set some boolians
  version = ctx.get("dataset_version");
  is_mc = ctx.get("dataset_type") == "MC";
  is_QCD = ctx.get("dataset_version").find("QCD") != std::string::npos;
  is_buggyPU = ctx.get("dataset_version").find("buggyPU") != std::string::npos;
  std::cout << "This sample has "<< (is_buggyPU ? "buggyPU" : "normalPU") << std::endl;
  const std::string& selection_ = ctx.get("selection", "");
  if     (selection_ == "ttbar") isTTbarSel = true;
  else if(selection_ == "vjets")   isVJetsSel = true;
  else throw runtime_error("JetMassModule: Select 'ttbar' or 'vjets' selection");

  std::string matching_selection_handlename("matching_selection");
  
  matching_selection_producer.reset(new MatchingSelectionProducer(ctx, matching_selection_handlename));
  handle_matching_selection = ctx.get_handle<MatchingSelection>(matching_selection_handlename);

  std::string gentopjet_matched_V_handlename("gentopjet_matched_V");
  nlo_weights.reset(new NLOWeights(ctx,gentopjet_matched_V_handlename));
  
  do_genStudies = string2bool(ctx.get("doGenStudies", "true"));
  
  // common modules
  // MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
  // ElectronId eleid = AndId<Electron>(ElectronID_Fall17_medium_noIso, PtEtaCut(55., 2.4));
  // // ElectronId eleid = AndId<Electron>(ElectronTagID(Electron::tagname2tag("cutBasedElectronID-Fall17-94X-V2-medium")), PtEtaCut(55., 2.4));

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));


  if(is_buggyPU){
    TString pu_file_path = (TString) ctx.get("dataset_version");
    pu_file_path = pu_file_path.ReplaceAll("MC","");
    pu_file_path = pu_file_path.ReplaceAll("_buggyPU","");
    pu_file_path = pu_file_path.ReplaceAll("_test","");
    pu_file_path = "common/data/2017/Pileup_QCD_PtBinned/MyMCPileupHistogram"+pu_file_path+".root";
    ctx.set("pileup_directory",(std::string) pu_file_path);
  }
  std::cout << "reweighting mc pileup using " << ctx.get("pileup_directory")<<" as mc profile dir" <<std::endl;

  common.reset(new CommonModules());
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->set_jet_id(jetid);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection(true);
  common->switch_jetPtSorter();
  common->set_HTjetid(jetid);
  if(is_mc) common->disable_metfilters();
  common->init(ctx);

  if(is_mc && isVJetsSel && is_QCD){
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


  ctx.undeclare_all_event_output();

  //HandleNames
  std::string reco_selection_handlename("pass_reco_selection");
  std::string gen_selection_handlename("pass_gen_selection");

  std::string ttbargen_handlename("ttbar_gen_system");

  std::string genHT_handlename("genHT");

  std::string gentopjet_handlename("gentopjet");
  std::string recotopjet_handlename("recotopjet");
  
  //Handles

  handle_reco_selection = ctx.declare_event_output<bool>(reco_selection_handlename);
  handle_gen_selection = ctx.declare_event_output<bool>(gen_selection_handlename);

  handle_gentopjet = ctx.get_handle<const GenTopJet*>(gentopjet_handlename);
  handle_recotopjet = ctx.get_handle<const TopJet*>(recotopjet_handlename);
  handle_gentopjet_matched_V = ctx.get_handle<const GenTopJet*>(gentopjet_matched_V_handlename);
  
  handle_gen_HT = ctx.declare_event_output<double>(genHT_handlename);
  
  //AnalysisModules
  ttgen_producer.reset(new TTbarGenProducer(ctx,ttbargen_handlename,false));
  
  topPtReweighting.reset(new TopPtReweight(ctx, 0.0615, -0.0005,"","",true));

  recojet_selector.reset(new JetSelector<TopJet>(ctx,recotopjet_handlename));
  genjet_selector.reset(new JetSelector<GenTopJet>(ctx,gentopjet_handlename));

  // AK8 JEC/JER
  topjetCorr.reset(new TopJetCorrections());
  topjetCorr->init(ctx);
  topjetCorr->disable_jersmear();

  // Jet cleaner
  AK4_Clean_pT = 30.0;
  AK4_Clean_eta = 2.4;
  AK8_Clean_pT = 170.0;
  AK8_Clean_eta = 2.4;

  ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  ak8cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(0.8));

  genak4cleaner.reset(new GenJetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
  genak8cleaner.reset(new GenTopJetCleaner(ctx,GenTopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));
  
  
  //Selections
  reco_selection_part1.reset(new AndSelection(ctx,"reco_selection_part1"));
  reco_selection_part2.reset(new AndSelection(ctx,"reco_selection_part2"));
  if(isTTbarSel){
    reco_selection_part1->add<TriggerSelection>("Trigger selection","HLT_Mu50_v*");
    reco_selection_part1->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 200 GeV", 1,-1,TopJetId(PtEtaCut(200.,100000.)));
    reco_selection_part1->add<NMuonSelection>("N_{#mu} #geq 1", 1,1);
    reco_selection_part1->add<NElectronSelection>("N_{e} = 0", 0,0);
    reco_selection_part1->add<METCut>("MET > 50 GeV", 50.,100000.);
    reco_selection_part1->add<TwoDCut>("2D-Cut",0.4,25.);
    reco_selection_part1->add<NMuonBTagSelection>("b-jet in muon hemisphere", 1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM));

    // reco_selection_part2->add<JetIdSelection<TopJet>>("selected jet - p_{T} > 200 GeV", ctx, TopJetId(PtEtaCut(200.,100000.)) ,recotopjet_handlename);
  }else if(isVJetsSel){
    reco_selection_part1->add<NElectronSelection>("ele-veto",0,0);
    reco_selection_part1->add<NMuonSelection>("muon-veto",0,0);
    reco_selection_part1->add<NTopJetSelection>("N_{AK8} #geq 1, p_{T} > 500 GeV", 1,-1,TopJetId(PtEtaCut(500.,100000.)));
    reco_selection_part1->add<HTCut>("H_{T} > 1000 GeV",ctx, 1000.);

    // reco_selection_part2->add<JetIdSelection<TopJet>>("selected jet - p_{T} > 500 GeV", ctx, TopJetId(PtEtaCut(500.,100000.)) ,recotopjet_handlename);
    reco_selection_part2->add<RhoCut<TopJet>>("rhocut",ctx, -6.0, -2.1,recotopjet_handlename);
  }


  gen_selection_part1.reset(new AndSelection(ctx,"gen_selection_part1"));
  gen_selection_part2.reset(new AndSelection(ctx,"gen_selection_part2"));
  if(isTTbarSel){
    // gen_selection->add<TriggerSelection>("Trigger selection","HLT_Mu50_v*");
    gen_selection_part1->add<NGenTopJetSelection>("N_{gen,AK8} #geq 1, p_{T} > 200 GeV", 1,-1,GenTopJetId(PtEtaCut(200.,100000.)));
    gen_selection_part1->add<TTbarGenSemilepSelection>("gen semilep-selection",ctx, ttbargen_handlename,55.);
    // gen_selection->add<NMuonSelection>("N_{#mu} #geq 1", 1,1);
    // gen_selection->add<NElectronSelection>("N_{e} = 0", 0,0);
    gen_selection_part1->add<METCut>("MET > 50 GeV", 50.,100000.,true);
    // gen_selection->add<TwoDCut>("2D-Cut",0.4,25.);
    // gen_selection->add<NMuonBTagSelection>("b-jet in muon hemisphere", 1, 999, DeepJetBTag(DeepJetBTag::WP_MEDIUM));
    // gen_selection_part2->add<JetIdSelection<GenTopJet>>("selected jet - p_{T} > 200 GeV", ctx, GenTopJetId(PtEtaCut(200.,100000.)) ,gentopjet_handlename); // used to make sure we only get jets over pt-threshold (only useful if we do not select leading jet in pT, e.g. for the test where i picked the trailing jet in N2)
  }else if(isVJetsSel){
    // gen_selection->add<GenParticleIdSelection>("genele-veto",GenParticleId(GenParticlePDGIdId(11)),0,0);
    // gen_selection->add<GenParticleIdSelection>("genmuon-veto",GenParticleId(GenParticlePDGIdId(13)),0,0);
    // uhh2::Event::Handle<std::vector<Jet>> gentopjet_handle = ctx.get_handle<std::vector<GenTopJet>> (ctx.get("GenTopJetCollection"));
    gen_selection_part1->add<NGenTopJetSelection>("N_{gen,AK8} #geq 1, p_{T} > 500 GeV", 1,-1,GenTopJetId(PtEtaCut(500.,100000.)));
    // gen_selection_part1->add<HTCut>("H_{T} > 1000 GeV",ctx, 1000., infinity, genHT_handlename);
    
    // gen_selection_part2->add<JetIdSelection<GenTopJet>>("selected jet - p_{T} > 500 GeV", ctx, GenTopJetId(PtEtaCut(500.,100000.)) ,gentopjet_handlename);
    gen_selection_part2->add<RhoCut<GenTopJet>>("rhocut",ctx, -6.0, -2.1,gentopjet_handlename);
  }
  
  
  // HISTOGRAMS
  hists.emplace_back(new ElectronHists(ctx, "ElectronHists"));
  hists.emplace_back(new EventHists(ctx, "EventHists"));
  hists.emplace_back(new MuonHists(ctx, "MuonHists"));
  hists.emplace_back(new JetHists(ctx, "JetHists"));
  hists.emplace_back(new TopJetHists(ctx, "TopJetHists"));

  h_pfhists_200to500.reset(new PFHists(ctx, "PFHists_200to500"));
  h_pfhists_500to1000.reset(new PFHists(ctx, "PFHists_500to1000"));
  h_pfhists_1000to2000.reset(new PFHists(ctx, "PFHists_1000to2000"));
  h_pfhists_2000to3000.reset(new PFHists(ctx, "PFHists_2000to3000"));
  h_pfhists_3000to4000.reset(new PFHists(ctx, "PFHists_3000to4000"));
  h_pfhists_4000to5000.reset(new PFHists(ctx, "PFHists_4000to5000"));

  h_pfhists_inclusive.reset(new PFHists(ctx, "PFHists_inclusive"));
  h_pfhists_500to550.reset(new PFHists(ctx, "PFHists_500to550"));
  h_pfhists_550to600.reset(new PFHists(ctx, "PFHists_550to600"));
  h_pfhists_600to675.reset(new PFHists(ctx, "PFHists_600to675"));
  h_pfhists_675to800.reset(new PFHists(ctx, "PFHists_675to800"));
  h_pfhists_800to1200.reset(new PFHists(ctx, "PFHists_800to1200"));
  h_pfhists_1200toInf.reset(new PFHists(ctx, "PFHists_1200toInf"));

  h_gen_hists_commonmodules.reset(new JetMassGenHists(ctx,"GenHistsCommonModules",ttbargen_handlename,genHT_handlename));
  h_gen_hists_gensel.reset(new JetMassGenHists(ctx, "GenHistsGenSel",ttbargen_handlename,genHT_handlename));


  const std::vector<double> pt_edges = {200.00,300.00,400.00,500.00,650.00,800.00,1200.00};

  const std::vector<double> msd_edges = {50.00,60.00,70.00,80.00,90.00,100.00,110.00,120.00,130.00,140.00,150.00,160.00,170.00,180.00,190.00,200.00,210.00,220.00,230.00,240.00,250.00,260.00,270.00,280.00,290.00,300.00};

  const std::vector<double> pt_edges_fine = {200.00,250.00,300.00,350.00,400.00,450.00,500.00,550.00,600.00,650.00,700.00,750.00,800.00,850.00,900.00,950.00,1000.00,1050.00,1100.00,1150.00,1200.00};
  
  const std::vector<double> msd_edges_fine = {50.00,55.00,60.00,65.00,70.00,75.00,80.00,85.00,90.00,95.00,100.00,105.00,110.00,115.00,120.00,125.00,130.00,135.00,140.00,145.00,150.00,155.00,160.00,165.00,170.00,175.00,180.00,185.00,190.00,195.00,200.00,205.00,210.00,215.00,220.00,225.00,230.00,235.00,240.00,245.00,250.00,255.00,260.00,265.00,270.00,275.00,280.00,285.00,290.00,295.00,300.00};

  h_unfolding_hists_sel_part1.reset(new UnfoldingHists(ctx,"unfolding_hists_sel_part1",
                                             msd_edges,
                                             pt_edges,
                                             reco_selection_handlename,
                                             gen_selection_handlename,
                                             matching_selection_handlename,
                                             recotopjet_handlename,
                                             gentopjet_handlename));

  h_unfolding_hists_sel_part2.reset(new UnfoldingHists(ctx,"unfolding_hists_sel_part2",
                                             msd_edges,
                                             pt_edges,
                                             reco_selection_handlename,
                                             gen_selection_handlename,
                                             matching_selection_handlename,
                                             recotopjet_handlename,
                                             gentopjet_handlename));

  h_unfolding_hists.reset(new UnfoldingHists(ctx,"unfolding_hists",
                                             msd_edges,
                                             pt_edges,
                                             reco_selection_handlename,
                                             gen_selection_handlename,
                                             matching_selection_handlename,
                                             recotopjet_handlename,
                                             gentopjet_handlename));
  
  h_unfolding_hists_fine.reset(new UnfoldingHists(ctx,"unfolding_hists_fine",
                                                  msd_edges_fine,
                                                  pt_edges_fine,
                                                  reco_selection_handlename,
                                                  gen_selection_handlename,
                                                  matching_selection_handlename,
                                                  recotopjet_handlename,
                                                  gentopjet_handlename));

  

  writer.reset(new WriteOutput(ctx,matching_selection_handlename));

  std::string N2DDT_file_path = "JetMass/Histograms/ddtmaps/QCD_2017_PFMass_smooth_gaus4p00sigma.root";
  std::string N2DDT_hist_name = "N2_v_pT_v_rho_0p05_smooth_gaus4p00sigma_maps_cleaner_PFMass";
  reco_n2ddt_computer.reset(new N2DDTComputer<TopJet>(N2DDT_file_path,N2DDT_hist_name));
  gen_n2ddt_computer.reset(new N2DDTComputer<GenTopJet>(N2DDT_file_path,N2DDT_hist_name));
  // TFile* f_N2DDT = new TFile(locate_file(N2DDT_file_path).c_str(),"READ");  
  // hist_N2DDT = (TH2D*)f_N2DDT->Get(N2DDT_hist_name.c_str());
}

bool JetMassModule::process(Event & event) {
  if(EXTRAOUT){
    cout << "JetMassModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  }

  if(is_mc)ttgen_producer->process(event);

  matching_selection_producer->process(event);
  // MatchingSelection matching_selection = event.get(handle_matching_selection);
  
  double genHT(-1.0);
  if(is_mc)genHT = calculateGenJetPtSum(event);
  event.set(handle_gen_HT,genHT);    
  if(EXTRAOUT) std::cout << "JetMassModule: ExtraHandles done!" << std::endl;
  
  // Throw away events with NVTX in buggy area for buggy samples
  if(is_mc && is_buggyPU){
    float n_true = event.genInfo->pileup_TrueNumInteractions();
    if(n_true < 10. || n_true > 72.) return false;
    if(EXTRAOUT) std::cout << "JetMassModule: BuggyPU Treatment done!" << std::endl;
  }

  // COMMON MODULES
  bool pass_common=common->process(event);
  if(!pass_common) return false;

  h_gen_hists_commonmodules->fill(event);
  if(EXTRAOUT) std::cout << "JetMassModule: CommonModules done!" << std::endl;

  //remove MC Events with very large unphysical weights
  if(is_mc && isVJetsSel && is_QCD){
    if(!mcSpikeKiller->passes(event)) return false;
    if(EXTRAOUT) std::cout << "JetMassModule: SpikeKiller done!" << std::endl;
  }

  if(is_mc)topPtReweighting->process(event);

  // AK8 JEC
  topjetCorr->process(event);

  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);
  if(is_mc)sort_by_pt<GenTopJet>(*event.gentopjets);
  if(EXTRAOUT) std::cout << "JetMassModule: TopJetCorrections done!" << std::endl;
  
  // CLEANER
  // ak4cleaner15->process(event);
  ak4cleaner->process(event);
  ak8cleaner->process(event);
  ak8cleaner_dRlep->process(event);
  genak8cleaner->process(event);
  genak4cleaner->process(event);

  if(EXTRAOUT) std::cout << "JetMassModule: Cleaner done!" << std::endl;

  if(is_mc){
    //find gen(top)jet from gen V boson and save it to handle for use in nloweight application
    // const Particle gen_V_particle = event.genparticles->at(get_V_index(*event.genparticles));
    // std::cout << "N jets: " << event.gentopjets->size() << " " << event.genjets->size() << std::endl;
    // const GenTopJet * v_gentopjet = closestParticle(gen_V_particle,*event.gentopjets);
    // const GenJet * v_genjet = closestParticle(gen_V_particle,*event.genjets);
    // const GenTopJet * v_gentopjet = find_Vmatched_jet<GenTopJet>(*event.gentopjets, matching_selection);
    // const GenJet * v_genjet = find_Vmatched_jet<GenJet>(*event.genjets, matching_selection);
    // std::cout << "genV jets pt: " << (v_gentopjet?v_gentopjet->pt():-1.0) << " " << (v_genjet?v_genjet->pt():-1.0) <<" " << gen_V_particle.pt() << std::endl; 
    // std::cout << "genV jets mass: " << v_gentopjet->softdropmass() << " " << v_genjet->v4().M() <<" " << gen_V_particle.v4().M() << std::endl; 
    // event.set(handle_gentopjet_matched_V,v_gentopjet);
    // nlo_weights->process(event);
  }
  
  if(event.topjets->size()>0){
    //PFHists
    float AK8_pt = event.topjets->at(0).pt();
    h_pfhists_inclusive->fill(event);
    if(AK8_pt>200 && AK8_pt<500)h_pfhists_200to500->fill(event);
    if(AK8_pt>500 && AK8_pt<1000)h_pfhists_500to1000->fill(event);
    if(AK8_pt>1000 && AK8_pt<2000)h_pfhists_1000to2000->fill(event);
    if(AK8_pt>2000 && AK8_pt<3000)h_pfhists_2000to3000->fill(event);
    if(AK8_pt>3000 && AK8_pt<4000)h_pfhists_3000to4000->fill(event);
    if(AK8_pt>4000 && AK8_pt<5000)h_pfhists_4000to5000->fill(event);
    if(AK8_pt>500 && AK8_pt<550)h_pfhists_500to550->fill(event);
    if(AK8_pt>550 && AK8_pt<600)h_pfhists_550to600->fill(event);
    if(AK8_pt>600 && AK8_pt<675)h_pfhists_600to675->fill(event);
    if(AK8_pt>675 && AK8_pt<800)h_pfhists_675to800->fill(event);
    if(AK8_pt>800 && AK8_pt<1200)h_pfhists_800to1200->fill(event);
    if(AK8_pt>1200)h_pfhists_1200toInf->fill(event);
  }
  
  bool pass_reco_selection_part1 = reco_selection_part1->passes(event);
  bool pass_reco_selection_part2(true);
  if(pass_reco_selection_part1){
    recojet_selector->process(event);
    pass_reco_selection_part2 = reco_selection_part2->passes(event);
  }
  if(EXTRAOUT) std::cout << "JetMassModule: RecoSelection done!" << std::endl;

    
  bool pass_gen_selection_part1(true),pass_gen_selection_part2(true);
  if(is_mc)pass_gen_selection_part1 = gen_selection_part1->passes(event);
  if(pass_gen_selection_part1 && is_mc){
    genjet_selector->process(event);
    pass_gen_selection_part2 = gen_selection_part2->passes(event);
  }
  if(EXTRAOUT) std::cout << "JetMassModule: GenSelection done!" << std::endl;

  
  //set reco selection to fail if topjet-gentopjet matching within R/2=0.4 fails.  
  bool pass_dR_reco_gen(true);
  if( (pass_gen_selection_part1 && pass_reco_selection_part1 )){
    double dR(numeric_limits<double>::infinity());
    if(event.is_valid(handle_gentopjet) && event.is_valid(handle_recotopjet)){
      const GenTopJet* genjet = event.get(handle_gentopjet);
      const TopJet* recojet = event.get(handle_recotopjet);
      dR = deltaR(*genjet,*recojet); 
    }
    if(dR > 0.4) pass_dR_reco_gen = false;
  }

  if(EXTRAOUT) std::cout << "JetMassModule: Reco-Gen Matching done!" << std::endl;

  bool pass_substructure_cut(true);  
  if(pass_reco_selection_part1){
    const TopJet * recotopjet = NULL;
    if(event.is_valid(handle_recotopjet))recotopjet=event.get(handle_recotopjet);
    if(isVJetsSel){
      double n2ddt = reco_n2ddt_computer->computeDDTValue(recotopjet);
      pass_substructure_cut = n2ddt < 0;
    }
    if(isTTbarSel){
      double tau32 = safe_tau32(recotopjet);
      pass_substructure_cut = tau32<0.5;
    }
  }
  if(EXTRAOUT) std::cout << "JetMassModule: SubstructureSelection done!" << std::endl;

  //fill first round of unfolding hists with part1 of selections
  //IMPORTANT: setting selection bits before filling UnfoldingHists, since those rely on correct setting of handles!
  event.set(handle_reco_selection,pass_reco_selection_part1);
  event.set(handle_gen_selection,pass_gen_selection_part1);
  h_unfolding_hists_sel_part1->fill(event);

  
  //fill second round of unfolding hists with part2 of selections
  bool pass_reco_selection = pass_reco_selection_part1 && pass_reco_selection_part2 && pass_substructure_cut;
  bool pass_gen_selection = pass_gen_selection_part1 && pass_gen_selection_part2;
  event.set(handle_reco_selection, pass_reco_selection);
  event.set(handle_gen_selection, pass_gen_selection);
  h_unfolding_hists_sel_part2->fill(event);  
  if(pass_gen_selection) h_gen_hists_gensel->fill(event);

  //last round of unfolding hists with dR matching
  event.set(handle_reco_selection, pass_reco_selection && pass_dR_reco_gen);
  h_unfolding_hists->fill(event);
  h_unfolding_hists_fine->fill(event);

  if(EXTRAOUT) std::cout << "JetMassModule: UnfoldingHistFilling done!" << std::endl;

  //Write everything used for JetMassCalibration to Tree if first part of reco-selection passes
  if(pass_reco_selection_part1){  
    // FILL HISTS
    for(auto & h: hists) h->fill(event);
    // STORE EVENT
    writer->process(event);
    
    if(EXTRAOUT) std::cout << "JetMassModule: OutputTreeWriter done!" << std::endl;
    return true;
  }else return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the JetMassModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(JetMassModule)
