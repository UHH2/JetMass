#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/Utils.h" //mainly for runPeriods vectors
#include "UHH2/core/include/TopJet.h"

#include "UHH2/JetMass/include/H3DDTHist.h"
#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/JetMass/include/JetMassHists_central.h"
#include "UHH2/JetMass/include/JetMassHists_central_noConst.h"
#include "UHH2/JetMass/include/SubstructureSelections.h"
#include "UHH2/JetMass/include/MatchingSelections.h"

#include "UHH2/JetMass/include/PFHists.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"


#include <unistd.h>


#define EXTRAOUT false


using namespace std;
using namespace uhh2;


	Particle GetGeneratorW(Event & event);

	class WMassPreSelModule: public AnalysisModule {
	public:

    explicit WMassPreSelModule(Context & ctx);
    virtual bool process(Event & event) override;

	private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCleaner> ak4cleaner;
    std::unique_ptr<TopJetCleaner> ak8cleaner;

		std::unique_ptr<Selection> n2_physical_sel,N_AK8_500_sel;
		std::unique_ptr<Selection> rho_sel;
		std::unique_ptr<MatchingSelection> matching_sel;

		std::unique_ptr<Hists> h_nocuts,h_jec,h_cleaner,h_nak8_500sel;
		std::unique_ptr<Hists> h_nocuts_noConst,h_jec_noConst,h_cleaner_noConst,h_nak8_500sel_noConst;
		std::unique_ptr<Hists> h_rhoCut,h_rhoCut_noConst;

		std::unique_ptr<Hists> h_N2_3D;

		std::unique_ptr<YearSwitcher> AK4_JEC_MC,AK8_JEC_MC;
		std::unique_ptr<YearSwitcher> AK4_JEC_data,AK8_JEC_data;
		std::shared_ptr<RunSwitcher> AK4_JEC_Switcher_16,AK4_JEC_Switcher_17,AK4_JEC_Switcher_18;
		std::shared_ptr<RunSwitcher> AK8_JEC_Switcher_16,AK8_JEC_Switcher_17,AK8_JEC_Switcher_18;

		std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
    std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;

		std::unique_ptr<AnalysisModule> pfparticles_jec_corrector,pf_applyPUPPI;

		std::unique_ptr<JetResolutionSmearer> AK4_jet_smearer;
		std::unique_ptr<GenericJetResolutionSmearer> AK8_jet_smearer;

		bool is_mc,matchW,is_WSample;

		Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;
	};


	WMassPreSelModule::WMassPreSelModule(Context & ctx){
		is_mc = ctx.get("dataset_type") == "MC";

		std::string version=ctx.get("dataset_version");
		is_WSample=version.find("WJets") != std::string::npos;
		matchW=version.find("WMatched") != std::string::npos;

    common.reset(new CommonModules());

		common->disable_jec();
		common->disable_jersmear();

		common->init(ctx);

		const std::string JEC_tag_2016="Summer16_07Aug2017";
		const std::string JEC_version_2016="11";
		const std::string JEC_tag_2017="Fall17_17Nov2017";
		const std::string JEC_version_2017="32";
		const std::string JEC_tag_2018="Autumn18";
		const std::string JEC_version_2018="8";
		const std::string AK4_jetcoll="AK4PFPuppi";
		const std::string AK8_jetcoll="AK8PFPuppi";

		if(is_mc){
      AK4_JEC_MC.reset(new YearSwitcher(ctx));
      AK4_JEC_MC->setup2016(std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2016, JEC_version_2016, AK4_jetcoll)));
      AK4_JEC_MC->setup2017(std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2017, JEC_version_2017, AK4_jetcoll)));
      AK4_JEC_MC->setup2018(std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2018, JEC_version_2018, AK4_jetcoll)));

      AK8_JEC_MC.reset(new YearSwitcher(ctx));
      AK8_JEC_MC->setup2016(std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2016, JEC_version_2016, AK8_jetcoll)));
      AK8_JEC_MC->setup2017(std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2017, JEC_version_2017, AK8_jetcoll)));
      AK8_JEC_MC->setup2018(std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesMC(JEC_tag_2018, JEC_version_2018, AK8_jetcoll)));

			// const JERSmearing::SFType1 AK8_JER_sf=JERSmearing::SF_13TeV_Autumn18_RunABCD_V4;
			// const TString resFilename="2018/Autumn18_V4_MC_PtResolution_AK8PFPuppi.txt";
			// AK4_jet_smearer.reset(new JetResolutionSmearer(ctx));
			// if(EXTRAOUT)std::cout << "AK4jetER_smearer set up!" << std::endl;
			// AK8_jet_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets","gentopjets", & AK8_JER_sf, resFilename));
			// if(EXTRAOUT)std::cout << "AK8jetER_smearer set up!" << std::endl;
		}else{
			AK4_JEC_Switcher_16.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriods2016) {
        AK4_JEC_Switcher_16->setupRun(runItr, std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2016, JEC_version_2016, AK4_jetcoll, runItr)));
      }

			AK4_JEC_Switcher_17.reset(new RunSwitcher(ctx, "2017"));
      for (const auto & runItr : runPeriods2017) {
        AK4_JEC_Switcher_17->setupRun(runItr, std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2017, JEC_version_2017, AK4_jetcoll, runItr)));
      }

      AK4_JEC_Switcher_18.reset(new RunSwitcher(ctx, "2018"));
      for (const auto & runItr : runPeriods2018) {
        AK4_JEC_Switcher_18->setupRun(runItr, std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2018, JEC_version_2018, AK4_jetcoll, runItr)));
      }

      AK4_JEC_data.reset(new YearSwitcher(ctx));
      AK4_JEC_data->setup2016(AK4_JEC_Switcher_16);
      AK4_JEC_data->setup2017(AK4_JEC_Switcher_17);
      AK4_JEC_data->setup2018(AK4_JEC_Switcher_18);

			AK8_JEC_Switcher_16.reset(new RunSwitcher(ctx, "2016"));
      for (const auto & runItr : runPeriods2016) {
        AK8_JEC_Switcher_16->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2016, JEC_version_2016, AK8_jetcoll, runItr)));
      }

			AK8_JEC_Switcher_17.reset(new RunSwitcher(ctx, "2017"));
      for (const auto & runItr : runPeriods2017) {
        AK8_JEC_Switcher_17->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2017, JEC_version_2017, AK8_jetcoll, runItr)));
      }

      AK8_JEC_Switcher_18.reset(new RunSwitcher(ctx, "2018"));
      for (const auto & runItr : runPeriods2018) {
        AK8_JEC_Switcher_18->setupRun(runItr, std::make_shared<TopJetCorrector>(ctx, JERFiles::JECFilesDATA(JEC_tag_2018, JEC_version_2018, AK8_jetcoll, runItr)));
      }

      AK8_JEC_data.reset(new YearSwitcher(ctx));
      AK8_JEC_data->setup2016(AK8_JEC_Switcher_16);
      AK8_JEC_data->setup2017(AK8_JEC_Switcher_17);
      AK8_JEC_data->setup2018(AK8_JEC_Switcher_18);
		}

		if(is_mc){
			pfparticles_jec_corrector.reset(new CorrectParticles());
			pf_applyPUPPI.reset(new ApplyPuppiToPF());
		}

		AK4_Clean_pT=30.0;
		AK4_Clean_eta=5.0;
		AK8_Clean_pT=200.0;
		AK8_Clean_eta=2.5;
		ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
    ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));

		// SELECTIONS
		// n2_physical_sel.reset(new N2Selection(1.0,0.0));


		N_AK8_500_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(500.,100000.))));
		rho_sel.reset(new RhoSelection(-6.0,-2.1));
		matching_sel.reset(new MatchingSelection(ctx,MatchingSelection::oIsLeadingGenW));
		// HISTOGRAMS
	  double variation = 0.1;

		if(is_mc){
			h_nocuts.reset(new JetMassHists(ctx, "JetMass_NoCuts", variation, "SD"));
			h_jec.reset(new JetMassHists(ctx, "JetMass_JetCorrections", variation, "SD"));
			h_cleaner.reset(new JetMassHists(ctx, "JetMass_AfterCleaner", variation, "SD"));
			h_nak8_500sel.reset(new JetMassHists(ctx, "JetMass_Nak8_500sel", variation, "SD"));
			h_N2_3D.reset(new H3DDTHist(ctx,"N2_3DMap"));
			h_rhoCut.reset(new JetMassHists(ctx, "JetMass_RhoCut_sel", variation, "SD"));

				}
		h_nocuts_noConst.reset(new JetMassHists_central_noConst(ctx, "JetMass_NoCuts_noConst"));
		h_jec_noConst.reset(new JetMassHists_central_noConst(ctx, "JetMass_JetCorrections_noConst"));
		h_cleaner_noConst.reset(new JetMassHists_central_noConst(ctx, "JetMass_AfterCleaner_noConst"));
		h_nak8_500sel_noConst.reset(new JetMassHists_central_noConst(ctx, "JetMass_Nak8_500sel_noConst"));
		h_rhoCut_noConst.reset(new JetMassHists(ctx, "JetMass_RhoCut_sel_noConst", variation, "SD"));

	}


	bool WMassPreSelModule::process(Event & event) {
    if(EXTRAOUT){
			cout << "WMassPreSelModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
		}
		if(is_mc && is_WSample){
			auto GeneratorW = GetGeneratorW(event);
			if(event.topjets->size() < 1) return false;
			bool Wmatched=deltaR(event.topjets->at(0),GeneratorW) < 0.8;
			if(matchW){
				if(!Wmatched) return false;
			}else{
				if(Wmatched) return false;
			}
		}
    bool pass_common=common->process(event);
		if(!pass_common) return false;
		if(is_mc){
    	h_nocuts->fill(event);
		}
		h_nocuts_noConst->fill(event);

		//Beginning of JEC
		// Applying JEC to AK4 (event.jets) and AK8 (event.topjets)
		bool printJEC=false;
		if(printJEC) std::cout << "-------------------------------------------------------" <<std::endl;
		if(event.jets->size()>0 && printJEC){
			std::cout << " jet pT * JEC_factor_raw -> raw_pT" << std::endl;
			std::cout << event.jets->at(0).pt() << " * " <<event.jets->at(0).JEC_factor_raw()<< " -> " << event.jets->at(0).pt()*event.jets->at(0).JEC_factor_raw()<<std::endl;
		}
		if(event.topjets->size()>0 && printJEC){
			std::cout << " topjet pT * JEC_factor_raw -> raw_pT" << std::endl;
			std::cout << event.topjets->at(0).pt() << " * " <<event.topjets->at(0).JEC_factor_raw()<< " -> " << event.topjets->at(0).pt()*event.topjets->at(0).JEC_factor_raw()<<std::endl;
		}
		if(is_mc){
			AK4_JEC_MC->process(event);
			AK8_JEC_MC->process(event);
			// AK4_jet_smearer->process(event);
			// AK8_jet_smearer->process(event);
		}else{
			AK4_JEC_data->process(event);
			AK8_JEC_data->process(event);
		}

		if(is_mc)pfparticles_jec_corrector->process(event);
		if(is_mc)pf_applyPUPPI->process(event);

		sort_by_pt<Jet>(*event.jets);
		sort_by_pt<TopJet>(*event.topjets);

		//End of JEC
		if(is_mc){
			h_jec->fill(event);
		}
		h_jec_noConst->fill(event);
    ak4cleaner->process(event);
    ak8cleaner->process(event);
		// if(!n2_physical_sel->passes(event)) return false;

		if(is_mc){
			h_cleaner->fill(event);
			h_N2_3D->fill(event);
		}
		h_cleaner_noConst->fill(event);

		bool N_AK8_500_selection=N_AK8_500_sel->passes(event);
		if(!N_AK8_500_selection) return false;
		if(is_mc)h_nak8_500sel->fill(event);
 		h_nak8_500sel_noConst->fill(event);

		bool rhoCut = rho_sel->passes(event);
		if(!rhoCut) return false;
		if(is_mc)h_rhoCut->fill(event);
		h_rhoCut_noConst->fill(event);
		
		TopJet probejet = event.topjets->at(0);
		if(matching_sel->passes_matching(event,probejet)) return true;

		return true;
	}

	// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
	// make sure the WMassPreSelModule is found by class name. This is ensured by this macro:
	UHH2_REGISTER_ANALYSIS_MODULE(WMassPreSelModule)

	Particle GetGeneratorW(Event & event ){
		int GenW_index=-1;
		for(unsigned int i = 0 ; i < event.genparticles->size() ; ++i){
			if(abs(event.genparticles->at(i).pdgId()) == 24){
				if(GenW_index == -1){
					GenW_index=i;
				}else{
					std::cout << "WARNING: More than one W on generator level. choosing W with higher pT!" << std::endl;
					GenW_index= i ? event.genparticles->at(i).pt() > event.genparticles->at(GenW_index).pt() : GenW_index;
				}
			}
		}
		return event.genparticles->at(GenW_index);
	}
