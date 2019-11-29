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

#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MuonHists.h"


#include "UHH2/JetMass/include/TopJetCorrections.h"
#include "UHH2/JetMass/include/CorrectParticles.h"
#include "UHH2/JetMass/include/ApplyPuppiToPF.h"


#include <unistd.h>


#define EXTRAOUT false


using namespace std;
using namespace uhh2;

namespace uhh2examples {

	Particle GetGeneratorW(Event & event);

	class PreSelModule: public AnalysisModule {
	public:

		explicit PreSelModule(Context & ctx);
		virtual bool process(Event & event) override;

	private:

		std::unique_ptr<CommonModules> common;
		std::unique_ptr<TopJetCorrections> topjetCorr;

		std::unique_ptr<JetCleaner> ak4cleaner;
		std::unique_ptr<TopJetCleaner> ak8cleaner;

		std::unique_ptr<Selection> N_AK8_300_sel, N_AK8_500_sel;

		std::vector<std::unique_ptr<uhh2::Hists>> hists;

		std::unique_ptr<AnalysisModule> pfparticles_jec_corrector,pf_applyPUPPI;

		bool is_mc,matchW,is_WSample;
		bool isTopSel = false;
		bool isWSel = false;

		Double_t AK4_Clean_pT,AK4_Clean_eta,AK8_Clean_pT,AK8_Clean_eta;
	};


	PreSelModule::PreSelModule(Context & ctx){

		// Set some boolians
		is_mc = ctx.get("dataset_type") == "MC";

		std::string version=ctx.get("dataset_version");
		is_WSample=version.find("WJets") != std::string::npos;
		matchW=version.find("WMatched") != std::string::npos;

		const std::string& channel = ctx.get("channel", "");
		if     (channel == "top") isTopSel = true;
		else if(channel == "W") isWSel = true;

		// common modules
		common.reset(new CommonModules());
		common->init(ctx);

		// AK8 JEC/JER
		topjetCorr.reset(new TopJetCorrections());
		topjetCorr->init(ctx);

		// PF correctors
		if(is_mc || !isWSel){
			pfparticles_jec_corrector.reset(new CorrectParticles());
			pf_applyPUPPI.reset(new ApplyPuppiToPF());
		}

		// Jet cleaner
		AK4_Clean_pT = 30.0;
		AK4_Clean_eta = 5.0;
		AK8_Clean_pT = 200.0;
		AK8_Clean_eta = 2.5;
		ak4cleaner.reset(new JetCleaner(ctx, AK4_Clean_pT, AK4_Clean_eta));
		ak8cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(AK8_Clean_pT,AK8_Clean_eta))));

		// SELECTIONS
		N_AK8_300_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(300.,100000.))));
		N_AK8_500_sel.reset(new NTopJetSelection(1,-1,TopJetId(PtEtaCut(500.,100000.))));

		// HISTOGRAMS
		hists.emplace_back(new ElectronHists(ctx, "ElectronHists"));
		hists.emplace_back(new EventHists(ctx, "EventHists"));
		hists.emplace_back(new MuonHists(ctx, "MuonHists"));
		hists.emplace_back(new JetHists(ctx, "JetHists"));
		hists.emplace_back(new TopJetHists(ctx, "TopJetHists"));


	}


	bool PreSelModule::process(Event & event) {
		if(EXTRAOUT){
			cout << "PreSelModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
		}
		// MATCHING

		// COMMON MODULES
		bool pass_common=common->process(event);
		if(!pass_common) return false;

		// AK8 JEC
		topjetCorr->process(event);

		if(is_mc || !isWSel) pfparticles_jec_corrector->process(event);
		if(is_mc || !isWSel) pf_applyPUPPI->process(event);

		sort_by_pt<Jet>(*event.jets);
		sort_by_pt<TopJet>(*event.topjets);

		// CLEANER
		ak4cleaner->process(event);
		ak8cleaner->process(event);

		// SELECTIONS
		bool passTOP = true;
		bool passW = true;

		if(isTopSel) passW = false;
		else if(isWSel) passTOP = false;

		if(!N_AK8_300_sel->passes(event)) passTOP = false;
		if(!N_AK8_500_sel->passes(event)) passW = false;

		double rho=2*TMath::Log(event.topjets->at(0).softdropmass()/event.topjets->at(0).pt());
		bool rhoCut= ( rho < -2.1 ) && ( rho > -6.0);
		if(!rhoCut) passW = false;

		if(passTOP && isTopSel)  for(auto & h: hists) h->fill(event);
		else if(passW && isWSel) for(auto & h: hists) h->fill(event);

		if(!passTOP && !passW) return false;
		else return true;
	}

	// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
	// make sure the PreSelModule is found by class name. This is ensured by this macro:
	UHH2_REGISTER_ANALYSIS_MODULE(PreSelModule)


}
