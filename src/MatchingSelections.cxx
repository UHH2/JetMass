#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/core/include/Event.h"

using namespace std;
using namespace uhh2;
// namespace uhh2{
  MatchingSelection::MatchingSelection(uhh2::Context & ctx){
    TString version = ctx.get("dataset_version", "");
    version.ToLower();
    is_VJets = (version.Contains("wjets") || version.Contains("zjets") );
    // is_TTbar = ((version.Contains("ttbar") || version.Contains("ttjets") || version.Contains("ttto")) && version.Contains("semilep")) || version.Contains("tt_mtt");
    is_TTbar = ((version.Contains("ttbar") || version.Contains("ttjets") || version.Contains("ttto")));
    is_HiggsToWW = (version.Contains("higgs"));
    is_valid = true;
  }

  void MatchingSelection::init(const uhh2::Event &event){
    //reset is_valid to true
    is_valid = true;
    if(event.isRealData){
      is_valid = false;
      return;
    }

    if(!is_VJets && !is_TTbar && !is_HiggsToWW){
      // cout << "This is not a TTbar or V+Jets sample. Matching is not implemented. Skipping MatchingSelection" << endl;
      is_valid = false;
      return;
    }
    extra_partons.clear();
    assert(event.genparticles);
    if(is_VJets){
      unsigned int n_V = 0;
      for(unsigned int i=0; i<event.genparticles->size(); ++i) {
        const GenParticle & genp = event.genparticles->at(i);
        if(abs(genp.pdgId()) == 23 || abs(genp.pdgId()) == 24){
          ++n_V;
          genV = genp;
          genQ1 = *genp.daughter(event.genparticles,1);
          genQ2 = *genp.daughter(event.genparticles,2);
        }
        if(
          (genp.mother1() != genV.index()) & (genp.status() ==23) &
          (abs(genp.pdgId()) == 1 || abs(genp.pdgId()) == 2 || abs(genp.pdgId()) == 3 ||
           abs(genp.pdgId()) == 4 || abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 21)
          ){
          extra_partons.push_back(genp);
        }
      }
      if(n_V > 1){
        cout << "MatchinSelection (V+jets): There are more than 1 W/Z Boson among generator-particles! (event:"<<event.event << ")" << endl;
        is_valid = false;
      }
    }else if(is_TTbar){
      TTbarGen ttgen(*event.genparticles);
      if(ttgen.IsSemiLeptonicDecay()){
        genTop = ttgen.TopHad();
        genB = ttgen.BHad();
        genV = ttgen.WHad();
        genQ1 = ttgen.Q1();
        genQ2 = ttgen.Q2();
      }else{
        is_valid = false;
      }
    }else if(is_HiggsToWW){
      for(unsigned int i=0; i<event.genparticles->size(); ++i) {
        const GenParticle & genp = event.genparticles->at(i);
        if(abs(genp.pdgId())==25){
          genH = genp;
          GenParticle W1 = *genp.daughter(event.genparticles,1);
          genQ1 = *W1.daughter(event.genparticles,1);
          genQ2 = *W1.daughter(event.genparticles,2);
          GenParticle W2 = *genp.daughter(event.genparticles,2);
          genQ3 = *W2.daughter(event.genparticles,1);
          genQ4 = *W2.daughter(event.genparticles,2);          
        }        
      }
    }else{
      is_valid = false; 
    }
  }

  bool MatchingSelection::passes_matching(const FlavorParticle &probe_jet, matchingOpt opt, float radius) {
    if(!is_VJets && !is_TTbar && !is_HiggsToWW) return false;

    if(!is_valid){
      if(opt == oIsNotMerged) return true;
      else return false;
    }

    bool b_in_Jet = is_TTbar ? deltaR(probe_jet, genB) < radius : false;
    bool q1_in_Jet = dR_Q1(probe_jet) < radius;
    bool q2_in_Jet = dR_Q2(probe_jet) < radius;

    bool q3_in_Jet = is_HiggsToWW ? deltaR(probe_jet, genQ3)< radius : false;
    bool q4_in_Jet = is_HiggsToWW ? deltaR(probe_jet, genQ4)< radius : false;

    if(is_HiggsToWW && (opt==oIsMergedHiggs)) return (q1_in_Jet && q2_in_Jet && q3_in_Jet && q4_in_Jet);
    
    if(b_in_Jet && q1_in_Jet && q2_in_Jet ){
      if( (opt == oIsMergedTop)) return true;
    }else if( ((b_in_Jet && !q1_in_Jet && q2_in_Jet) || (b_in_Jet && q1_in_Jet && !q2_in_Jet)) ){
      if( ((opt == oIsMergedQB) || (opt == oIsSemiMergedTop)) ) return true;
    }else if( (!b_in_Jet && q1_in_Jet && q2_in_Jet) ){
      if( ((opt == oIsMergedV) || (opt == oIsMergedW) || (opt == oIsMergedZ) || (opt == oIsSemiMergedTop)) ) return true;
    }else if(opt == oIsNotMerged) return true;
    return false;
  }

  int MatchingSelection::FlavourQ1() const{
    if(!is_VJets && !is_TTbar) return 10000;
    return genQ1.pdgId();
  }

  int MatchingSelection::FlavourQ2() const{
    if(!is_VJets && !is_TTbar) return 10000;
    return genQ2.pdgId();
  }

  double MatchingSelection::dR_Q1(const FlavorParticle &probe_jet) {return deltaR(probe_jet, genQ1);}
  double MatchingSelection::dR_Q2(const FlavorParticle &probe_jet) {return deltaR(probe_jet, genQ2);}

  int MatchingSelection::n_merged_partons(const FlavorParticle &probe_jet,float radius) {
    int n(0);
    
    // std::cout << "pt W: " << genV.pt()  << std::endl;
    // std::cout << "pt W daughters: " << genQ1.pt() << " " << genQ2.pt() << std::endl;
    // std::cout << "n extra partons " << extra_partons.size() << std::endl;
    for(unsigned int i=0; i<extra_partons.size(); i++){
      GenParticle genp = extra_partons[i];
      bool in_jet = deltaR(probe_jet, genp) <= radius;
      
      if(in_jet){
        // std::cout << genp.pt() << " " << probe_jet.pt() << std::endl;
        // for(unsigned int isj=0; isj < probe_jet.subjets().size(); isj++) {
        //   FlavorParticle subjet = probe_jet.subjets().at(isj);
        //   std::cout << "subjet " <<isj << " " << subjet.pt() << std::endl; 
        // }
        n++;
      }
    }
    return n;
    
  }

  MatchingSelectionProducer::MatchingSelectionProducer(uhh2::Context & ctx, const std::string & name){
    matching_selection.reset(new MatchingSelection(ctx));
    h_matching_selection = ctx.get_handle<MatchingSelection>(name);
  }
  
  bool MatchingSelectionProducer::process(Event & event){
    matching_selection->init(event);
    event.set(h_matching_selection, *matching_selection);
    return true;
  }

// }
