#include "UHH2/JetMass/include/MatchingSelections.h"
#include "UHH2/core/include/Event.h"

using namespace std;
// using namespace uhh2;
namespace uhh2{
  MatchingSelection::MatchingSelection(uhh2::Context & ctx, matchingOpt opt_):opt(opt_){
    string version = ctx.get("dataset_version", "");
    std::cout << "MatchingSelection: Using Selection with option: " << opt << std::endl;
  }

  bool MatchingSelection::passes_matching(const uhh2::Event &event, const TopJet &probe_jet) {
    if(opt == oIsAnyGenW || opt == oIsLeadingGenW){
      //finding all W in genparticles  
      assert(event.genparticles);
      std::vector<GenParticle> GenWs={};
      for(unsigned int i = 0 ; i < event.genparticles->size() ; i++){
          GenParticle thisGen = event.genparticles->at(i); 
          if(abs(thisGen.pdgId()) == 24 ) GenWs.push_back(thisGen);
          //if matching with ANY W is wanted lets do it here
          if(opt == oIsAnyGenW && deltaR(thisGen,probe_jet) < 0.8)return true;
      }
      //sort all found W by its pT
      sort_by_pt<GenParticle>(GenWs);
      //if there is any W match the leading
      if(GenWs.size()>0){
        if(opt == oIsLeadingGenW && deltaR(GenWs[0],probe_jet) < 0.8) return true;
      }
    }
    return false;
  }
}
