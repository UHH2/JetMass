#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include <vector>
#include "boost/algorithm/string.hpp"
// namespace uhh2{
using namespace uhh2;

class MatchingSelection{
public:
  // enum matchingOpt{oIsAnyGenW,oIsLeadingGenW,oIsMergedGenW,oIsAnyGenZ,oIsLeadingGenZ,oIsMergedGenZ,oIsMergedW,};
  enum matchingOpt{oIsMergedV,oIsMergedW,oIsMergedZ,oIsMergedQB,oIsMergedTop,oIsSemiMergedTop,oIsNotMerged,oIsMergedHiggs};
  
  explicit MatchingSelection(uhh2::Context & ctx);
  
  //this will not be used since we want to pass a (top)jet to be probed.
  //In order to ommit the warning that 'event' is not used we print this useless warning.
  // virtual bool passes(const uhh2::Event & event) override{
  //   std::cout << "MatchingSelection: Using empty passes function on event "<< event.event << ". Are you sure you want to do that?" << std::endl;
  //   return true;
  // };

  void init(const uhh2::Event & event);
  
  bool passes_matching(const FlavorParticle &probe_jet, matchingOpt opt, float radius = 0.8);

  int FlavourQ1() const;
  int FlavourQ2() const;

  double dR_Q1(const FlavorParticle &probe_jet);
  double dR_Q2(const FlavorParticle &probe_jet);

  int n_merged_partons(const FlavorParticle &probe_jet, float radius = 0.8);

  const GenParticle * get_genV(){return &genV;};
  
private:
  bool initialzed, is_VJets, is_TTbar, is_HiggsToWW, is_valid;
  GenParticle genV,genTop,genB,genQ1,genQ2;

  GenParticle genH,genQ3,genQ4;
  
  std::vector<GenParticle> extra_partons;
};



class MatchingSelectionProducer: public uhh2::AnalysisModule {
 public:
  explicit MatchingSelectionProducer(uhh2::Context & ctx, const std::string &name = "matching_selection");
  virtual bool process(uhh2::Event & event) override;

 private:
  uhh2::Event::Handle<MatchingSelection> h_matching_selection;
  std::unique_ptr<MatchingSelection> matching_selection;
};

// }
