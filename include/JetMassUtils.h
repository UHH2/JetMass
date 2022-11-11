#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/GenTopJet.h"
#include "UHH2/JetMass/include/MatchingSelections.h"

#include <TFile.h>
#include <TH2D.h>

// const GenParticle hardestParticle(uhh2::Event & event, const GenTopJet & jet);

double calculateGenTopPtSum(uhh2::Event & event);

double calculateGenJetPtSum(uhh2::Event & event);

int get_V_index(std::vector<GenParticle> particles);

template<typename GenericTopJet>
class JetSelector: public uhh2::AnalysisModule{
public:
  explicit JetSelector(uhh2::Context & ctx, const std::string & jet_handlename, const std::string & matching_selection_handlename=""):
    jet_handle(ctx.get_handle<const GenericTopJet*>(jet_handlename)),
    matching_selection_handle(ctx.get_handle<MatchingSelection>(matching_selection_handlename)){
    v_matching = (ctx.get("dataset_version").find("WJets") != std::string::npos ) && (ctx.get("selection").find("vjets") != std::string::npos);
  };
  virtual bool process(uhh2::Event & event) override;
  
private:
  uhh2::Event::Handle<const GenericTopJet*> jet_handle;
  uhh2::Event::Handle<MatchingSelection> matching_selection_handle;
  bool v_matching;
};


class NLOWeights: public uhh2::AnalysisModule{
public:
  NLOWeights(uhh2::Context & ctx, const std::string & boson_pt_handlename="V_pt");
  virtual bool process(uhh2::Event & event) override;
  
private:
  uhh2::Event::Handle<double> boson_pt_handle;  
  TH1F *h_kfactor, *h_ewcorr;
  bool apply_nloweights;
};

template<typename GenericTopJet>
float safe_tau21(const GenericTopJet * jet){
  if(jet == nullptr) throw std::runtime_error("safe_tau21: passed jet is nullptr!");
  float tau1 = jet->tau1();
  float tau2 = jet->tau2();
  return  tau1>0 ? tau2/tau1 : 0.0;
};


template<typename GenericTopJet>
float safe_tau32(const GenericTopJet * jet){
  if(jet == nullptr) throw std::runtime_error("safe_tau32: passed jet is nullptr!");
  float tau2 = jet->tau2();
  float tau3 = jet->tau3();
  return  tau2>0 ? tau3/tau2 : 0.0;
};

template<typename GenericTopJet>
class N2DDTComputer{
public:
  N2DDTComputer(const std::string & map_file_path, const std::string map_name);
  double computeDDTValue(const GenericTopJet * jet);

private:
  std::unique_ptr<TH2D> h_n2ddt_map;
};

template class N2DDTComputer<TopJet>;
template class N2DDTComputer<GenTopJet>;


template<typename GenericJet>
const GenericJet * find_Vmatched_jet(std::vector<GenericJet> & jets, MatchingSelection & matching_selection){
  std::vector<GenericJet> matched_jets;
  for(auto thisjet:jets){
    if(matching_selection.passes_matching(thisjet,MatchingSelection::oIsMergedV)) matched_jets.push_back(thisjet);
  }
  std::cout << "found " << matched_jets.size() << " matching jets. picking closest to gen particle" << std::endl;
  const GenericJet *jet(NULL);
  if(matched_jets.size()>0){
    jet = closestParticle(*matching_selection.get_genV(),matched_jets);
  }else{
    jet = closestParticle(*matching_selection.get_genV(),jets);
  }
  return jet;
};

int PFMultiplicity(std::vector<PFParticle>,int);

float PFEnergy(std::vector<PFParticle>,int);

float PFTransverseMomentum(std::vector<PFParticle>,int);


