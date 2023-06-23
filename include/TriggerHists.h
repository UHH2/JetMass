#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"


using namespace std;

class TriggerHists: public uhh2::Hists {
public:
  TriggerHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~TriggerHists();

private:
  bool isMC;
  uhh2::Event::Handle<double> h_ht;
  JetId jet_tight_id;
  bool is_SingleMuon; 
  std::unique_ptr<NMuonSelection> nmuon_sel;
  MuonId muid_fortriggereff = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(30., 2.4));
  std::vector<TString> probe_triggers;
};
