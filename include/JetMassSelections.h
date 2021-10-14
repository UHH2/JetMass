#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include <vector>


template<typename GenericJet>
class JetIdSelection: public uhh2::Selection {
public:
  explicit JetIdSelection(uhh2::Context & ctx, const std::function<bool (const GenericJet &, const uhh2::Event &)> jet_id_, const std::string jet_handlename=""):jet_handle(ctx.get_handle<const GenericJet*>(jet_handlename)),jet_id(jet_id_){};
  virtual bool passes(const uhh2::Event & event) override;  
private:
  uhh2::Event::Handle<const GenericJet*> jet_handle;
  const std::function<bool (const GenericJet &, const uhh2::Event &)> jet_id;
};

template class JetIdSelection<Jet>;
template class JetIdSelection<TopJet>;
template class JetIdSelection<GenJet>;
template class JetIdSelection<GenTopJet>;



template<typename GenericTopJet>
class RhoCut: public uhh2::Selection {
public:
  explicit RhoCut(uhh2::Context & ctx, double rho_min_ = - infinity, double rho_max_ = infinity, const std::string & jet_handlename=""): rho_min(rho_min_),rho_max(rho_max_),jet_handle(ctx.get_handle<const GenericTopJet*>(jet_handlename)){};
  virtual bool passes(const uhh2::Event & event) override;
  
private:
  double rho_min,rho_max;
  uhh2::Event::Handle<const GenericTopJet*> jet_handle;
};
template class RhoCut<TopJet>;
template class RhoCut<GenTopJet>;

class MassCut: public uhh2::Selection {
public:
    MassCut();
    virtual bool passes(const uhh2::Event & event) override;
};

class TwoDCut : public uhh2::Selection {
 public:
  explicit TwoDCut(float min_deltaR_, float min_pTrel_): min_deltaR(min_deltaR_), min_pTrel(min_pTrel_) {}
  virtual bool passes(const uhh2::Event & event) override;

 private:
  float min_deltaR, min_pTrel;
};

class METCut : public uhh2::Selection {
 public:
  explicit METCut(float minMet_, float maxMet_=infinity, bool genMet_=false);
  virtual bool passes(const uhh2::Event & event) override;

 private:
  float minMet, maxMet;
  bool genMet;
};

class NMuonBTagSelection: public uhh2::Selection {
  public:

    explicit NMuonBTagSelection(int min_nbtag, int max_nbtag=999, JetId btag=CSVBTag(CSVBTag::WP_LOOSE), double ptmin=0., double etamax=infinity );

    virtual bool passes(const uhh2::Event &) override;

  private:
    int m_min_nbtag;
    int m_max_nbtag;
    JetId m_btag;
    double m_ptmin;
    double m_etamax;
};

class HTCut: public uhh2::Selection {
 public:

  explicit HTCut(uhh2::Context & ctx, float min_ht_=-1,float max_ht_=infinity, std::string ht_handlename="HT");

  virtual bool passes(const uhh2::Event &) override;

 private:
  float min_ht, max_ht;
  uhh2::Event::Handle<double> h_ht;
};

class WToMuNuSelection: public uhh2::Selection {
 public:

  explicit WToMuNuSelection(float min_pt_=-1,float max_pt_=infinity);

  virtual bool passes(const uhh2::Event &) override;

 private:
  float min_pt, max_pt;
};

// //adapted from https://github.com/UHH2/VLQSemiLepPreSel/blob/master/include/VLQCommonModules.h#L573
class GenParticlePDGIdId{
 public:
  GenParticlePDGIdId(const int pdgId, const int status = -1):pdgId_(pdgId),status_(status) {}
  bool operator()(const GenParticle & p, const Event &) const { return (p.pdgId() == pdgId_) && (p.status()==status_);}
 private:
  int pdgId_;
  int status_;
};

class TTbarGenSemilepSelection: public uhh2::Selection {
 public:

  explicit TTbarGenSemilepSelection(uhh2::Context & ctx, std::string ttgen_handlename, float mu_pt_min = -1., float mu_pt_max=infinity);

  virtual bool passes(const uhh2::Event &) override;

 private:
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  float mu_pt_min,mu_pt_max;
};


