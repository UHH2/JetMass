#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"
#include <vector>


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
  explicit METCut(float minMet_, float maxMet_=infinity);
  virtual bool passes(const uhh2::Event & event) override;

 private:
  float minMet, maxMet;
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

  explicit HTCut(uhh2::Context & ctx, float min_ht_=-1,float max_ht_=infinity);

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
