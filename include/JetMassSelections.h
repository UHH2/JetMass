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
