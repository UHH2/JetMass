#pragma once
#ifndef UHH2_SUBSTRUCTURESEL_INCL
#define UHH2_SUBSTRUCTURESEL_INCL

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "TH2F.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/PFParticle.h"

namespace uhh2 {
  double computeDDTValue(TopJet ak8jet,TH2F* ddtMap,bool useRho=true,std::vector<PFParticle> particles={});
  double computeTau21DDT(TopJet ak8jet, double slope);

  class N2ddtSelection: public uhh2::Selection {
  public:
    N2ddtSelection(TH2F* ddtMap,bool useRho=false,bool usePFMass=false);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    TH2F* ddtMap;
    bool useRho;
    bool usePFMass;
  };

  class DeepBoosted_WvsQCD_ddtSelection: public uhh2::Selection {
  public:
    DeepBoosted_WvsQCD_ddtSelection(TH2F* ddtMap,bool useRho=false,bool usePFMass=false);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    TH2F* ddtMap;
    bool useRho;
    bool usePFMass;
  };
  
  class N2Selection: public uhh2::Selection {
  public:
    N2Selection(double n2_cutValue,double n2_minValue=0.0);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double n2_cutValue,n2_minValue;
  };

  class Tau21DDTSelection: public uhh2::Selection {
  public:
    Tau21DDTSelection(double cutValue,double slope);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double cutValue,slope;
  };
    
  class RhoSelection: public uhh2::Selection {
  public:
    RhoSelection(double rho_min, double rho_max);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double rho_min, rho_max;    
  };
}
#endif
