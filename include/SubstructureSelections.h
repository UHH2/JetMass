#pragma once
#ifndef UHH2_SUBSTRUCTURESEL_INCL
#define UHH2_SUBSTRUCTURESEL_INCL

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "TH2F.h"
#include "UHH2/core/include/TopJet.h"

namespace uhh2 {
  double computeDDTValue(TopJet ak8jet,TH2F* ddtMap,bool useRho=true);
  double computeTau21DDT(TopJet ak8jet, double slope);

  class N2ddtSelection: public uhh2::Selection {
  	public:
  		N2ddtSelection(TH2F* ddtMap,bool useRho=false);
  		virtual bool passes(const uhh2::Event &event) override;
  	private:
      TH2F* ddtMap;
      bool useRho;
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
}
#endif
