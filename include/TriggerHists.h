#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


using namespace std;

class TriggerHists: public uhh2::Hists {
public:
  TriggerHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual ~TriggerHists();

private:
  bool isMC;
  uhh2::Event::Handle<double> h_ht;
};
