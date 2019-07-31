#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include <vector>

class MassCut: public uhh2::Selection {
public:
    MassCut();
    virtual bool passes(const uhh2::Event & event) override;
};
