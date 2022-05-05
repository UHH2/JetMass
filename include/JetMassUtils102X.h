#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

using namespace uhh2;

std::string UHH2_Release() { return "RunII_102X_v2";};


MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
ElectronId eleid = AndId<Electron>(ElectronID_Fall17_medium_noIso, PtEtaCut(55., 2.4));

