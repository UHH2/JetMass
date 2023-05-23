#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

using namespace uhh2;

std::string UHH2_Release() { return "RunII_106X_v1";};


MuonId muid = AndId<Muon>(MuonID(Muon::CutBasedIdTight), PtEtaCut(55., 2.4));
ElectronId eleid = AndId<Electron>(ElectronTagID(Electron::tagname2tag("cutBasedElectronID-Fall17-94X-V2-medium")), PtEtaCut(55., 2.4));

ElectronId eleid_veto = AndId<Electron>(ElectronTagID(Electron::tagname2tag("cutBasedElectronID-Fall17-94X-V2-veto")), PtEtaCut(10.,2.4));
ElectronId eleid_loose = AndId<Electron>(ElectronTagID(Electron::tagname2tag("cutBasedElectronID-Fall17-94X-V2-loose")), PtEtaCut(10.,2.4));

MuonId muid_loose = AndId<Muon>(MuonID(Muon::CutBasedIdLoose), PtEtaCut(10., 2.4));
