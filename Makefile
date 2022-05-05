LIBRARY := SUHH2JetMass
LHAPDFINC=$(shell scram tool tag lhapdf INCLUDE)
LHAPDFLIB=$(shell scram tool tag LHAPDF LIBDIR)
FJINC=$(shell scram tool tag FASTJET INCLUDE) #added
FJLIB=$(shell scram tool tag FASTJET LIBDIR) #added
USERCXXFLAGS := -I${LHAPDFINC}
USERLDFLAGS := -lUnfold -lSUHH2core -lSUHH2common -lGenVector -lSUHH2JetMETObjects -L${LHAPDFLIB} -lLHAPDF

USERCXXFLAGS += -I${FJINC} #added
USERLDFLAGS += -L${FJLIB} -lfastjettools -lfastjet #added 

USERCXXFLAGS += -D$$CMSSW_VERSION

# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
