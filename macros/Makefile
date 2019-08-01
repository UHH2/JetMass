CC = g++
ROOTFLAGS = `root-config --cflags --libs` -L $ROOTSYS/lib -lRooFit -lRooFitCore

INCDIR   = include
SRCDIR   = src

SOURCES  := $(wildcard $(SRCDIR)/*.cc)
INCLUDES := $(wildcard $(INCDIR)/*.h)

MACROS_temp := $(basename $(SOURCES))
MACROS := $(notdir $(MACROS_temp) )

all: $(MACROS)


$(MACROS): % : $(SRCDIR)/%.cc $(INCLUDES)
	@echo 'building' $@
	@$(CC) -o $@ $(SRCDIR)/$@.cc $(ROOTFLAGS)


clean:
	@rm -f $(MACROS)
