CXX = g++
CXXFLAGS = -g -O2 -Wall -fPIC

# --- ROOT --------------------------------------------------------------
CXXFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lThread -lSpectrum

# --- EXOSoftware -------------------------------------------------------
CXXFLAGS +=-I/nfs/slac/g/exo/software/hudson/builds-rhel5/current/include
EXOLIBS = -L/nfs/slac/g/exo/software/hudson/builds-rhel5/current/lib -lEXOUtilities

all: DoubleEscape

DoubleEscape: DoubleEscape.o
	@ echo "Linking $@..."
	@ ${CXX} ${CXXFLAGS} -o $@ $^ ${ROOTLIBS} ${EXOLIBS}

EnergyCalib: EnergyCalib.o
	@ echo "Linking $@..."
	@ ${CXX} ${CXXFLAGS} -o $@ $^ ${ROOTLIBS} ${EXOLIBS}

MinimizePlot: MinimizePlot.o
	@ echo "Linking $@..."
	@ ${CXX} ${CXXFLAGS} -o $@ $^ ${ROOTLIBS} ${EXOLIBS}

DoubleEscape.o: DoubleEscape.cc
	@ ${CXX} ${CXXFLAGS} -c $^ -o $@

EnergyCalib.o: EnergyCalib.cc
	@ ${CXX} ${CXXFLAGS} -c $^ -o $@

MinimizePlot.o: MinimizePlot.cc
	@ ${CXX} ${CXXFLAGS} -c $^ -o $@

.PHONY : clean

clean:
	@echo "Cleaning up..."
	@ rm -f DoubleEscape DoubleEscape.o EnergyCalib EnergyCalib.o MinimizePlot MinimizePlot.o
