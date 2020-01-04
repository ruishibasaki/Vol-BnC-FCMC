##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of Osi 							     #
#									     #
#
##############################################################################

libOsi_suffix = \
	OsiAuxInfo  \
	OsiBranchingObject  \
	OsiChooseVariable  \
	OsiColCut  \
	OsiCut  \
	OsiCuts  \
	OsiNames \
	OsiPresolve  \
	OsiRowCut  \
	OsiRowCutDebugger  \
	OsiSolverBranch  \
	OsiSolverInterface  

includeOSI_HEADERS = \
	OsiAuxInfo.hpp \
	OsiBranchingObject.hpp \
	OsiChooseVariable.hpp \
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.hpp \
	OsiCuts.hpp \
	OsiPresolve.hpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.hpp \
	OsiSolverBranch.hpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp \
	
	

LINK = -Wl,-rpath,../../lib -L../../lib -lCoinUtils

COININC = -I../CoinUtils/include/ -I../CoinUtils/config/
# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

#STATIC LIBRARY
#OSILIB = -L$(OSIDIR)build/ -lOsi 
#libOsi_ = $(OSIDIR)build/libOsi.a 

#SHARED LIBRARY
libOsi_ = ../../lib/libOsi.so 

OSIOBJ = $(libOsi_suffix:%=build/%.o)
OSIINC = -Iinclude/ -Iconfig/
OSIH   = include/$(includeOSI_HEADERS)

# dependencies: every .o from its .C + every recursively included .h- - - - -

all: $(libOsi_suffix) $(libOsi_)  

$(libOsi_suffix):
	$(CC) -fPIC -c src/$@.cpp -o build/$@.o  $(COININC) $(OSIINC) -O3

#SHARED LIBRARY
$(libOsi_):
	g++ -shared  -Wl,-install_name,/Users/ruishibasaki/Documents/thesis/development/mcnd/bcp/lib/libOsi.so -o $@ $(LINK) $(OSIOBJ)
	
# STATIC LIBRARY
#$(libOsi_): $(OSIOBJ)
#	ar rcs $@ $+
#	ranlib $@
#	rm -f $(OSIDIR)*.o $(OSIDIR)*~

clean::
	rm -f build/*.o build/*~ ../../lib/libOsi.so


############################ End of makefile #################################
