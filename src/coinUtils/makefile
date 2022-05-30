##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of CoinUtils 							     #
#	
#
##############################################################################

libCoinUtils_suffix = \
	CoinAlloc   \
	CoinBuild   \
	CoinDenseVector   \
	CoinError   \
	CoinFactorization1 \
	CoinFactorization2 \
	CoinFactorization3 \
	CoinFactorization4 \
	CoinSimpFactorization \
	CoinDenseFactorization \
	CoinOslFactorization \
	CoinOslFactorization2 \
	CoinOslFactorization3 \
	CoinFileIO   \
	CoinFinite   \
	CoinIndexedVector   \
	CoinLpIO   \
	CoinMessage   \
	CoinMessageHandler   \
	CoinModel   \
	CoinStructuredModel   \
	CoinModelUseful   \
	CoinModelUseful2 \
	CoinMpsIO   \
	CoinPackedMatrix   \
	CoinPackedVector   \
	CoinPackedVectorBase   \
	CoinParam CoinParamUtils   \
	CoinPostsolveMatrix \
	CoinPrePostsolveMatrix \
	CoinPresolveDoubleton   \
	CoinPresolveDual   \
	CoinPresolveDupcol   \
	CoinPresolveEmpty   \
	CoinPresolveFixed   \
	CoinPresolveForcing   \
	CoinPresolveHelperFunctions \
	CoinPresolveImpliedFree   \
	CoinPresolveIsolated   \
	CoinPresolveMatrix   \
	CoinPresolvePsdebug   \
	CoinPresolveMonitor   \
	CoinPresolveSingleton   \
	CoinPresolveSubst   \
	CoinPresolveTighten   \
	CoinPresolveTripleton   \
	CoinPresolveUseless   \
	CoinPresolveZeros   \
	CoinRational   \
	CoinSearchTree   \
	CoinShallowPackedVector   \
	CoinSnapshot   \
	CoinWarmStartBasis   \
	CoinWarmStartVector   \
	CoinWarmStartDual   \
	CoinWarmStartPrimalDual 

includecoin_HEADERS = \
	CoinOslC.h \
	CoinAlloc.hpp \
	CoinBuild.hpp \
	CoinDenseVector.hpp \
	CoinDistance.hpp \
	CoinError.hpp \
	CoinFactorization.hpp \
	CoinSimpFactorization.hpp \
	CoinDenseFactorization.hpp \
	CoinOslFactorization.hpp \
	CoinFileIO.hpp \
	CoinFinite.hpp \
	CoinFloatEqual.hpp \
	CoinHelperFunctions.hpp \
	CoinIndexedVector.hpp \
	CoinLpIO.hpp \
	CoinMessage.hpp \
	CoinMessageHandler.hpp \
	CoinModel.hpp \
	CoinStructuredModel.hpp \
	CoinModelUseful.hpp \
	CoinMpsIO.hpp \
	CoinPackedMatrix.hpp \
	CoinPackedVector.hpp \
	CoinPackedVectorBase.hpp \
	CoinParam.hpp \
	CoinPragma.hpp \
	CoinPresolveDoubleton.hpp \
	CoinPresolveDual.hpp \
	CoinPresolveDupcol.hpp \
	CoinPresolveEmpty.hpp \
	CoinPresolveFixed.hpp \
	CoinPresolveForcing.hpp \
	CoinPresolveImpliedFree.hpp \
	CoinPresolveIsolated.hpp \
	CoinPresolveMatrix.hpp \
	CoinPresolveMonitor.hpp \
	CoinPresolvePsdebug.hpp \
	CoinPresolveSingleton.hpp \
	CoinPresolveSubst.hpp \
	CoinPresolveTighten.hpp \
	CoinPresolveTripleton.hpp \
	CoinPresolveUseless.hpp \
	CoinPresolveZeros.hpp \
	CoinRational.hpp \
	CoinSearchTree.hpp \
	CoinShallowPackedVector.hpp \
	CoinSignal.hpp \
	CoinSmartPtr.hpp \
	CoinSnapshot.hpp \
	CoinSort.hpp \
	CoinTime.hpp \
	CoinTypes.hpp \
	CoinUtility.hpp \
	CoinWarmStart.hpp \
	CoinWarmStartBasis.hpp \
	CoinWarmStartVector.hpp \
	CoinWarmStartDual.hpp \
	CoinWarmStartPrimalDual.hpp

	


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

#STATIC LIBRARY
#COINLIB = -L$(COINDIR)build/ -lCoinUtils 
#libCoinUtils_ = $(COINDIR)build/libCoinUtils.a 

#SHARED LIBRARY
libCoinUtils_ = ../../lib/libCoinUtils.so

COINOBJ = $(libCoinUtils_suffix:%=build/%.o)
COININC = -Iinclude/ -Iconfig/
COINH   = include/$(includecoin_HEADERS)

# dependencies: every .o from its .C + every recursively included .h- - - - -

all: $(libCoinUtils_suffix) $(libCoinUtils_) 

$(libCoinUtils_suffix):
	g++ -fPIC -c ./src/$@.cpp -o build/$@.o $(COININC) -O3 

#SHARED LIBRARY
#MAC
#$(libCoinUtils_):
#	g++ -shared  -Wl,-install_name,/Users/ruishibasaki/Documents/thesis/development/mcnd/bcp/lib/libCoinUtils.so -o $@ $(COINOBJ) -lbz2 -lz -lm -llapack -O3

#LIMOS
$(libCoinUtils_):
	g++ -shared  -Wl,-soname,/home/shibasaki/development/bcp/Vol-BnC-FCMC-master/lib/libCoinUtils.so -o $@ $(COINOBJ) -lbz2 -lz -llapack -lblas

# STATIC LIBRARY
#$(libCoinUtils_): $(COINOBJ)
#	g++ -shared

clean::
	rm -f build/*.o build/*~  ../../lib/libCoinUtils.so

############################ End of makefile #################################
