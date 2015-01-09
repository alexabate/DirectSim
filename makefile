################################################################################
##  makefile                                                                  ##
################################################################################

include $(SOPHYABASE)/include/sophyamake.inc

#### Automatically set up GALSIM variable if it doesn't exist
ifndef GALSIM
	GALSIM := ${PWD}
	echo "GALSIM = ${GALSIM}"
endif


#define GEN3D_FLOAT

OBJ = ${GALSIM}/objs

EXE = ${GALSIM}/exe

ROOTOUT = ${GALSIM}/root

TESTS =${GALSIM}/testfiles

MYCL = ${GALSIM}/classes

PROGS = ${GALSIM}/progs

KCORR = ${GALSIM}/kCorrections

LIBH := $(MYCL)/cosmocalcs.h $(MYCL)/geneutils.h $(MYCL)/gftdist.h \
$(MYCL)/schechter.h $(MYCL)/sinterp.h $(MYCL)/simdata.h $(MYCL)/reddening.h \
$(MYCL)/sedfilter.h $(MYCL)/genefluct3d.h  $(MYCL)/pkspectrum.h \
$(MYCL)/mass2gal.h $(MYCL)/powerspec.h $(MYCL)/matrix.h $(MYCL)/igm.h \
$(MYCL)/hpoly.h $(MYCL)/shapelets.h $(MYCL)/em.h $(MYCL)/igmstatistics.h
#$(MYCL)/constcosmo.h
#$(MYCL)/root_plots.h

LIBO := $(OBJ)/cosmocalcs.o $(OBJ)/geneutils.o $(OBJ)/gftdist.o \
$(OBJ)/schechter.o $(OBJ)/sinterp.o $(OBJ)/simdata.o $(OBJ)/reddening.o \
$(OBJ)/sedfilter.o $(OBJ)/genefluct3d.o  $(OBJ)/pkspectrum.o $(OBJ)/mass2gal.o \
$(OBJ)/powerspec.o $(OBJ)/matrix.o $(OBJ)/igm.o $(OBJ)/hpoly.o $(OBJ)/shapelets.o\
$(OBJ)/em.o $(OBJ)/igmstatistics.o
#$(OBJ)/root_plots.o

# root libraries
ROOTLIB   = $(shell root-config --libs)
ROOTINC = $(shell root-config --incdir)

MYLIB = -lfftw3_threads
MINUIT = -lMinuit


################################################################################

all : progs tests test

progs : analyzeBPZ baseSimulation calculateKcorrections colorDistributions \
convertSEDS fitLSSTspectra lineOfSightLymanAlpha lineOfSightMagnitude \
pcaTemplates photoZdist priorFitter projectTemplates rdlss sdssElColors simdensity \
simulateAbsorberLinesOfSight simulateLSSTobsFromTruth simulateLSSTobs 


tests : test2Dinterp testbasesim testErrors testEMalgorithm testgoodsmagsim    \
testKcorrColors testKcorrMethod testLF testMadau testMeiksin \
 testTemplateFitting testsimdensity 
#testSimReadKcorr testsimulateIGM testSimulation

clean : 
	rm  $(OBJ)/* $(EXE)/*
	
# MAIN PROGS

analyzeBPZ : $(EXE)/analyzeBPZ
	@echo 'makefile : analyzeBPZ made'

baseSimulation : $(EXE)/baseSimulation
	@echo 'makefile : baseSimulation made'
	
calculateKcorrections	: $(EXE)/calculateKcorrections
	@echo 'makefile : calculateKcorrections made'
	
colorDistributions	: $(EXE)/colorDistributions
	@echo 'makefile : colorDistributions made'

convertSEDS	: $(EXE)/convertSEDS
	@echo 'makefile : convertSEDS made'

fitLSSTspectra : $(EXE)/fitLSSTspectra
	@echo 'makefile : fitLSSTspectra made'
	
lineOfSightLymanAlpha : $(EXE)/lineOfSightLymanAlpha
	@echo 'makefile : lineOfSightLymanAlpha made'

lineOfSightMagnitude : $(EXE)/lineOfSightMagnitude
	@echo 'makefile : lineOfSightMagnitude made'
	
pcaTemplates : $(EXE)/pcaTemplates
	@echo 'makefile : pcaTemplates made'
	
photoZdist : $(EXE)/photoZdist
	@echo 'makefile : photoZdist made'

priorFitter : $(EXE)/priorFitter
	@echo 'makefile : priorFitter made'
	
projectTemplates : $(EXE)/projectTemplates
	@echo 'makefile : projectTemplates made'
	
rdlss : $(EXE)/rdlss
	@echo 'makefile : rdlss made'
	
sdssElColors : $(EXE)/sdssElColors
	@echo 'makefile : sdssElColors made'
	
simdensity : $(EXE)/simdensity
	@echo 'makefile : simdensity made'

simulateAbsorberLinesOfSight : $(EXE)/simulateAbsorberLinesOfSight
	@echo 'makefile : simulateAbsorberLinesOfSight made'
	
simulateLSSTobs : $(EXE)/simulateLSSTobs
	@echo 'makefile : simulateLSSTobs made'	
	
simulateLSSTobsFromTruth : $(EXE)/simulateLSSTobsFromTruth
	@echo 'makefile : simulateLSSTobsFromTruth made'


# TESTING PROGS

test2Dinterp : $(EXE)/test2Dinterp 
	@echo 'makefile :test2Dinterp made'
	
testbasesim : $(EXE)/testbasesim 
	@echo 'makefile :testbasesim made'
	
test : $(EXE)/test
	@echo 'makefile : test made'
	
testEMalgorithm : $(EXE)/testEMalgorithm
	@echo 'makefile : testEMalgorithm made'
	
testErrors : $(EXE)/testErrors
	@echo 'makefile : testErrors made'
	
testgoodsmagsim : $(EXE)/testgoodsmagsim
	@echo 'makefile : testgoodsmagsim made'

testKcorrColors : $(EXE)/testKcorrColors 
	@echo 'makefile : testKcorrColors made'
	
testKcorrMethod : $(EXE)/testKcorrMethod 
	@echo 'makefile : testKcorrMethod made'

testLF : $(EXE)/testLF 
	@echo 'makefile :testLF made'

testMadau : $(EXE)/testMadau
	@echo 'makefile : testMadau made'

testMeiksin : $(EXE)/testMeiksin
	@echo 'makefile : testMeiksin made'
	
testSimReadKcorr : $(EXE)/testSimReadKcorr
	@echo 'makefile : testSimReadKcorr made'
	
testsimulateIGM : $(EXE)/testsimulateIGM
	@echo 'makefile : testsimulateIGM made'
	
testSimulation : $(EXE)/testSimulation
	@echo 'makefile : testSimulation made'
	
testTemplateFitting : $(EXE)/testTemplateFitting
	@echo 'makefile : testTemplateFitting made'

## programs below here have not been CHECKED or maybe even finished...

lsstSbIgmColors : $(EXE)/lsstSbIgmColors
	@echo 'makefile : lsstSbIgmColors made'

testpsdenscube : $(EXE)/testpsdenscube
	@echo 'makefile :  testpsdenscube made'

testsimdensity : $(EXE)/testsimdensity
	@echo 'makefile :  testsimdensity made'
	
###################### MAIN PROGRAMS ###########################################

# ANALYZE A BPZ CATALOG
$(EXE)/analyzeBPZ : $(OBJ)/analyzeBPZ.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/analyzeBPZ $(OBJ)/analyzeBPZ.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/analyzeBPZ.o : $(PROGS)/analyzeBPZ.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/analyzeBPZ.o $(PROGS)/analyzeBPZ.cc 
	
# BASE SIMULATION
$(EXE)/baseSimulation : $(OBJ)/baseSimulation.o $(LIBO) 
	mkdir -p $(EXE)
	$(CXXLINK) -o $(EXE)/baseSimulation $(OBJ)/baseSimulation.o $(LIBO) \
	$(SOPHYAEXTSLBLIST)	$(MYLIB) $(ROOTLIB)

$(OBJ)/baseSimulation.o : $(PROGS)/baseSimulation.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/baseSimulation.o $(PROGS)/baseSimulation.cc
	
# CALCULATE K CORRECTION TABLES
$(EXE)/calculateKcorrections : $(OBJ)/calculateKcorrections.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/calculateKcorrections $(OBJ)/calculateKcorrections.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/calculateKcorrections.o : $(PROGS)/calculateKcorrections.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/calculateKcorrections.o $(PROGS)/calculateKcorrections.cc
	
# COLOR DISTRIBUTIONS
$(EXE)/colorDistributions : $(OBJ)/colorDistributions.o $(LIBO) 
	mkdir -p $(EXE)
	$(CXXLINK) -o $(EXE)/colorDistributions $(OBJ)/colorDistributions.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/colorDistributions.o : $(PROGS)/colorDistributions.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/colorDistributions.o $(PROGS)/colorDistributions.cc

# CONVERT SED UNITS
$(EXE)/convertSEDS : $(OBJ)/convertSEDS.o
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/convertSEDS $(OBJ)/convertSEDS.o $(SOPHYAEXTSLBLIST)

$(OBJ)/convertSEDS.o : $(PROGS)/convertSEDS.cc
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -o $(OBJ)/convertSEDS.o $(PROGS)/convertSEDS.cc 
	
# FIT LSST SPECTRA TO CWWK
$(EXE)/fitLSSTspectra : $(OBJ)/fitLSSTspectra.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitLSSTspectra $(OBJ)/fitLSSTspectra.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitLSSTspectra.o : $(PROGS)/fitLSSTspectra.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitLSSTspectra.o $(PROGS)/fitLSSTspectra.cc
	
# SIMULATE LINE OF SIGHT LYMAN ALPHA TRANSMISSION
$(EXE)/lineOfSightLymanAlpha : $(OBJ)/lineOfSightLymanAlpha.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lineOfSightLymanAlpha $(OBJ)/lineOfSightLymanAlpha.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightLymanAlpha.o : $(PROGS)/lineOfSightLymanAlpha.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lineOfSightLymanAlpha.o $(PROGS)/lineOfSightLymanAlpha.cc
	
# U,G MAGNITUDES WITH IGM LINE OF SIGHT
$(EXE)/lineOfSightMagnitude : $(OBJ)/lineOfSightMagnitude.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lineOfSightMagnitude $(OBJ)/lineOfSightMagnitude.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightMagnitude.o : $(PROGS)/lineOfSightMagnitude.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lineOfSightMagnitude.o $(PROGS)/lineOfSightMagnitude.cc
	
# TEMPLATE PCA
$(EXE)/pcaTemplates : $(OBJ)/pcaTemplates.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/pcaTemplates $(OBJ)/pcaTemplates.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/pcaTemplates.o : $(PROGS)/pcaTemplates.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/pcaTemplates.o $(PROGS)/pcaTemplates.cc 
	
# PHOTO-Z DISTRIBUTION
$(EXE)/photoZdist : $(OBJ)/photoZdist.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/photoZdist $(OBJ)/photoZdist.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/photoZdist.o : $(PROGS)/photoZdist.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/photoZdist.o $(PROGS)/photoZdist.cc 
	
# PRIOR FITTER
$(EXE)/priorFitter : $(OBJ)/priorFitter.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/priorFitter $(OBJ)/priorFitter.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(MINUIT)

$(OBJ)/priorFitter.o : $(PROGS)/priorFitter.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/priorFitter.o $(PROGS)/priorFitter.cc 
	
# PROJECT TEMPLATES
$(EXE)/projectTemplates : $(OBJ)/projectTemplates.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/projectTemplates $(OBJ)/projectTemplates.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/projectTemplates.o : $(PROGS)/projectTemplates.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/projectTemplates.o $(PROGS)/projectTemplates.cc 
	
# SIMULATE CATALOG OF BASIC GALAXY PROPERTIES FROM OVERDENSITY GRID
$(EXE)/rdlss : $(OBJ)/rdlss.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/rdlss $(OBJ)/rdlss.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/rdlss.o : $(PROGS)/rdlss.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/rdlss.o $(PROGS)/rdlss.cc 
	
# CALCULATE SDSS COLORS OF ELLIPTICAL GALAXY
$(EXE)/sdssElColors : $(OBJ)/sdssElColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sdssElColors $(OBJ)/sdssElColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssElColors.o : $(PROGS)/sdssElColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sdssElColors.o $(PROGS)/sdssElColors.cc 
	
# SIMULATE OVERDENSITY GRID
$(EXE)/simdensity : $(OBJ)/simdensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simdensity $(OBJ)/simdensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simdensity.o : $(PROGS)/simdensity.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simdensity.o $(PROGS)/simdensity.cc 
	
# SIMULATE LINE OF SIGHT ABSORBER DISTRIBUTIONS 
$(EXE)/simulateAbsorberLinesOfSight : $(OBJ)/simulateAbsorberLinesOfSight.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateAbsorberLinesOfSight $(OBJ)/simulateAbsorberLinesOfSight.o \
	$(LIBO) $(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateAbsorberLinesOfSight.o : $(PROGS)/simulateAbsorberLinesOfSight.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateAbsorberLinesOfSight.o \
	$(PROGS)/simulateAbsorberLinesOfSight.cc 
	
# SIMULATE LSST OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateLSSTobs : $(OBJ)/simulateLSSTobs.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateLSSTobs $(OBJ)/simulateLSSTobs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobs.o : $(PROGS)/simulateLSSTobs.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateLSSTobs.o \
	$(PROGS)/simulateLSSTobs.cc 
	
# SIMULATE LSST OBSERVATIONS FROM INPUT IMSIM CATALOG OF TRUE MAGS AND REDSHIFTS
$(EXE)/simulateLSSTobsFromTruth : $(OBJ)/simulateLSSTobsFromTruth.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateLSSTobsFromTruth $(OBJ)/simulateLSSTobsFromTruth.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobsFromTruth.o : $(PROGS)/simulateLSSTobsFromTruth.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateLSSTobsFromTruth.o \
	$(PROGS)/simulateLSSTobsFromTruth.cc 
	
###################### INCOMPLETE PROGRAMS ####################################


$(EXE)/lsstSbIgmColors : $(OBJ)/lsstSbIgmColors.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/lsstSbIgmColors $(OBJ)/lsstSbIgmColors.o \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(LIBO)

$(OBJ)/lsstSbIgmColors.o : $(PROGS)/lsstSbIgmColors.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lsstSbIgmColors.o $(PROGS)/lsstSbIgmColors.cc 


	
###################### TESTING PROGRAMS ########################################

# TEST 2D INTERPOLATION
$(EXE)/test2Dinterp : $(OBJ)/test2Dinterp.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/test2Dinterp $(OBJ)/test2Dinterp.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/test2Dinterp.o : $(PROGS)/test2Dinterp.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/test2Dinterp.o $(PROGS)/test2Dinterp.cc

# TEST BASESIM
$(EXE)/testbasesim : $(OBJ)/testbasesim.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testbasesim $(OBJ)/testbasesim.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testbasesim.o : $(PROGS)/testbasesim.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testbasesim.o $(PROGS)/testbasesim.cc
	
# TEST EXPECTATION-MAXIMIZATION ALGORITHM
$(EXE)/testEMalgorithm : $(OBJ)/testEMalgorithm.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testEMalgorithm $(OBJ)/testEMalgorithm.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testEMalgorithm.o : $(PROGS)/testEMalgorithm.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testEMalgorithm.o $(PROGS)/testEMalgorithm.cc
	
# TEST LSST ERRORS
$(EXE)/testErrors : $(OBJ)/testErrors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testErrors $(OBJ)/testErrors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testErrors.o : $(PROGS)/testErrors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testErrors.o $(PROGS)/testErrors.cc
	
# TEST GOODS SIM
$(EXE)/testgoodsmagsim : $(OBJ)/testgoodsmagsim.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testgoodsmagsim $(OBJ)/testgoodsmagsim.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testgoodsmagsim.o : $(PROGS)/testgoodsmagsim.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testgoodsmagsim.o $(PROGS)/testgoodsmagsim.cc

# TEST K CORRECTION 
$(EXE)/testKcorrColors : $(OBJ)/testKcorrColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testKcorrColors $(OBJ)/testKcorrColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrColors.o : $(PROGS)/testKcorrColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testKcorrColors.o $(PROGS)/testKcorrColors.cc
	
# TEST K CORRECTION INTERPOLATION
$(EXE)/testKcorrMethod : $(OBJ)/testKcorrMethod.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testKcorrMethod $(OBJ)/testKcorrMethod.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrMethod.o : $(PROGS)/testKcorrMethod.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testKcorrMethod.o $(PROGS)/testKcorrMethod.cc

# TEST LF 
$(EXE)/testLF : $(OBJ)/testLF.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testLF $(OBJ)/testLF.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLF.o : $(PROGS)/testLF.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/testLF.o $(PROGS)/testLF.cc 
	
# TEST MADAU
$(EXE)/testMadau : $(OBJ)/testMadau.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testMadau $(OBJ)/testMadau.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMadau.o : $(PROGS)/testMadau.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testMadau.o $(PROGS)/testMadau.cc
	
# TEST MEIKSIN
$(EXE)/testMeiksin : $(OBJ)/testMeiksin.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testMeiksin $(OBJ)/testMeiksin.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMeiksin.o : $(PROGS)/testMeiksin.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testMeiksin.o $(PROGS)/testMeiksin.cc
	
# TEST SIMULATE DATA USING K CORRECTION TABLES
$(EXE)/testSimReadKcorr : $(OBJ)/testSimReadKcorr.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/testSimReadKcorr $(OBJ)/testSimReadKcorr.o \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(LIBO)

$(OBJ)/testSimReadKcorr.o : $(PROGS)/testSimReadKcorr.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testSimReadKcorr.o $(PROGS)/testSimReadKcorr.cc 

# TEST SIMULATE IGM
$(EXE)/testsimulateIGM : $(OBJ)/testsimulateIGM.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testsimulateIGM $(OBJ)/testsimulateIGM.o \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(LIBO)

$(OBJ)/testsimulateIGM.o : $(PROGS)/testsimulateIGM.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testsimulateIGM.o $(PROGS)/testsimulateIGM.cc 
	
# TEST SIMULATION
$(EXE)/testSimulation : $(OBJ)/testSimulation.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/testSimulation $(OBJ)/testSimulation.o \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(LIBO)

$(OBJ)/testSimulation.o : $(PROGS)/testSimulation.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testSimulation.o $(PROGS)/testSimulation.cc
	
# TEST TEMPLATE FITTING
$(EXE)/testTemplateFitting : $(OBJ)/testTemplateFitting.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testTemplateFitting $(OBJ)/testTemplateFitting.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testTemplateFitting.o : $(PROGS)/testTemplateFitting.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testTemplateFitting.o $(PROGS)/testTemplateFitting.cc
	
	
	
# TEMPORARY TEST CODE
$(EXE)/test : $(OBJ)/test.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/test $(OBJ)/test.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/test.o : $(PROGS)/test.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/test.o $(PROGS)/test.cc

## programs below here have not been CHECKED or maybe even finished...

# TEST POWER SPECTRUM FROM OVER-DENSITY CUBE
$(EXE)/testpsdenscube : $(OBJ)/testpsdenscube.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testpsdenscube $(OBJ)/testpsdenscube.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testpsdenscube.o : $(PROGS)/testpsdenscube.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testpsdenscube.o $(PROGS)/testpsdenscube.cc 
	
# TEST DENSITY SIMULATION
$(EXE)/testsimdensity : $(OBJ)/testsimdensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testsimdensity $(OBJ)/testsimdensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testsimdensity.o : $(PROGS)/testsimdensity.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testsimdensity.o $(PROGS)/testsimdensity.cc 

############################# UTILITIES ########################################

#$(OBJ)/root_plots.o : $(MYCL)/root_plots.cpp
#	$(CXXCOMPILE) -o $(OBJ)/root_plots.o -I$(ROOTINC) -I$(MYCL) -c $(MYCL)/root_plots.cpp


################################  THE CLASSES ##################################

$(OBJ)/%.o : $(MYCL)/%.cc $(MYCL)/%.h 
	mkdir -p $(OBJ)
	mkdir -p $(TESTS)
	$(CXXCOMPILE) -I$(ROOTINC) -I$(MYCL) -o $@ $<

#$(OBJ)genefluct3df.o: $(MYCL)genefluct3d.cc $(MYCL)genefluct3d.h
#	$(CXXCOMPILE) -I$(MYCL) -DGEN3D_FLOAT -o $@ $(MYCL)genefluct3d.cc



