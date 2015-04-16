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

BAOPROGS = ${GALSIM}/baoprogs

ROOTOUT = ${GALSIM}/root

TESTS =${GALSIM}/testfiles

MYCL = ${GALSIM}/classes

PROGS = ${GALSIM}/progs

KCORR = ${GALSIM}/kCorrections

LIBH := $(MYCL)/cosmocalcs.h $(MYCL)/geneutils.h $(MYCL)/gftdist.h \
$(MYCL)/schechter.h $(MYCL)/sinterp.h $(MYCL)/simdata.h $(MYCL)/reddening.h \
$(MYCL)/sedfilter.h $(MYCL)/genefluct3d.h  $(MYCL)/pkspectrum.h \
$(MYCL)/mass2gal.h $(MYCL)/powerspec.h $(MYCL)/matrix.h $(MYCL)/igm.h \
$(MYCL)/hpoly.h $(MYCL)/shapelets.h $(MYCL)/em.h $(MYCL)/cat2grid.h $(MYCL)/igmstatistics.h \
$(MYCL)/fitkbaoscale.h $(MYCL)/chisqstats.h $(MYCL)/massfunc.h $(MYCL)/testclass.h

#$(MYCL)/constcosmo.h
#$(MYCL)/root_plots.h

LIBO := $(OBJ)/cosmocalcs.o $(OBJ)/geneutils.o $(OBJ)/gftdist.o \
$(OBJ)/schechter.o $(OBJ)/sinterp.o $(OBJ)/simdata.o $(OBJ)/reddening.o \
$(OBJ)/sedfilter.o $(OBJ)/genefluct3d.o  $(OBJ)/pkspectrum.o $(OBJ)/mass2gal.o \
$(OBJ)/powerspec.o $(OBJ)/matrix.o $(OBJ)/igm.o $(OBJ)/hpoly.o $(OBJ)/shapelets.o\
$(OBJ)/em.o $(OBJ)/igmstatistics.o $(OBJ)/cat2grid.o $(OBJ)/fitkbaoscale.o $(OBJ)/chisqstats.o \
$(OBJ)/massfunc.o $(OBJ)/testclass.o

#$(OBJ)/root_plots.o

# root libraries
ROOTLIB   = $(shell root-config --libs)
ROOTINC = $(shell root-config --incdir)

MYLIB = -lfftw3_threads
MINUIT = -lMinuit


################################################################################

all : progs tests 
#bao

progs : addIGMToSED analyzeBPZ baseSimulation calculateKcorrections cfhtColors \
colorDistributions convertSEDS fitLSSTspectra lineOfSightLymanAlpha lineOfSightMagnitude \
lsstPicklesLibrary lymanAlphaToDensity pcaTemplates photoZdist priorFitter \
projectTemplates rdlss sdssElColors sdssPicklesLibrary simdensity  \
simulateAbsorberLinesOfSight simulateLSSTobsFromTruth \
cluster_mass_function sim_with_clusters quickSimMags restFrameColors generateTemplSet fitSEDsToColors \
getSEDcolors

tests : test2Dinterp testbasesim testEMalgorithm testErrors \
testKcorrColors testLF testLymanAlphaAbs testMadau testMeiksin \
testsimulateIGM 
#testTemplateFitting -> defunct but might need to look at again
# testsimdensity 

bao : addGausszerr computepsfromarray fitkbao getpzconvf getsf grid_data \
sim_mcgrids subfromfull 

clean : 
	rm  $(OBJ)/* $(EXE)/*
	
# MAIN PROGS

addIGMToSED : $(EXE)/addIGMToSED
	@echo 'makefile : addIGMToSED made'

analyzeBPZ : $(EXE)/analyzeBPZ
	@echo 'makefile : analyzeBPZ made'

baseSimulation : $(EXE)/baseSimulation
	@echo 'makefile : baseSimulation made'
	
calculateKcorrections : $(EXE)/calculateKcorrections
	@echo 'makefile : calculateKcorrections made'
	
cfhtColors : $(EXE)/cfhtColors
	@echo 'makefile : cfhtColors made'
	
checkmagsim : $(EXE)/checkmagsim
	@echo 'makefile : checkmagsim made'
	
cluster_mass_function : $(EXE)/cluster_mass_function
	@echo 'makefile : cluster_mass_function made' 
	
colorDistributions	: $(EXE)/colorDistributions
	@echo 'makefile : colorDistributions made'

convertSEDS	: $(EXE)/convertSEDS
	@echo 'makefile : convertSEDS made'
	
findCWWinLSSTspectra : $(EXE)/findCWWinLSSTspectra
	@echo 'makefile : findCWWinLSSTspectra made'

fitLSSTspectra : $(EXE)/fitLSSTspectra
	@echo 'makefile : fitLSSTspectra made'
	
fitSEDsToColors : $(EXE)/fitSEDsToColors
	@echo 'makefile : fitSEDsToColors made'
	
generateTemplSet : $(EXE)/generateTemplSet
	@echo 'makefile : generateTemplSet made'
	
getSEDcolors : $(EXE)/getSEDcolors
	@echo 'makefile : getSEDcolors made'
	
lineOfSightLymanAlpha : $(EXE)/lineOfSightLymanAlpha
	@echo 'makefile : lineOfSightLymanAlpha made'

lineOfSightMagnitude : $(EXE)/lineOfSightMagnitude
	@echo 'makefile : lineOfSightMagnitude made'

lsstPicklesLibrary : $(EXE)/lsstPicklesLibrary
	@echo 'makefile : lsstPicklesLibrary made'
	
lymanAlphaToDensity : $(EXE)/lymanAlphaToDensity
	@echo 'makefile : lymanAlphaToDensity made'
	
pcaTemplates : $(EXE)/pcaTemplates
	@echo 'makefile : pcaTemplates made'
	
photoZdist : $(EXE)/photoZdist
	@echo 'makefile : photoZdist made'

priorFitter : $(EXE)/priorFitter
	@echo 'makefile : priorFitter made'
	
produceCatalogForBPZ : $(EXE)/produceCatalogForBPZ
	@echo 'makefile : produceCatalogForBPZ made'
	
projectTemplates : $(EXE)/projectTemplates
	@echo 'makefile : projectTemplates made'
	
quickSimMags : $(EXE)/quickSimMags
	@echo 'makefile : quickSimMags made'
	
restFrameColors : $(EXE)/restFrameColors
	@echo 'makefile : restFrameColors made'

rdlss : $(EXE)/rdlss
	@echo 'makefile : rdlss made'
	
sdssElColors : $(EXE)/sdssElColors
	@echo 'makefile : sdssElColors made'

sdssPicklesLibrary : $(EXE)/sdssPicklesLibrary
	@echo 'makefile : sdssPicklesLibrary made'

sim_with_clusters : $(EXE)/sim_with_clusters
	@echo 'makefile : sim_with_clusters made'
	
simdensity : $(EXE)/simdensity
	@echo 'makefile : simdensity made'

simulateAbsorberLinesOfSight : $(EXE)/simulateAbsorberLinesOfSight
	@echo 'makefile : simulateAbsorberLinesOfSight made'

simulateCFHTobs : $(EXE)/simulateCFHTobs
	@echo 'makefile : simulateCFHTobs made'

# DEFUNCT	
#simulateLSSTobs : $(EXE)/simulateLSSTobs
#	@echo 'makefile : simulateLSSTobs made'	
	
simulateLSSTobsFromTruth : $(EXE)/simulateLSSTobsFromTruth
	@echo 'makefile : simulateLSSTobsFromTruth made'

# BAO PROGS

addGausszerr : $(EXE)/addGausszerr
	@echo 'makefile : addGausszerr made'
	
computepsfromarray : $(EXE)/computepsfromarray
	@echo 'makefile : computepsfromarray made'

fitkbao : $(EXE)/fitkbao
	@echo 'makefile : fitkbao made'
	
getpzconvf : $(EXE)/getpzconvf
	@echo 'makefile : getpzconvf made'
	
getsf : $(EXE)/getsf
	@echo 'makefile : getsf made'
	
grid_data : $(EXE)/grid_data
	@echo 'makefile : grid_data made'
	
sim_mcgrids : $(EXE)/sim_mcgrids
	@echo 'makefile : sim_mcgrids made'
	
subfromfull : $(EXE)/subfromfull
	@echo 'makefile : subfromfull made'

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

# DEFUNCT	
#testgoodsmagsim : $(EXE)/testgoodsmagsim
#	@echo 'makefile : testgoodsmagsim made'

testKcorrColors : $(EXE)/testKcorrColors 
	@echo 'makefile : testKcorrColors made'

testLF : $(EXE)/testLF 
	@echo 'makefile :testLF made'
	
testLymanAlphaAbs : $(EXE)/testLymanAlphaAbs 
	@echo 'makefile :testLymanAlphaAbs made'

testMadau : $(EXE)/testMadau
	@echo 'makefile : testMadau made'

testMeiksin : $(EXE)/testMeiksin
	@echo 'makefile : testMeiksin made'

# DEFUNCT	
#testSimReadKcorr : $(EXE)/testSimReadKcorr
#	@echo 'makefile : testSimReadKcorr made'
	
testsimulateIGM : $(EXE)/testsimulateIGM
	@echo 'makefile : testsimulateIGM made'

# DEFUNCT	
#testSimulation : $(EXE)/testSimulation
#	@echo 'makefile : testSimulation made'
	
testTemplateFitting : $(EXE)/testTemplateFitting
	@echo 'makefile : testTemplateFitting made'
	
time_test_getmag : $(EXE)/time_test_getmag
	@echo 'makefile : time_test_getmag made'

## programs below here have not been CHECKED or maybe even finished...

lsstSbIgmColors : $(EXE)/lsstSbIgmColors
	@echo 'makefile : lsstSbIgmColors made'

testpsdenscube : $(EXE)/testpsdenscube
	@echo 'makefile :  testpsdenscube made'

testsimdensity : $(EXE)/testsimdensity
	@echo 'makefile :  testsimdensity made'
	
###################### MAIN PROGRAMS ###########################################

# ADD LINE OF SIGHT TRANSMISSON TO SEDS IN A LIBRARY
$(EXE)/addIGMToSED : $(OBJ)/addIGMToSED.o $(LIBO) 
	mkdir -p $(EXE)
	$(CXXLINK) -o $(EXE)/addIGMToSED $(OBJ)/addIGMToSED.o $(LIBO) \
	$(SOPHYAEXTSLBLIST)	$(MYLIB) $(ROOTLIB)

$(OBJ)/addIGMToSED.o : $(PROGS)/addIGMToSED.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/addIGMToSED.o $(PROGS)/addIGMToSED.cc

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
	
# CALCULATE U-G, I-Z CFHT colors w/ and wo/ host galaxy reddening
$(EXE)/cfhtColors : $(OBJ)/cfhtColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/cfhtColors $(OBJ)/cfhtColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/cfhtColors.o : $(PROGS)/cfhtColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/cfhtColors.o $(PROGS)/cfhtColors.cc
	
# CHECK MAG SIM
$(EXE)/checkmagsim : $(OBJ)/checkmagsim.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/checkmagsim $(OBJ)/checkmagsim.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/checkmagsim.o : $(PROGS)/checkmagsim.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/checkmagsim.o $(PROGS)/checkmagsim.cc
	
# CLUSTER MASS FUNCTION
$(EXE)/cluster_mass_function : $(OBJ)/cluster_mass_function.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/cluster_mass_function $(OBJ)/cluster_mass_function.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/cluster_mass_function.o : $(PROGS)/cluster_mass_function.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/cluster_mass_function.o \
	$(PROGS)/cluster_mass_function.cc
	
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

# FIND CWWK SPECTRA IN LSST SPECTRA 
$(EXE)/findCWWinLSSTspectra : $(OBJ)/findCWWinLSSTspectra.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/findCWWinLSSTspectra $(OBJ)/findCWWinLSSTspectra.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/findCWWinLSSTspectra.o : $(PROGS)/findCWWinLSSTspectra.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/findCWWinLSSTspectra.o $(PROGS)/findCWWinLSSTspectra.cc
	
# FIT LSST SPECTRA TO CWWK
$(EXE)/fitLSSTspectra : $(OBJ)/fitLSSTspectra.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitLSSTspectra $(OBJ)/fitLSSTspectra.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitLSSTspectra.o : $(PROGS)/fitLSSTspectra.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitLSSTspectra.o $(PROGS)/fitLSSTspectra.cc
	
# FIT SEDS TO COLORS
$(EXE)/fitSEDsToColors : $(OBJ)/fitSEDsToColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitSEDsToColors $(OBJ)/fitSEDsToColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitSEDsToColors.o : $(PROGS)/fitSEDsToColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitSEDsToColors.o $(PROGS)/fitSEDsToColors.cc

# GENERATE A NEW TEMPLATE SET
$(EXE)/generateTemplSet : $(OBJ)/generateTemplSet.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/generateTemplSet $(OBJ)/generateTemplSet.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/generateTemplSet.o : $(PROGS)/generateTemplSet.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/generateTemplSet.o $(PROGS)/generateTemplSet.cc
	
# GET SED COLORS
$(EXE)/getSEDcolors : $(OBJ)/getSEDcolors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/getSEDcolors $(OBJ)/getSEDcolors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/getSEDcolors.o : $(PROGS)/getSEDcolors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/getSEDcolors.o $(PROGS)/getSEDcolors.cc
	
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

# CALCULATE LSST COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/lsstPicklesLibrary : $(OBJ)/lsstPicklesLibrary.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lsstPicklesLibrary $(OBJ)/lsstPicklesLibrary.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lsstPicklesLibrary.o : $(PROGS)/lsstPicklesLibrary.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lsstPicklesLibrary.o $(PROGS)/lsstPicklesLibrary.cc
	
# LYMAN ALPHA ALONG LINE OF SIGHT CONVERTED TO DENSITY
$(EXE)/lymanAlphaToDensity : $(OBJ)/lymanAlphaToDensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lymanAlphaToDensity $(OBJ)/lymanAlphaToDensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lymanAlphaToDensity.o : $(PROGS)/lymanAlphaToDensity.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lymanAlphaToDensity.o $(PROGS)/lymanAlphaToDensity.cc
	
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
	
# PRODUCE BPZ CATALOG
$(EXE)/produceCatalogForBPZ : $(OBJ)/produceCatalogForBPZ.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/produceCatalogForBPZ $(OBJ)/produceCatalogForBPZ.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/produceCatalogForBPZ.o : $(PROGS)/produceCatalogForBPZ.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/produceCatalogForBPZ.o \
	$(PROGS)/produceCatalogForBPZ.cc 
	
# PROJECT TEMPLATES
$(EXE)/projectTemplates : $(OBJ)/projectTemplates.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/projectTemplates $(OBJ)/projectTemplates.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/projectTemplates.o : $(PROGS)/projectTemplates.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/projectTemplates.o $(PROGS)/projectTemplates.cc 
	
# QUICK SIM MAGS
$(EXE)/quickSimMags : $(OBJ)/quickSimMags.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/quickSimMags $(OBJ)/quickSimMags.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/quickSimMags.o : $(PROGS)/quickSimMags.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/quickSimMags.o $(PROGS)/quickSimMags.cc 
	
# RESTFRAME COLORS
$(EXE)/restFrameColors : $(OBJ)/restFrameColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/restFrameColors $(OBJ)/restFrameColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/restFrameColors.o : $(PROGS)/restFrameColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/restFrameColors.o $(PROGS)/restFrameColors.cc 

# SIMULATE CATALOG OF BASIC GALAXY PROPERTIES FROM OVER-DENSITY GRID
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

# CALCULATE SDSS COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/sdssPicklesLibrary : $(OBJ)/sdssPicklesLibrary.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sdssPicklesLibrary $(OBJ)/sdssPicklesLibrary.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssPicklesLibrary.o : $(PROGS)/sdssPicklesLibrary.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sdssPicklesLibrary.o $(PROGS)/sdssPicklesLibrary.cc

# SIMULATE GALAXY DISTRIBUTION WITH CLUSTERS
$(EXE)/sim_with_clusters : $(OBJ)/sim_with_clusters.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sim_with_clusters $(OBJ)/sim_with_clusters.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sim_with_clusters.o : $(PROGS)/sim_with_clusters.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sim_with_clusters.o $(PROGS)/sim_with_clusters.cc
	
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
	
# SIMULATE CFHT OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateCFHTobs : $(OBJ)/simulateCFHTobs.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateCFHTobs $(OBJ)/simulateCFHTobs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateCFHTobs.o : $(PROGS)/simulateCFHTobs.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateCFHTobs.o \
	$(PROGS)/simulateCFHTobs.cc 
	
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

# Calculate magnitudes for IGM lines of sight
$(EXE)/lsstSbIgmColors : $(OBJ)/lsstSbIgmColors.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/lsstSbIgmColors $(OBJ)/lsstSbIgmColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)  
	
$(OBJ)/lsstSbIgmColors.o : $(PROGS)/lsstSbIgmColors.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lsstSbIgmColors.o $(PROGS)/lsstSbIgmColors.cc 


# ADD GAUSSIAN Z ERROR TO CATALOG
$(EXE)/addGausszerr : $(OBJ)/addGausszerr.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/addGausszerr $(OBJ)/addGausszerr.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/addGausszerr.o : $(BAOPROGS)/addGausszerr.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/addGausszerr.o \
	$(BAOPROGS)/addGausszerr.cc 
	
# COMPUTE POWER SPECTRA
$(EXE)/computepsfromarray : $(OBJ)/computepsfromarray.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/computepsfromarray $(OBJ)/computepsfromarray.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/computepsfromarray.o : $(BAOPROGS)/computepsfromarray.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/computepsfromarray.o \
	$(BAOPROGS)/computepsfromarray.cc 
	
# FIT K BAO
$(EXE)/fitkbao : $(OBJ)/fitkbao.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitkbao $(OBJ)/fitkbao.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitkbao.o : $(BAOPROGS)/fitkbao.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitkbao.o \
	$(BAOPROGS)/fitkbao.cc
	
# CALCULATE PHOTO-Z CONVOLUTION FUNCTION
$(EXE)/getpzconvf : $(OBJ)/getpzconvf.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/getpzconvf $(OBJ)/getpzconvf.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/getpzconvf.o : $(BAOPROGS)/getpzconvf.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/getpzconvf.o \
	$(BAOPROGS)/getpzconvf.cc 
	
# CALCULATE SELECTION FUNCTION OF OBSERVED CATALOG
$(EXE)/getsf : $(OBJ)/getsf.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/getsf $(OBJ)/getsf.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/getsf.o : $(BAOPROGS)/getsf.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/getsf.o \
	$(BAOPROGS)/getsf.cc 
	
# GRID GALAXY DATA
$(EXE)/grid_data : $(OBJ)/grid_data.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/grid_data $(OBJ)/grid_data.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/grid_data.o : $(BAOPROGS)/grid_data.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/grid_data.o \
	$(BAOPROGS)/grid_data.cc
	
# SIMULATE RANDOM CATALOGS
$(EXE)/sim_mcgrids : $(OBJ)/sim_mcgrids.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sim_mcgrids $(OBJ)/sim_mcgrids.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sim_mcgrids.o : $(BAOPROGS)/sim_mcgrids.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sim_mcgrids.o \
	$(BAOPROGS)/sim_mcgrids.cc
	
# GET DATA SUB GRID
$(EXE)/subfromfull : $(OBJ)/subfromfull.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/subfromfull $(OBJ)/subfromfull.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/subfromfull.o : $(BAOPROGS)/subfromfull.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/subfromfull.o \
	$(BAOPROGS)/subfromfull.cc


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

# TEST LF 
$(EXE)/testLF : $(OBJ)/testLF.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testLF $(OBJ)/testLF.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLF.o : $(PROGS)/testLF.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/testLF.o $(PROGS)/testLF.cc 

# TEST Lyman-alpha absorption calculation parts 
$(EXE)/testLymanAlphaAbs : $(OBJ)/testLymanAlphaAbs.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testLymanAlphaAbs $(OBJ)/testLymanAlphaAbs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLymanAlphaAbs.o : $(PROGS)/testLymanAlphaAbs.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/testLymanAlphaAbs.o $(PROGS)/testLymanAlphaAbs.cc 
	
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
	$(CXXLINK) -o $(EXE)/testSimReadKcorr $(OBJ)/testSimReadKcorr.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testSimReadKcorr.o : $(PROGS)/testSimReadKcorr.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testSimReadKcorr.o $(PROGS)/testSimReadKcorr.cc 

# TEST SIMULATE IGM
$(EXE)/testsimulateIGM : $(OBJ)/testsimulateIGM.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testsimulateIGM $(OBJ)/testsimulateIGM.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testsimulateIGM.o : $(PROGS)/testsimulateIGM.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testsimulateIGM.o $(PROGS)/testsimulateIGM.cc 
	
# TEST SIMULATION
$(EXE)/testSimulation : $(OBJ)/testSimulation.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/testSimulation $(OBJ)/testSimulation.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

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
	

# TIME TEST OF GETMAG
$(EXE)/time_test_getmag : $(OBJ)/time_test_getmag.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/time_test_getmag $(OBJ)/time_test_getmag.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/time_test_getmag.o : $(PROGS)/time_test_getmag.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/time_test_getmag.o $(PROGS)/time_test_getmag.cc
	
	
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
	mkdir -p $(OBJ)/usr/include/c++/4.6/bits/stl_vector.h:142:9: error: invalid use of incomplete type ‘struct SOPHYA::FunRan’

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



