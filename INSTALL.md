Installation Instructions
=========================

    
    * Install the external libraries 
             - CFITSIO 
             - FFTW
             - LessTif (if you want to use spiapp; see below)

    Perform the usual ./configure, make, make install. 
    Type ./configure -h to check the configure options.  Make sure --prefix=[INSTALLDIR] 
    if you don't want to use the default install directory.
    If you want to use double precision with FFTW (which you will!) you need to do this twice,
    a second time adding the --enable-float argument to ./configure


    * Install SOPHYA

        * Configure Sophya 
        
                 - cd BuildMgr/  
                 - Create the target installation directory (if not already existing):  mkdir $WHEREVER/SObjs 
                 - Define the SOPHYABASE environment variable:  setenv $SOPHYABASE $WHEREVER/SObjs/ 
                 - Run the configure script: ./configure -sbase $SOPHYABASE -scxx g++ -extp /usr/local -noext lapack -noext astro 
                   (or include lapack by removing -noext lapack [recommended])

        * Compile and build the libraries 

                 - make libs extlibs slb slbext  
                 - Or if you want to also compile spiapp (need LessTif):  make all slball 

        * If compiling spiapp ... 

                 - You should change the path for Motif/Lesstif include files and libraries in the file 
                   $SOPHYABASE/include/sophyamake.inc. 
                 - Locate the two lines: PIINC = -I/usr/X11R6/include/ -I/sw/include and PILIBS = -L/sw/lib -lXm -L/usr/X11R6/lib  -lXt -lX11
                 - Change these to: PIINC = -I/usr/X11R6/include/ -I/usr/local/Sophya/ExtLibs/include and 
                   PILIBS = -L/usr/local/Sophya/ExtLibs/lib -lXm -L/usr/X11R6/lib -lXt -lX11 
                 - spiapp should then compile.

        * Compile some test and utility programs 

                 - make basetests prgutil 
                 - To compile spiapp: make piapp

        * Clean .o files if not going to recompile after: make cleanobj

    
    * SOPHYA set up 

        - Set the Sophya environment variable, the repository where the libraries and the programs stand, 
           for example:  setenv SOPHYABASE /home/.../SObjs/
        - Using Sophya shared libraries: as you would do for any software, add the repository name of
          Sophya shared libraries to LD_LIBRARY_PATH: setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${SOPHYABASE}/slb
        - Using Sophya programs: as you would do for any software, add the repository name of 
          Sophya programs to PATH: setenv PATH ${PATH}:${SOPHYABASE}/exe
        - Making your own programs: if you are writing code using Sophya you may like to put
          "include $(SOPHYABASE)/include/sophyamake.inc" in your own Makefile. That will provide you with 
          default variables for the libraries repositories and names as well as standard compilation and 
          link options for various current OS and compilers. (please see:  more $SOPHYABASE/include/sophyamake.inc)
    
    * Download the git repo 

        - [instructions to get repo] 
        - cd DirectSim 
        - type make


