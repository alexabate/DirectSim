/** \mainpage DirectSim (aka Franzona)
*
* @authors Alexandra Abate (UArizona), Reza Ansari (LAL), Lidens Cheng (UArizona), 
* Matthew Kirby (UArizona), Christophe Magneville (Saclay)
*
* \tableofcontents
*
* <hr>
* @section intro Introduction
* Simulation of galaxy catalogs for cosmological analysis. Please see the 
* documentation in the "Related Pages" tab above on how to run 
* the programs and hopefully get what you want!
*
* This code repository is hosted on <A HREF="https://github.com/alexabate/DirectSim"> github </A> 
*
* <hr>
* @section install Installation
*
* For installation instructions, please see the file <CODE>INSTALL.md</CODE> in the main 
* repository directory, or see the install page under the "Related Pages" tab.
*
* <hr>
* @section dirs Directory Structure
* <CODE> $DIRECTSIM </CODE> refers to the root directory of DirectSim
* 
* The source code:
* - <B>classes</B> Contains all the definitions and source code for the classes
* - <B>progs</B> Contains the programs
* - <B>baoprogs</B> Contains the programs related to power spectrum estimation and BAO
*
* Compiled code (directories created on first compilation):
* - <B>objs</B> Object files 
* - <B>exe</B> Executables
*
* Also note the existence of the following directories:
* - <B>filters</B>: Contains filter sets and filter transmissions.  Filter sets are in
*            files of the form: <CODE>[SETNAME-ALL-IN-CAPS].filters</CODE>.  Filter
*            transmissions are in files of the form: <CODE>[band]_[instrument]_[source].txt</CODE>.
*            If adding more filters sets please keep to this same convention.
* - <B>SEDs</B>: Contains SED libraries and spectra. SED libraries are in files of the 
*         form: <CODE>[LIBNAME-ALL-IN-CAPS].list</CODE>. Spectra are in other files
*         with extentions like <CODE>.txt</CODE>, <CODE>.sed</CODE> or <CODE>.spec</CODE>.
*         If adding more SED libraries please keep to this same convention.
* - <B>kCorrections</B>: If <CODE>calculateKcorrections</CODE> has been run, contains pre-calculated
*                 k-corrections for different SED-filter combinations.
* - <B>LFs</B>: Contains luminosity functions
* - <B>dat</B>: Contains data required for parts of the simulation
* - <B>docs</B>: Where the doxygen documation is generated
* - <B>picklesStellarLib</B>: Pickles stellar library, relating to star simulation
*
* Try to keep to something like LSST coding standards:
* <A HREF="https://dev.lsstcorp.org/trac/wiki/CodingStandardsDesirable"> Coding Standards </A>
*
* For info on how to use git see:<BR>
* <A HREF="http://www.slideshare.net/saharabeara/advanced-git-tutorial"> Advanced Git Tutorial </A> by Sarah Sharp <BR>
* The <A HREF="http://git-scm.com/book"> Pro Git book </A> <BR>
*
* <hr>
* @section notes Release notes
* Details about current release ...
* <hr>
* @section requirements Requirements
* @verbinclude requirements ...
* <hr> 
* @todo Click here for a list of updates to be made
*
*/
