DirectSim
=========

simulation of galaxy catalogs for cosmological analysis


Installation
------------

For full installation instructions, please see the file `INSTALL.md` in the main repository directory.

The "easiest" way to install DirectSim is to ask for an account on cosmo.physics.arizona.edu and then clone the DirectSim master branch into your home directory via either 1) or 2):


1) ```git clone git@github.com:alexabate/DirectSim.git``` You will need to generate ssh key according to these instructions: https://help.github.com/articles/generating-ssh-keys and create the ~/.ssh directory if it doesn't exist


2) ```git clone https://github.com/alexabate/DirectSim.git```


Then add the following lines to your ~/.bash_profile file (or whatever it's called)

```
# for SOPHYA
export SOPHYABASE=/software/lsstlib/Import/SObjs/
export PATH=${PATH}:${SOPHYABASE}/exe
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SOPHYABASE}/slb


# for ROOT
. /software/root/v5.32.00/bin/thisroot.sh

```

Reference documentation
-----------------------

Please see: http://www.u.arizona.edu/~abate/docs/html/index.html

Under the Related Pages tab there is documentation on how to use various programs and a description of what they do.

Warning! The webpages may not be completely up to date!


Repository directory structure
------------------------------

Here is a guide to the contents of the repo's various subdirectories:

* classes/ :      definitions and source code for the classes
* dat/ :    data files
* docs/ :      doxygen files for program documentations (see Related Pages tab)
* filters/ : filter sets and filter transmissions. Filter sets are in files of the form: `[SETNAME-ALL-IN-CAPS].filters`. Filter transmissions are in files of the form: `[band]_[instrument]_[source].txt`. If adding more filters sets please keep to this same convention
* LFs/ :   luminosity function data
* progs/ : program source codes
* SEDs/ :  SED libraries and spectra. SED libraries are in files of the form: `[LIBNAME-ALL-IN-CAPS].list`. Spectra are in other files with extentions like .txt, .sed or .spec. If adding more SED libraries please keep to this same convention.
