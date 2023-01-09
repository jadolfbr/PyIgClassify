#!/usr/bin/env bash

#Run this to setup any auxilliary python packages.  Please use Sudo in order to allow installation.
#Requires installation of curl package.  This is shipped with newer mac OS.  Linux: sudo apt-get install curl
#Weekly updates require Linux machine.  To use pip, setuptools needs to be installed as well.  sudo apt-get install python-setuptools

echo Installing HMMER3.0
curl -O ftp://selab.janelia.org/pub/software/hmmer3/3.0/hmmer-3.0.tar.gz
tar xvf hmmer-3.0.tar.gz
cd hmmer-3.0
./configure
make
make install
make check
cd -

echo Installing Pip
curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
python get-pip.py

echo Installing Numpy
pip install numpy

echo Installing scipy
pip install scipy

echo Installing BioPython
pip install biopython

echo Installing weblogo
pip install weblogo


#echo Setting up R packages
#Rscript setup_R_packages.R

