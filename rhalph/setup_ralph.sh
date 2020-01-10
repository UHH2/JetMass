#!/bin/bash

if [ ! -d "rhalphalib" ]; then
  echo 'rhalphalib does not exists'
  git clone https://github.com/nsmith-/rhalphalib
fi

if [ ! -d "CMSSW_10_2_13" ]; then
  echo 'There is no CMSSW_10_2_13 installation in this directory'
  echo 'setting up CMSSW_10_2_13 and combine'
  export SCRAM_ARCH=slc6_amd64_gcc700
  cmsrel CMSSW_10_2_13
  cd CMSSW_10_2_13/src

  cmsenv

  git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
  cd HiggsAnalysis/CombinedLimit
  git fetch origin
  git checkout v8.0.1
  # scramv1 b clean; scramv1 b
  
  echo 'installing combineHarvester'
  git clone https://github.com/cms-analysis/CombineHarvester.git 

  scram b clean; scram b -j8
  
  cd ../../../../
  echo "installed CMSSW_10_2_13 with HiggsAnalysis/CombinedLimit"
  echo "to actually use rhalphalib & co. you should use a fresh shell"
fi

if uname -r | grep -q el6; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-slc6-gcc8-opt/setup.sh
else
  source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
fi

export PYTHONPATH=~/.local/lib/python3.6/site-packages:$PYTHONPATH
