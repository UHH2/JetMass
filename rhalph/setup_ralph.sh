#!/bin/bash

if [ ! -d "rhalphalib" ]; then
  echo 'rhalphalib does not exists'
	git clone https://github.com/nsmith-/rhalphalib
fi

if uname -r | grep -q el6; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-slc6-gcc8-opt/setup.sh
else
  source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
fi

export PYTHONPATH=~/.local/lib/python3.6/site-packages:$PYTHONPATH
