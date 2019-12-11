# rhalphalib

![Ralph](https://upload.wikimedia.org/wikipedia/en/thumb/1/14/Ralph_Wiggum.png/220px-Ralph_Wiggum.png)

Requires:
  - Python 3 (to change...)
  - ROOT (& pyROOT)
  - coffea

## Requirements
Standalone model creation requires:
  - Python 2.7+ or 3.6+
  - `numpy >= 1.14`

RooFit+combine rendering requires:
  - `ROOT < 6.18` (i.e. LCG96 is too recent, CMSSW 8 combine cannot handle it.  LCG95a is fine)

Use in combine requires, well, [combine](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit).
The CMSSW 10 (CC7) [recipe](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#cc7-release-cmssw_10_2_x-recommended-version)
satisfies the requirements, however the CMSSW 8 recipe has a too old version of numpy.

There is a python 3 compatible standalone fork of combine [available](https://github.com/guitargeek/combine).
It is also possible to render the model folder using the quickstart recipe, and then move the folder or switch
environments to a CMSSW+combine environment and proceed from there.

## Setup
  * Rhalphalib environment:
  ```bash
  mkdir DAZSLE/rhalphabet
  cd DAZSLE/rhalphabet
  #git clone git@github.com:DryRun/coffeandbacon.git # For environment setup script
  wget https://raw.githubusercontent.com/DryRun/coffeandbacon/master/setup_lcg.sh  
  wget https://raw.githubusercontent.com/DryRun/coffeandbacon/master/env_lcg.sh
  source setup_lcg.sh
  (source env_lcg.sh)
  git clone git@github.com:stalbrec/rhalphalib.git
  cd rhalphalib
  python test_rhalphabet.py
  ```
  * Combine environment (new shell):
  (following instructions from [Combined Wiki](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) )
  ```bash
  export SCRAM_ARCH=slc6_amd64_gcc530
  cmsrel CMSSW_8_1_0
  cd CMSSW_8_1_0/src
  git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
  cd HiggsAnalysis/CombinedLimit
  cmsenv
  git fetch origin
  git checkout v7.0.13
  scramv1 b clean; scramv1 b # always make a clean build

  cmsenv
  cd ../../../../rhalphalib/testModel
  source build.sh
  combine -M FitDiagnostics testModel_combined.txt --plots
  ```

## UHH producer
The producer script reads the model from a .json file, builds data cards and runs combine.
Also paths to the combine installation, root files and a 2D grid are taken from the .json.
The grid is created before running the analysis code containing the bins in pt and eta that are used to vary PF particles.
There is one grid per particle category (charged, neutral,...). The grid is constructed as root file and contains all necessary information to construct nuisance parameters that are consistent with the varied histograms. If QCD is part of the given samples, rhalphabet takes care of the the background estimation in the pass region.

How to run:
```bash
python uhh_producer.py Model.json
```
