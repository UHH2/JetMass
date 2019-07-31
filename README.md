# rhalphalib

![Ralph](https://upload.wikimedia.org/wikipedia/en/thumb/1/14/Ralph_Wiggum.png/220px-Ralph_Wiggum.png)

Requires:
  - Python 3 (to change...)
  - ROOT (& pyROOT)
  - coffea

  * Rhalphalib environment:
  ```bash
  mkdir DAZSLE/rhalphabet
  cd DAZSLE/rhalphabet
  git clone git@github.com:DryRun/coffeandbacon.git # For environment setup script
  source coffeandbacon/setup_lcg.sh
  git clone git@github.com:DryRun/rhalphalib.git
  (source coffeandbacon/env_lcg.sh)
  cd rhalphalib
  python test_rhalphabet.py
  ```
  * Combine environment (new shell):
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
