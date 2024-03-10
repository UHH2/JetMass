# JetMass
[![flake8](https://github.com/UHH2/JetMass/actions/workflows/lint.yml/badge.svg)](https://github.com/UHH2/JetMass/actions/workflows/lint.yml)

This repository is home of the code for two W & top jet mass related analyses: the measurement of jet mass scale (JMS) data-to-simulation scale factors in the hadronic $W(q#bar{q}$)+jets and semileptonic channel ($t #bar{t} #rightarrow #mu #nu$+ jets  and the 2D likelihood-based unfolding of the jet mass in the hadronic channel.

Both analyses run off of the custom UHH2 MiniAOD N-tuples in the [UHH2 framework](https://github.com/UHH2/UHH2).
The first step is to run the pre-selection in the UHH2 framework. The config-files for all UL-years and the two samples can be found in the `config/` sub-directory. The cpp-modules can be found in the `include/` and `src/` sub-directories.


## W/Z+jets EW-Correction and QCD k-factors
- make sure to get ROOT files from 
[UHHNtupleConverter](https://github.com/Diboson3D/UHHNtupleConverter):

```
mkdir NLOWeights; cd NLOWeights
wget https://github.com/Diboson3D/UHHNtupleConverter/raw/master/NLOweights/WJetsCorr.root
wget https://github.com/Diboson3D/UHHNtupleConverter/raw/master/NLOweights/ZJetsCorr.root
```


## setting up python `venv`
- TODO
