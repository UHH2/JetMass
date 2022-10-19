[![flake8](https://github.com/UHH2/JetMass/actions/workflows/lint.yml/badge.svg?branch=master)](https://github.com/UHH2/JetMass/actions/workflows/lint.yml)
# JetMass
Create histograms for jet mass calibration
- make sure to get ROOT files with EW-Correction and QCD k-factors from [UHHNtupleConverter](https://github.com/Diboson3D/UHHNtupleConverter):

```
mkdir NLOWeights; cd NLOWeights
wget https://github.com/Diboson3D/UHHNtupleConverter/raw/master/NLOweights/WJetsCorr.root
wget https://github.com/Diboson3D/UHHNtupleConverter/raw/master/NLOweights/ZJetsCorr.root
```
