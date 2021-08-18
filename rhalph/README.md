# Jet mass fit

## Workspace creation and running the fit

Both the workspace creation and fit should be done in the CMSSW environment that is installed with combine (using the `setup_rhalph.sh` script)
So to setup the suitable environment simply do a `cmsenv` inside the CMSSW directory and you should be good to go.


The workspace is created using `python jetmass.py CONFIG.json`, where the `CONFIG.json` file specifies the region setup and where the histograms are located.
Examples for "single sample" configs are `WJets.json`,`ZbbJets.json` and `TTbar.json`, while `CombinedFit.json` combines all three of those to one simultaneous fit.
 
WARNING: With the current combine version depending on a CMSSW-version, where ROOT is not available in python3 you should not try to run `jetmass.py` in python3