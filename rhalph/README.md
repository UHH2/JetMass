# Jet mass fit

## Installation
### setup CMSSW environment with `combine` installation

Follow the most recent [installation instructions](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#installation-instructions):

```shell
# create CMSSW release are and clone combine
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

# update to recommended tag and compile
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.0.0
scramv1 b clean; scramv1 b # always make a clean build
```

and install CombineHarvester for additional utilities:

```shell
# clone; sparse checkout
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/main/CombineTools/scripts/sparse-checkout-https.sh)

# build
scram b
```
#### setup environment
Once installed the environment can be activated by:

```shell
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_11_3_4/src
cmsenv
```

### Installing `python3` venv for additional utilities
make sure to setup the virtual environment as described in the main [README.md](../README.md)

## Workspace creation and running the fit

### create the workspace and run the fit

The workspace is created using `python jetmass.py CONFIG.<json,py>`, where the `CONFIG.<json,py>` file specifies the region setup and where the histograms are located.

```shell
usage: jetmass.py [-h] [--justplots] [--build] [--job_index JOB_INDEX]
                  [--minimalModel] [--skipTemplatePlots]
                  [--customCombineWrapper] [--noNuisances] [--JMRparameter]
                  [--prefitAsimov] [--pTdependetJMRParameter] [--noNormUnc]
                  [--skipExtArgRender] [--seed SEED]
                  [--combineOptions COMBINEOPTIONS] [--verbose VERBOSE] -M
                  {unfolding,jms,default}
                  [--uncertainty-breakdown UNCERTAINTY_BREAKDOWN [UNCERTAINTY_BREAKDOWN ...]]
                  config
```

the most important arguments:
```shell
positional arguments:
  config                path to json with config

optional arguments:
  --justplots           just make plots.
  --build               just build workspace for combine
  -M {unfolding,jms,default}, --mode {unfolding,jms,default}
                        choose what fit should be run. 'unfolding' performs
                        MaxLik. unfolding, jms measures JetMassScale
                        variations,'default' measures signal cross section.
```

 Per default this will both create and execute some fitting procedure. 

 ## JMS fits

 First perform some fits, e.g. all UL years:

 ```shell
ls configs/nominal/*UL*.py | xargs -n1 -P10 -I% python jetmass.py -M jms %
ls */config.json | sed -r 's|/[^/]+$||'  | xargs echo "python extractMassScales.py --fits"
 ```