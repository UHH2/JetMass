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
ls configs/nominal/*UL*.py | xargs -n1 -P10 -I% python jetmass.py -M jms % -workdir <workdir>
ls <workdir>/*/config.json | sed -r 's|/[^/]+$||'  | xargs echo "python3 extractMassScales.py --fits"
 ```

For the JMS comparison plots with different treatment of JEC Uncertainties we need 4 fits per year and per sample(-combination) (TTBar, VJets, Combined).
To do this all together use `JECStudyFits.py`:

```shell
#submit condor jobs for each fit (120 in total)
python JECStudyFits.py --submit

#extract fitResults and write to json files (4 in total)
python JECStudyFits.py --extract
```



### Impacts

- inital fit, and fits
```shell
combineTool.py -M Impacts -d model_combined.root --doInitialFit --robustFit 1 -m 0 --exclude rgx{qcdparam*} --cminDefaultMinimizerStrategy 0
# make sure to have enough disk space, some fit might give problems which will result in massive logs. Better move to separate workdir to execute this.
combineTool.py -M Impacts -d model_combined.root -m 0 --robustFit 1 --doFits --exclude rgx{qcdparam*} --parallel 100 --job-mode condor --cminDefaultMinimizerStrategy 0
```

- plot impacts for each POI and use renaming json to texify POI name

```shell
combineTool.py -M Impacts -d model_combined.root -m 0 -o impacts.json

for ipt in 0 1 2 3; do echo $ipt ;printf "%s\n" 0 1 2 3 |xargs -I{} -P 5 -n1 plotImpacts.py -i impacts.json -o impacts_r_ptgen${ipt}_msdgen{} --POI r_ptgen${ipt}_msdgen{} --per-page 35 --translate ../../POI_rename.json --height 600;done
```
