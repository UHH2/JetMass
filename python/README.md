
Some of these scripts are accessed using the `__init__.py` generated by SCRAM. 
If you run into problems with scripts that import a package/script from this directory, you can try to generate the init by doing a simple `scram b`.

# Workflow
1. first create flat trees that are skimmed using some pre-selection. See [JetMass/README.md](../README.md) for details.
2. fill 3D maps that are needed for DDT maps. You can use the coffea processor
   1. create 2D percentile maps for each DDT tagger for each year using `build_ddt_map.py`
3. fill templates using the coffea processor `JetMassTemplateProcessor.py`
   1. export suitable temporary directory path into `TMPDIR`, .e.g. `export TMPDIR=/tmp` (the coffea processor will be executed here to avoid bottlenecks such as slow afs/nfs.)
   1. adjust path to ddt-maps if necessary (e.g. see `self._n2ddtmaps`)
   2. adjust other pre-compiled data (i.e. trigger-sf, k-factors) if needed.
   3. right now the histograms have to be filled for each variation (`nominal jec_up jec_down triggersf_up triggersf_down isr_up isr_down fsr_up fsr_down pu_up pu_down toppt_off`) separetely and once for each of the tagger approaches (Substructure & ParticleNet)
      1. the script `submit_jms_templates.sh <number-of-jobs> <variations> <tagger>` will run variations sequentially for all years and a given tagging approach. So for all variations one can do `./submit_jms_templates.sh 50 all substructure` and `./submit_jms_templates.sh 50 all particlenet`. This will spawn a maximum of 100 dask-jobs on condor at any time, so will put minimal stress on the system. Be careful not to submit too many jobs that access the same files simultaneously (i.e. calling the script for each variation individually for one tagger). This can/will put stress on the filesystem the UHH2 flat tree ROOT files are stored on (e.g. DUST@DESY).
4. flatten output from coffea processor into 1d ROOT (and coffea) hists.
   1. `./flatten_templates.sh` and `./flatten_templates.sh _particlenet`
   2. create control plots
      1. `flatten_control_plots.sh` to create ROOT files with hists from coffea processor
      2. `python hadd_control_plots --preselection --output UHH_controlplots_presel` (before this make sure you hadded the output files from UHH2)
      3. finally `./data_mc_plots.sh` to make plots
5. templates are ready to be used by ralph (see [JetMass/rhalp/README.md](../rhalph/README.md) for details)
  
For the unfolding some intermediate steps are necessary once the UHH2 trees are ready and before the templates are created:

1. convert WJets trees to parquet for faster studies: `./tiny_tree.py --year <UL16preVFP,UL16postVFP,UL17,UL18> --tagger <n2ddt, pNetddt>`
   1. for parallel use something like: `printf "%s\n" UL16preVFP UL16postVFP UL17 UL18 | xargs -I{} -n1 -P4 ./tiny_tree.py --year {} --tagger n2ddt`
2. using the trees derive reco to gen level correction factors for soft drop mass: `./JMS_from_MC.py --tagger <n2ddt, pNetddt>`
   1. This will create one correctionset json including all years for one given tagger. Update the filename accordingly in `utils.py` - the filename includes part of the `sha512` sum to directly spot versions of sets.
3. using the corrections check if binning seems reasonable with `printf "%s\n" UL16preVFP UL16postVFP UL17 UL18 | xargs -I{} -n1 -P4 ./unfolding_binning.py --input WJetsToQQ_tinyTree_{}_<n2ddt,pNetddt>.parquet --tagger <n2ddt,pNetddt>`
   1. update binning in `JetMassTemplateProcessor.py` if necessary.
