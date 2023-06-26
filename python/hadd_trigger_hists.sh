#!/bin/bash

YEARS=(UL16preVFP UL16postVFP UL17 UL18)

export WORKDIR=${1:-/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/ttbarTrees/ForTriggerEff}

function hadd_sample {
  YEAR=${1:-UL17}
  SAMPLE=${2:-Data}
  hadd -T -f ${WORKDIR}/${SAMPLE}_${YEAR}.root ${WORKDIR}//uhh2*${SAMPLE}*${YEAR}*.root
  
}

export -f hadd_sample

SAMPLES=(Data QCD)
for SAMPLE in ${SAMPLES[@]}
do
  printf "%s\n" ${YEARS[@]} | xargs -I{} bash -c 'hadd_sample "$@"' _ {} $SAMPLE
done
