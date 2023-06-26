#!/bin/bash

function submit_templates_parallel {
  SCALEOUT=$1
  VAR=${2:-nominal}
  TAGGER=${3:-substructure}
  OUTDIR="coffea_hists/"
  if [ ! -d ${OUTDIR} ]; then
    mkdir -p $OUTDIR
  fi
  echo $SCALEOUT $VAR
  if [ "$TAGGER" == "substructure" ]; then
    NAMESUFFIX=""
  elif [ "$TAGGER" == "particlenet" ]; then
    NAMESUFFIX="_particlenet"
  else
    echo "You did not provide a valid tagger"
  fi
  
  for YEAR in UL16preVFP UL16postVFP UL17 UL18
  do
    if [ ${VAR} == "nominal" ]; then
      echo "default -> ${YEAR}" 
      nohup ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}.coffea --year ${YEAR} --tagger ${TAGGER} --scaleout ${SCALEOUT} > ${YEAR}${NAMESUFFIX}.stdout &
    elif [[ ${VAR} == *"jec"* ]]; then
      echo "jec-variation (${VAR}) -> ${YEAR}" $VAR $YEAR
      DIRECTION=${VAR#*_}
      nohup ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --JEC ${DIRECTION} --scaleout ${SCALEOUT} > ${YEAR}${NAMESUFFIX}_${VAR}.stdout &
    else
      echo "variation ${VAR} -> ${YEAR}"
      nohup ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --variation ${VAR} --scaleout ${SCALEOUT} > ${YEAR}${NAMESUFFIX}_${VAR}.stdout &
    fi
  done
}

function submit_templates {
  SCALEOUT=$1
  VAR=${2:-nominal}
  TAGGER=${3:-substructure}
  OUTDIR="coffea_hists/"
  if [ ! -d ${OUTDIR} ]; then
    mkdir -p $OUTDIR
  fi
  echo $SCALEOUT $VAR
  if [ "$TAGGER" == "substructure" ]; then
    NAMESUFFIX=""
  elif [ "$TAGGER" == "particlenet" ]; then
    NAMESUFFIX="_particlenet"
  else
    echo "You did not provide a valid tagger"
  fi
  
  for YEAR in UL16preVFP UL16postVFP UL17 UL18
  do
    if [ ${VAR} == "nominal" ]; then
      echo "default -> ${YEAR}" 
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}.coffea --year ${YEAR} --tagger ${TAGGER} --scaleout ${SCALEOUT}
    elif [[ ${VAR} == *"jec"* ]]; then
      echo "jec-variation (${VAR}) -> ${YEAR}" $VAR $YEAR
      DIRECTION=${VAR#*_}
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --JEC ${DIRECTION} --scaleout ${SCALEOUT}
    else
      echo "variation ${VAR} -> ${YEAR}"
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --variation ${VAR} --scaleout ${SCALEOUT}
    fi
  done
}

export TMPDIR=/tmp/
SCALEOUT=$1
VARIATION=${2:-nominal}
VARIATIONS_ALL=(nominal jec_up jec_down jer_up jer_down isr_up isr_down fsr_up fsr_down toppt_off)
VARIATIONS_PART1=(nominal jec_up jec_down jer_up jer_down)
VARIATIONS_PART2=(isr_up isr_down fsr_up fsr_down toppt_off)
TAGGER=${3:-substructure}
if [ "$VARIATION" == "all" ];
then
  VARIATIONS=${VARIATIONS_ALL[@]}
elif [ "$VARIATION" == "part1" ];
then
  VARIATIONS=${VARIATIONS_PART1[@]}
elif [ "$VARIATION" == "part2" ];
then
  VARIATIONS=${VARIATIONS_PART2[@]}
else
  VARIATIONS=($VARIATION)
fi
echo ${VARIATIONS[@]}
  
for VAR in ${VARIATIONS[@]};
do
    submit_templates $SCALEOUT $VAR $TAGGER
done
