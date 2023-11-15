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
  elif [ "$TAGGER" == "particlenetDDT" ]; then
    NAMESUFFIX="_particlenetDDT"
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

  VJETSONLY=${4}
  if [[ ${VAR} == *"qcd"* ]]; then
    VJETSONLY="--VJetsOnly"
  elif [[ ${VAR} == *"ewk"* ]]; then
    VJETSONLY="--VJetsOnly"
  fi
  if [ ${SCALEOUT} -eq "0" ]; then
    echo "running locally"
  else
    SCALEOUTARG="--scaleout ${SCALEOUT}"
  fi
  if [ ! -d ${OUTDIR} ]; then
    mkdir -p $OUTDIR
  fi
  echo $SCALEOUT $VAR
  if [ "$TAGGER" == "substructure" ]; then
    NAMESUFFIX=""
  elif [ "$TAGGER" == "particlenet" ]; then
    NAMESUFFIX="_particlenet"
  elif [ "$TAGGER" == "particlenetDDT" ]; then
    NAMESUFFIX="_particlenetDDT"
  else
    echo "You did not provide a valid tagger"
  fi
  
  for YEAR in UL16preVFP UL16postVFP UL17 UL18
  do
    if [ ${VAR} == "nominal" ]; then
      echo "default -> ${YEAR}" 
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}.coffea --year ${YEAR} --tagger ${TAGGER} --scaleout ${SCALEOUT} ${VJETSONLY}
    elif [[ ${VAR} == *"jec"* ]]; then
      echo "jec-variation (${VAR}) -> ${YEAR}" $VAR $YEAR
      DIRECTION=${VAR#*_}
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --JEC ${DIRECTION} --scaleout ${SCALEOUT} ${VJETSONLY}
    elif [[ ${VAR} == *"trigger"* ]]; then
      echo "triggersf-variation (${VAR}) -> ${YEAR}" $VAR $YEAR
      DIRECTION=${VAR#*_}
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --triggersf ${DIRECTION} --scaleout ${SCALEOUT} ${VJETSONLY}
    else
      echo "variation ${VAR} -> ${YEAR}"
      ./JetMassTemplateProcessor.py -o ${OUTDIR}/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --variation ${VAR} ${SCALEOUTARG} ${VJETSONLY}
    fi
  done
}

export TMPDIR=/tmp/
SCALEOUT=$1
VARIATION=${2:-nominal}
VARIATIONS_ALL=(nominal jec_up jec_down triggersf_up triggersf_down isr_up isr_down fsr_up fsr_down pu_up pu_down toppt_off)
# VARIATIONS_PART1=(nominal jec_up jec_down jer_up jer_down triggersf_up triggersf_down)
VARIATIONS_PART1=(nominal jec_up jec_down triggersf_up triggersf_down)
VARIATIONS_PART2=(isr_up isr_down fsr_up fsr_down pu_up pu_down toppt_off)
VARIATIONS_PART3=(v_qcd_up v_qcd_down w_ewk_up w_ewk_down z_ewk_up z_ewk_down)
VARIATIONS_GENHISTS=(nominal v_qcd_up v_qcd_down w_ewk_up w_ewk_down)
TAGGER=${3:-substructure}
# VJETSONLY="--VJetsOnly"
VJETSONLY=""
if [ "$VARIATION" == "all" ];
then
  VARIATIONS=${VARIATIONS_ALL[@]}
elif [ "$VARIATION" == "part1" ];
then
  VARIATIONS=${VARIATIONS_PART1[@]}
elif [ "$VARIATION" == "part2" ];
then
  VARIATIONS=${VARIATIONS_PART2[@]}
elif [ "$VARIATION" == "part3" ];
then
  VARIATIONS=${VARIATIONS_PART3[@]}
elif [ "$VARIATION" == "mgen" ];
then
  VARIATIONS=${VARIATIONS_GENHISTS[@]}
else
  VARIATIONS=($VARIATION)
fi
echo ${VARIATIONS[@]}
  
for VAR in ${VARIATIONS[@]};
do
    submit_templates $SCALEOUT $VAR $TAGGER $VJETSONLY
done
