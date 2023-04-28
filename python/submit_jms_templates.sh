#!/bin/bash

function submit_templates {
  SCALEOUT=$1
  VAR=${2:-nominal}
  TAGGER=${3:-substructure}
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
      nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}${NAMESUFFIX}.coffea --year ${YEAR} --tagger ${TAGGER} --scaleout ${SCALEOUT} > ${YEAR}${NAMESUFFIX}.stdout &
    else
      echo "variation ${VAR} -> ${YEAR}"
      nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}${NAMESUFFIX}_${VAR}.coffea --year ${YEAR} --tagger ${TAGGER} --variation ${VAR} --scaleout ${SCALEOUT} > ${YEAR}_${VAR}.stdout &
    fi
  done
}

export TMPDIR=/tmp/
SCALEOUT=$1
VARIATION=${2:-nominal}
VARIATIONS=(nominal jec_up jec_down jer_up jer_down isr_up isr_down fsr_up fsr_down toppt_off)
TAGGER=${3:-substructure}
echo $VARIATION
if [ "$VARIATION" == "all" ];
then
  for VAR in ${VARIATIONS[@]};
  do
    submit_templates $SCALEOUT $VAR $TAGGER
  done
else
  submit_templates $SCALEOUT $VARIATION $TAGGER
fi
