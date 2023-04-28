#!/bin/bash

function submit_templates {
  SCALEOUT=$1
  VAR=${2:-nominal}
  echo $SCALEOUT $VAR
  for YEAR in UL16preVFP UL16postVFP UL17 UL18
  do
    if [ ${VAR} == "nominal" ]; then
      echo "default -> ${YEAR}" 
      nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}.coffea --year ${YEAR} --scaleout ${SCALEOUT} > ${YEAR}.stdout &
      # nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}_jec_up.coffea --year ${YEAR} --JEC up --scaleout ${SCALEOUT} > ${YEAR}_jec_up.stdout &
      # nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}_jec_down.coffea --year ${YEAR} --JEC down --scaleout ${SCALEOUT} > ${YEAR}_jec_down.stdout &
    else
      echo "variation ${VAR} -> ${YEAR}"
      nohup ./JetMassTemplateProcessor.py -o coffea_hists/templates_${YEAR}_${VAR}.coffea --year ${YEAR} --variation ${VAR} --scaleout ${SCALEOUT} > ${YEAR}_${VAR}.stdout &
    fi
  done
}

export TMPDIR=/tmp/
SCALEOUT=$1
VARIATION=${2:-nominal}
VARIATIONS=(nominal jec_up jec_down jer_up jer_down isr_up isr_down fsr_up fsr_down toppt_off)
echo $VARIATION
if [ "$VARIATION" == "all" ];
then
  for VAR in ${VARIATIONS[@]};
  do
    submit_templates $SCALEOUT $VAR
  done
else
  submit_templates $SCALEOUT $VARIATION
fi
