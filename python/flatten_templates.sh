#!/bin/bash

if [ ! -d $outdir ]; then
  mkdir -p  $outdir
fi


function flatten_eta_regions {
  for YEAR in ${YEARS[@]};
  do
    for eta_region in barrel endcap;
    do
      # JEC on pt
      ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${eta_region} --JEC "pt" --eta ${eta_region}
      ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_jec_up_${eta_region} --JEC "pt" --eta ${eta_region}
      ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_jecpt_jec_down_${eta_region} --JEC "pt" --eta ${eta_region}

      # JEC on pt&mJ
      ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_${eta_region} --JEC "pt&mJ" --eta ${eta_region}
      ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_jec_up_${eta_region} --JEC "pt&mJ" --eta ${eta_region}
      ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_jec_down_${eta_region} --JEC "pt&mJ" --eta ${eta_region}
    done
  done
}

flatten_templates () {
  YEARS=(UL16preVFP UL16postVFP UL17 UL18)
  outdir=flat_templates/
  indir=coffea_hists/

  VAR=${1:-none}
  MASS=$2
  if [ "$MASS" == "mPnet" ];then
  outdir="${outdir}/mPnet/"
  fi

  NAMEPREFIX=${3:-""}

  if [ "$VAR" == "none" ]; then
    echo "You have to provide a variation name!"
    exit -1
  fi
  
  for YEAR in ${YEARS[@]};
  do
    NAME=templates_${YEAR}${NAMEPREFIX}
    if [ "$VAR" == "nominal" ]; then
      echo "flattening ${YEAR} templates for variation of ${VAR}"
      # JEC on pt
      ./create_root_templates.py -i $indir/${NAME}.coffea -o $outdir/${NAME}_1d_jecpt --JEC "pt" --mass ${MASS}
      ./create_root_templates.py -i $indir/${NAME}.coffea -o $outdir/${NAME}_1d_jecpt_unfolding --JEC "pt" --unfolding --mass ${MASS}
      
      # JEC on pt&mJ
      ./create_root_templates.py -i $indir/${NAME}.coffea -o $outdir/${NAME}_1d --JEC "pt&mJ" --mass ${MASS}
      ./create_root_templates.py -i $indir/${NAME}.coffea -o $outdir/${NAME}_1d_unfolding --JEC "pt&mJ" --unfolding --mass ${MASS}
    elif [ "$VAR" == "toppt_off" ]; then
      echo "flattening ${YEAR} templates for variation of ${VAR}"
      # JEC on pt
      ./create_root_templates.py -i $indir/${NAME}_${VAR}.coffea -o $outdir/${NAME}_1d_jecpt_${VAR} --JEC "pt" --mass ${MASS}
      ./create_root_templates.py -i $indir/${NAME}_${VAR}.coffea -o $outdir/${NAME}_1d_jecpt_unfolding_${VAR} --JEC "pt" --unfolding --mass ${MASS}

      # JEC on pt&mJ
      ./create_root_templates.py -i $indir/${NAME}_${VAR}.coffea -o $outdir/${NAME}_1d_${VAR} --JEC "pt&mJ" --mass ${MASS}
      ./create_root_templates.py -i $indir/${NAME}_${VAR}.coffea -o $outdir/${NAME}_1d_unfolding_${VAR} --JEC "pt&mJ" --unfolding --mass ${MASS}
    else
      for DIR in up down;
      do
        echo "flattening ${YEAR} templates for variation of ${VAR}"
        # JEC on pt
        ./create_root_templates.py -i $indir/${NAME}_${VAR}_${DIR}.coffea -o $outdir/${NAME}_1d_jecpt_${VAR}_${DIR} --JEC "pt" --mass ${MASS}
        ./create_root_templates.py -i $indir/${NAME}_${VAR}_${DIR}.coffea -o $outdir/${NAME}_1d_jecpt_unfolding_${VAR}_${DIR} --JEC "pt" --unfolding --mass ${MASS}

        # JEC on pt&mJ
        ./create_root_templates.py -i $indir/${NAME}_${VAR}_${DIR}.coffea -o $outdir/${NAME}_1d_${VAR}_${DIR} --JEC "pt&mJ" --mass ${MASS}
        ./create_root_templates.py -i $indir/${NAME}_${VAR}_${DIR}.coffea -o $outdir/${NAME}_1d_unfolding_${VAR}_${DIR} --JEC "pt&mJ" --unfolding --mass ${MASS}
      done
    fi
  done
}

# run the things
VARS_ALL=(nominal jec fsr isr triggersf pu toppt_off)
# VARS=(${2:-${VARS_ALL[@]}})
VARS=(${VARS_ALL[@]})
MAXPROCS=5
# NAMEPREFIX=_particlenet
NAMEPREFIX=${1:-}
# MASS=mPnet
MASS=mjet
export -f flatten_templates
printf '%s\n' "${VARS[@]}" | xargs -n1 -P${MAXPROCS} -I{} bash -c 'flatten_templates "$@"' _ {} $MASS $NAMEPREFIX

