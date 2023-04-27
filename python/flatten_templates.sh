#!/bin/bash

YEARS=(UL16preVFP UL16postVFP UL17 UL18)
outdir=flat_templates/
indir=coffea_hists/

if [ ! -d $outdir ];then
mkdir -p  $outdir
fi

function flatten_nominal {
    for YEAR in ${YEARS[@]};
    do
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d --JEC "pt&mJ"
        # python plot_templates.py --input $outdir/templates_${YEAR}_1d.root --year ${YEAR}  
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_jec_up --JEC "pt&mJ"
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_jec_down --JEC "pt&mJ"

        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt --JEC "pt"
        # python plot_templates.py --input $outdir/templates_${YEAR}_1d_jecpt.root --year ${YEAR}
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_jec_up --JEC "pt"
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_jecpt_jec_down --JEC "pt"

        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_unfolding --JEC "pt&mJ" --unfolding
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding --JEC "pt" --unfolding
        
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_jec_up --JEC "pt" --unfolding
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_jec_down --JEC "pt" --unfolding

        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_up.coffea -o $outdir/templates_${YEAR}_1d_unfolding_jec_up --JEC "pt&mJ" --unfolding
        # ./create_root_templates.py -i $indir/templates_${YEAR}_jec_down.coffea -o $outdir/templates_${YEAR}_1d_unfolding_jec_down --JEC "pt&mJ" --unfolding
    done
}

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

function flatten_templates {
    VAR=${1:-none}

    if [ "$VAR" == "none" ];then
        echo "You have to provide a variation name!"
        exit -1
    fi

    for YEAR in ${YEARS[@]};
    do
    if [ "$VAR" == "nominal" ];then
        # JEC on pt
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt --JEC "pt"
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding --JEC "pt" --unfolding

        # JEC on pt&mJ
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d --JEC "pt&mJ"
        ./create_root_templates.py -i $indir/templates_${YEAR}.coffea -o $outdir/templates_${YEAR}_1d_unfolding --JEC "pt&mJ" --unfolding
    elif [ "$VAR" == "toppt_off" ];then
        echo "flattening ${YEAR} templates for variation of ${VAR}"
        # JEC on pt
        ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${VAR} --JEC "pt"
        ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_${VAR} --JEC "pt" --unfolding

        # JEC on pt&mJ
        ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}.coffea -o $outdir/templates_${YEAR}_1d_${VAR} --JEC "pt&mJ"
        ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}.coffea -o $outdir/templates_${YEAR}_1d_unfolding_${VAR} --JEC "pt&mJ" --unfolding
    else
        for DIR in up down;
        do
            echo "flattening ${YEAR} templates for variation of ${VAR}"
            # JEC on pt
            # ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${VAR}_up --JEC "pt"
            ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_${DIR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${VAR}_${DIR} --JEC "pt"

            # ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_${VAR}_up --JEC "pt" --unfolding
            ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_${DIR}.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_${VAR}_${DIR} --JEC "pt" --unfolding

            # JEC on pt&mJ
            # ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_${VAR}_up --JEC "pt&mJ"
            ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_${DIR}.coffea -o $outdir/templates_${YEAR}_1d_${VAR}_${DIR} --JEC "pt&mJ"

            # ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_unfolding_${VAR}_up --JEC "pt&mJ" --unfolding
            ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_${DIR}.coffea -o $outdir/templates_${YEAR}_1d_unfolding_${VAR}_${DIR} --JEC "pt&mJ" --unfolding
        done
    fi
    done
}

# run the things

# flatten_nominal

flatten_template jec
flatten_template fsr
flatten_templates isr
flatten_templates toppt_off
