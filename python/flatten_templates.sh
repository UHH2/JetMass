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

function flatten_variations {
    VAR=${1:-none}

    if [ "$VAR" == "none" ];then
    echo "You have to provide a variation name!"
    exit -1
    fi

    for YEAR in ${YEARS[@]};
    do
    echo "flattening ${YEAR} templates for variation of ${VAR}"
    # JEC on pt
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${VAR}_up --JEC "pt"
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_down.coffea -o $outdir/templates_${YEAR}_1d_jecpt_${VAR}_down --JEC "pt"
    
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_${VAR}_up --JEC "pt" --unfolding
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_down.coffea -o $outdir/templates_${YEAR}_1d_jecpt_unfolding_${VAR}_down --JEC "pt" --unfolding
                                    
    # JEC on pt&mJ
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_${VAR}_up --JEC "pt&mJ"
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_down.coffea -o $outdir/templates_${YEAR}_1d_${VAR}_down --JEC "pt&mJ"

    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_up.coffea -o $outdir/templates_${YEAR}_1d_unfolding_${VAR}_up --JEC "pt&mJ" --unfolding
    ./create_root_templates.py -i $indir/templates_${YEAR}_${VAR}_down.coffea -o $outdir/templates_${YEAR}_1d_unfolding_${VAR}_down --JEC "pt&mJ" --unfolding
    done
}

# run the things

flatten_nominal

flatten_variation jec
flatten_variation fsr
# flatten_variation isr
flatten_variation toppt_off
