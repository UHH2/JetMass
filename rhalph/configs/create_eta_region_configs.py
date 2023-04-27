#!/usr/bin/env python3
import os

years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
samples = ["WJets", "TTBar", "Combined"]
jec_vars = ["", "_JEC_down", "_JEC_up", "_noJEC"]

for eta_region in ["barrel", "endcap"]:
    for year in years:
        for sample in samples:
            for jec_var in jec_vars:
                config_name = f"{sample}_{year}{jec_var}.py"
                config_name_eta_region = f"eta_regions/{sample}_{year}{jec_var}_{eta_region}.py"
                os.system(f"ls -la {config_name}")
                os.system(f"cp {config_name} {config_name_eta_region}")
                os.system(f'sed -i "s/\\(templates.*\\).root/\\1_{eta_region}.root/g" {config_name_eta_region}')
                os.system(
                    f'sed -i "s/\\"\\(\(VJets\|TTBar\|Combined\){year}.*\\)\\"/\\"\\1{eta_region.upper()}\\"/g" {config_name_eta_region}' # noqa
                )
