#!/bin/bash

jec=${1:-nominal}
cmd_prefix="python jetmass.py -M jms "
cmd_postfix=${2:-"--exe"}
if [ "$jec" = "nominal" ];then
   echo "jec nominal"
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/nominal/TTBar_{YEAR}.py --window 1 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/nominal/WJets_{YEAR}.py --window 2 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/nominal/Combined_{YEAR}.py --window 3 $cmd_postfix
elif [ "$jec" = "up" ];then
   echo "up"
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/TTBar_{YEAR}_JEC_up.py --window 1 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/WJets_{YEAR}_JEC_up.py --window 2 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/Combined_{YEAR}_JEC_up.py --window 3 $cmd_postfix
elif [ "$jec" = "down" ]; then
  echo "down"
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/TTBar_{YEAR}_JEC_down.py --window 1 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/WJets_{YEAR}_JEC_down.py --window 2 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/JECVar/Combined_{YEAR}_JEC_down.py --window 3 $cmd_postfix
elif [ "$jec" = "none" ]; then
  echo "noJEC (on msd)"
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/noJEC/TTBar_{YEAR}_noJEC.py --window 1 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/noJEC/WJets_{YEAR}_noJEC.py --window 2 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix}" --pattern configs/noJEC/Combined_{YEAR}_noJEC.py --window 3 $cmd_postfix
elif [ "$jec" = "nuis" ]; then
  echo "JEC var as nuisance"
   ./tmux_util.py -p "${cmd_prefix} --JECVar " --pattern configs/nominal/TTBar_{YEAR}.py --window 1 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix} --JECVar " --pattern configs/nominal/WJets_{YEAR}.py --window 2 $cmd_postfix
   ./tmux_util.py -p "${cmd_prefix} --JECVar " --pattern configs/nominal/Combined_{YEAR}.py --window 3 $cmd_postfix
else
  echo "option ${jec} not implemented"
fi
     
