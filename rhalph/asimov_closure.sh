#!/bin/bash

MODELDIR=`realpath ${1:-.}`
MODELDIRBASE=`basename ${MODELDIR}`
echo "performing series of closure tests using asimov data in ${MODELDIRBASE}"

if [[ ! -d ${MODELDIR} ]]; then
  echo "provided modeldir (${MODELDIR}) does not exist!"
fi

cd $MODELDIR
source wrapper.sh env

# function join_by { local IFS="$1"; shift; echo "$*"; }
function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf %s "$f" "${@/#/$d}"
  fi
}

function printexec {
  echo $1
  eval $1
}

common_opt="--redefineSignalPOIs $(join_by , ${POIS[@]})"
common_opt+=" --cminDefaultMinimizerStrategy 0"
common_opt+=" --robustFit 1"
common_opt+=" --trackParam $(join_by , ${POIS[@]})"

# perform simple fit with different values for expectSignal
for S in 1 1.1 2; do
  expS=$(join_by =${S}, ${POIS[@]}=${S})
  name="AsimovExpectSignal$(echo $S | sed 's/\./p/g')"
  # printexec "combine -M FitDiagnostics -d model_combined.root -t -1 --setParameters ${expS} -n .${name} ${common_opt}"
  printexec "./../plot_pois.py -f fitDiagnostics.${name}.root -P $(join_by " " ${POIS[@]}) --expectSignal ${S}"
done

# perform fit with all POIs set to one except for one
expS=$(join_by =1, ${POIS[@]}=1)
# for POI in ${POIS[@]}; do
for S in 0.9, 1.1; do
  for POI in ${POIS[@]}; do
    name="AsimovExpectSignal${POI}_$(echo $S | sed 's/\./p/g')"
    printexec "combine -M FitDiagnostics -d model_combined.root -t -1 --setParameters ${expS},${POI}=${S} -n .${name} ${common_opt}"
    printexec "./../plot_pois.py -f fitDiagnostics.${name}.root -P $(join_by " " ${POIS[@]}) --expectSignal 1.0"
  done
done

for S in 2 -2; do
  for POI in ${POIS[@]}; do
    name="AsimovExpectSignal${POI}_$(echo $S | sed 's/-/n/g')"
    printexec "./../plot_pois.py -f fitDiagnostics.${name}.root -P $(join_by " " ${POIS[@]}) --expectSignal 1.0"
  done
done
