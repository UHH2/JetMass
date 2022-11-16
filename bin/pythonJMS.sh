#!/bin/bash

# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# VENVDIR=$SCRIPT_DIR/pythonJMS
# VENVPYTHON=$VENVDIR/bin/python
VENVPYTHON=/nfs/dust/cms/user/albrechs/python/coffea/bin/python

# #ToDo:just use miniconda to install python3.9/3.10
# PYTHON3=$(which python3)
# if [[ $(hostname -s) = naf-cms* ]]; then
#     PYTHON3_INSTALL=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/cms/cmssw/CMSSW_12_0_0/external/slc7_amd64_gcc900
#     PYTHON3=LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHON3_INSTALL/lib;$PYTHON3_INSTALL/bin/python3.9
# fi
# echo "running python3 from ${PYTHON3}"

# if [ ! -f $VENVPYTHON ]; then
#     echo "py3 venv for JetMass unfolding scripts execution does not exists. creating it now..."
#     env -u LD_LIBRARY_PATH $PYTHON3 -m venv $VENVDIR
#     source $VENVDIR/bin/activate
#     env -u LD_LIBRARY_PATH pip install pip --upgrade
#     env -u LD_LIBRARY_PATH pip install -r $SCRIPT_DIR/../python/requirements.txt
#     deactivate
# fi

# sometimes CMSSW appends to LD_LIBRARY_PATH and breaks local python3 install. 
# in case LD_LIBRARY is needed inside py3 venv execution we "rename" it here
export LD_LIBRARY_PATH_BACKUP=$LD_LIBRARY_PATH
env -u LD_LIBRARY_PATH "$VENVPYTHON" "$@"
