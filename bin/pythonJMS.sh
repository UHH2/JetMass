#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
VENVDIR=$SCRIPT_DIR/pythonJMS
VENVPYTHON=$VENVDIR/bin/python

if [ ! -f $VENVPYTHON ]; then
    echo "py3 venv for JetMass unfolding scripts execution does not exists. creating it now..."
    # sometimes CMSSW appends to LD_LIBRARY_PATH and breaks local python3 install. 
    env -u LD_LIBRARY_PATH /usr/bin/python3 -m venv $VENVDIR
    source $VENVDIR/bin/activate
    env -u LD_LIBRARY_PATH pip install -r $SCRIPT_DIR/requirements.txt
    deactivate
fi

# in case LD_LIBRARY is needed inside py3 venv execution we "rename" it here
export LD_LIBRARY_PATH_BACKUP=$LD_LIBRARY_PATH
env -u LD_LIBRARY_PATH "$VENVPYTHON" "$@"