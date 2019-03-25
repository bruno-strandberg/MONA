#! /bin/bash

#==============================================================
# check that JPP, aanet and OSCPROB directories are set
#==============================================================

if [ -z $JPP_DIR ] || [ -z $AADIR ]; then
    echo "ERROR! MONA/setenv.sh Jpp and aanet not found"
    return 1
fi

if [ -z $OSCPROBDIR ]; then
    echo "ERROR! MONA/setenv.sh environment variable OSCPROBDIR is not set"
    return 1
fi

#==============================================================
# MONADIR is used in the code to avoid the use of relative path.
# Directories with shares libs are added to LD_LIBRARY_PATH.
# CPATH allows one to do anywhere eg "#include "AtmFlux.h", useful in ROOT macros
#==============================================================

export MONADIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && readlink -e ./ )"

# only add to path if the directory is not already included

if ! echo $LD_LIBRARY_PATH | grep -q $MONADIR/common_software; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$MONADIR/common_software
    echo "MONA/setenv.sh::added $MONADIR/common_software to LD_LIBRARY_PATH"
fi

if ! echo $LD_LIBRARY_PATH | grep -q $MONADIR/fitter_software; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$MONADIR/fitter_software
    echo "MONA/setenv.sh::added $MONADIR/fitter_software to LD_LIBRARY_PATH" 
fi

if ! echo $LD_LIBRARY_PATH | grep -q $OSCPROBDIR; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$OSCPROBDIR
    echo "MONA/setenv.sh::added $OSCPROBDIR to LD_LIBRARY_PATH" 
fi

if ! echo $CPATH | grep -q $MONADIR/common_software; then
    export CPATH=${CPATH}:$MONADIR/common_software
    echo "MONA/setenv.sh::added $MONADIR/common_software to CPATH" 
fi

if ! echo $CPATH | grep -q $MONADIR/fitter_software; then
    export CPATH=${CPATH}:$MONADIR/fitter_software
    echo "MONA/setenv.sh::added $MONADIR/fitter_software to CPATH" 
fi

if ! echo $CPATH | grep -q $OSCPROBDIR; then
    export CPATH=${CPATH}:$OSCPROBDIR
    echo "MONA/setenv.sh::added $OSCPROBDIR to CPATH" 
fi

if ! echo $PATH | grep -q $MONADIR/apps/cmdapps; then
    export PATH=${PATH}:$MONADIR/apps/cmdapps
    echo "MONA/setenv.sh::added $MONADIR/apps/cmdapps to PATH" 
fi

echo "MONA/setenv.sh::finished" 
