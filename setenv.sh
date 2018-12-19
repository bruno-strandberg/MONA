#! /bin/bash

#==============================================================
# parsing from https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash-script/7069755
#==============================================================

if [ -z $JPP_DIR ] || [ -z $AADIR ]; then
    echo "ERROR! NMH/setenv.sh Jpp and aanet not found"
    return 1
fi

if [ -z $OSCPROBDIR ]; then
    echo "ERROR! NMH/setenv.sh environment variable OSCPROBDIR is not set"
    return 1
fi

#==============================================================
# NMHDIR is used in the code to avoid the use of relative path. NMHDIR extraction trick from https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within and Swim. Directories with shares libs are added to LD_LIBRARY_PATH. CPATH allows one to do anywhere eg "#include "AtmFlux.h"
#==============================================================

export NMHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && readlink -e ./ )"

# only add to path if the directory is not already included

if ! echo $LD_LIBRARY_PATH | grep -q $NMHDIR/common_software; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$NMHDIR/common_software
    echo "NMH/setenv.sh::added $NMHDIR/common_software to LD_LIBRARY_PATH"
fi

if ! echo $LD_LIBRARY_PATH | grep -q $NMHDIR/fitter_software; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$NMHDIR/fitter_software
    echo "NMH/setenv.sh::added $NMHDIR/fitter_software to LD_LIBRARY_PATH" 
fi

if ! echo $LD_LIBRARY_PATH | grep -q $OSCPROBDIR; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$OSCPROBDIR
    echo "NMH/setenv.sh::added $OSCPROBDIR to LD_LIBRARY_PATH" 
fi

if ! echo $CPATH | grep -q $NMHDIR/common_software; then
    export CPATH=${CPATH}:$NMHDIR/common_software
    echo "NMH/setenv.sh::added $NMHDIR/common_software to CPATH" 
fi

if ! echo $CPATH | grep -q $NMHDIR/fitter_software; then
    export CPATH=${CPATH}:$NMHDIR/fitter_software
    echo "NMH/setenv.sh::added $NMHDIR/fitter_software to CPATH" 
fi

if ! echo $CPATH | grep -q $OSCPROBDIR; then
    export CPATH=${CPATH}:$OSCPROBDIR
    echo "NMH/setenv.sh::added $OSCPROBDIR to CPATH" 
fi

echo "NMH/setenv.sh::finished" 
