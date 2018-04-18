#! /bin/bash

#==============================================================
# Determine whether to compile agains aanet or not
# parsing from https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash-script/7069755
#==============================================================

aanet='false'

while getopts 'a' flag; do
  case "${flag}" in
    a) aanet='true' ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

#==============================================================
# NMHDIR is used in the code to avoid the use of relative path, NMHDIR extraction trick from https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within and Swim
# LD_LIBRARY_PATH handy when compiled against libnmhsoft.so
# CPATH allows one to do anywhere eg "#include "AtmFlux.h"
#==============================================================

export NMHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export OSCPROBDIR="/home/bstrand/Code/OscProb"

# only add to path if the directory is not already included

if ! echo $LD_LIBRARY_PATH | grep -q $NMHDIR/common_software; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$NMHDIR/common_software
fi

if ! echo $CPATH | grep -q $NMHDIR/common_software; then
    export CPATH=${CPATH}:$NMHDIR/common_software
fi

if ! echo $CPATH | grep -q $OSCPROBDIR; then
    export CPATH=${CPATH}:$OSCPROBDIR
fi

if ($aanet -eq 'true'); then
    echo "Set variables for NMH analysis, set to link against AANET"
    if ! echo $CPATH | grep -q $AADIR; then
	export CPATH=${CPATH}:${AADIR}:${AADIR}evt
    fi
    export USE_AANET=true
else
    echo "Set variables for NMH analysis"
    export USE_AANET=false
fi

