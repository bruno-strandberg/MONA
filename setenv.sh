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
# NMHDIR is used in the code to avoid the use of relative path
# LD_LIBRARY_PATH may come handy when you wish to compile agains libnmhsoft.so
# CPATH is nice, so that you can do anywhere eg "#include "AtmFlux.h"
# thanks Simon and https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
#==============================================================

export NMHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$NMHDIR/common_software
export CPATH=${CPATH}:$NMHDIR/common_software

if ($aanet -eq 'true'); then
    echo "Set variables for NMH analysis, set to link against AANET"
    export CPATH=${CPATH}:${AADIR}:${AADIR}evt
    export USE_AANET=true
else
    echo "Set variables for NMH analysis"
    export USE_AANET=false
fi

