#! /bin/bash

#==============================================================
# Determine whether to compile agains aanet or not
# parsing from https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash-script/7069755
#==============================================================

aanet='false'
oscprob=''

while getopts 'ao:' flag; do
  case "${flag}" in
    a) aanet='true' ;;
    o) oscprob="${OPTARG}" ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

oscprob="$(readlink -e $oscprob)"

#==============================================================
# NMHDIR is used in the code to avoid the use of relative path, NMHDIR extraction trick from https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within and Swim
# LD_LIBRARY_PATH handy when compiled against libnmhsoft.so
# CPATH allows one to do anywhere eg "#include "AtmFlux.h"
#==============================================================

export NMHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export OSCPROBDIR="$oscprob"

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
    if ! echo $CPATH | grep -q $AADIR; then
	export CPATH=${CPATH}:${AADIR}:${AADIR}evt
    fi
    export USE_AANET=true
else
    export USE_AANET=false
fi

echo "NMH/setenv.sh::compiling against aanet --> $aanet" 
echo "NMH/setenv.sh::OscProb directory --> $oscprob" 
echo "NMH/setenv.sh::finished" 
