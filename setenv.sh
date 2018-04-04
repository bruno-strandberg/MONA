#! /bin/bash

# thanks Simon and
# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within

# the directory is used in the code to avoid the use of relative path
# the LD_LIBRARY_PATH may come handy when you wish to compile agains libnmhsoft.so
# CPATH is nice, so that you can do anywhere eg "#include "AtmFlux.h"

export NMHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$NMHDIR/common_software
export CPATH=${CPATH}:$NMHDIR/common_software

# add aanet to cpath
export CPATH=${CPATH}:${AADIR}:${AADIR}evt

echo "Set variables for NMH analysis"
