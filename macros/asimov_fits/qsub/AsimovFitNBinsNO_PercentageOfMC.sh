#!/bin/bash

inputnr=$1
directory="$MONADIR"

echo "input: ${inputnr}"

command="root -q ${directory}/macros/asimov_fits/cross_check_macros/AsimovFitNBinsNO_PercentageOfMC.C\+\(${inputnr}\)"
echo $command && eval $command

