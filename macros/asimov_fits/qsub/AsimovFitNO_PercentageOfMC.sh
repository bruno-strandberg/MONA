#!/bin/bash

inputnr=$1
directory="$NMHDIR"

echo "input: ${inputnr}"

command="root -q ${directory}/macros/asimov_fits/cross_check_macros/AsimovFitNO_PercentageOfMC.C\+\(${inputnr}\)"
echo $command && eval $command

