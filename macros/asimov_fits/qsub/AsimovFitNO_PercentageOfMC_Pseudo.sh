#!/bin/bash

inputnr=$1
directory="$MONADIR"

echo "input: ${inputnr}"

command="root -q ${directory}/macros/asimov_fits/cross_check_macros/AsimovFitNO_PercentageOfMC_Pseudo.C\+\(${inputnr}\)"
echo $command && eval $command

