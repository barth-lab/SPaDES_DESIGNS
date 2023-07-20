#!/bin/bash

read -p "WARNING: overwrites grid pdb files [CTRL+C to cancel] " null

# setup grid scan methods
origPwd=$(pwd)
input=${origPwd}/inputs
scripts_folder=${origPwd}/scripts

for x in active_des active_wt inactive_des inactive_wt; do
    state=$(echo $x | cut -d'_' -f1)
    cd $x
    ${origPwd}/cubic_grid_scan.sh ${input}/${x}.pdb ${input}/sodium_ion_${state}.pdb
    cd ${origPwd}/
done

