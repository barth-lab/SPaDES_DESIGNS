#!/bin/bash

read -p "WARNING: overwrites grid pdb files [CTRL+C to cancel] " null
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

rm -Rf */grid_refine
origPwd=$(pwd)
input=${origPwd}/inputs

for x in active_des active_wt inactive_des inactive_wt; do
    cd $x
    echo $x

    label=0
    for pdb in rescore_scan_selected/0*.pdb
    do
        label=$(basename $pdb .pdb)
        echo $label
        ${SCRIPT_DIR}/cubic_grid_refine.sh ${input}/${x}.pdb ${pdb} ${label}
    done
    cd $origPwd
done
