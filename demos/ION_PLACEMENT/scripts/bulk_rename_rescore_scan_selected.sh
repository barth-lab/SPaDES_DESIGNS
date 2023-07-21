#!/bin/bash

read -p "WARNING: creates files names [CTRL+C to cancel] " null

for dir in *_{wt,des}/rescore_scan_selected
do
    count=0
    for pdb in ${dir}/*.pdb #/cluster.cluster*.pdb
    do
        cp ${pdb} ${dir}/$(printf "%04d" ${count}).pdb
        count=$((count+1))
    done
done

