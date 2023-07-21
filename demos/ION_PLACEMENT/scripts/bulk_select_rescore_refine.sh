#!/bin/bash

read -p "WARNING: will overwrite old files [CTRL+C to cancel] " null

origPwd=$(pwd)
for dir in *_{wt,des}/
do
    cd ${dir}
    counter=0
    touch rescore_refine_selected
    rm -Rf rescore_refine_selected
    mkdir -p rescore_refine_selected
    cd ./rescore_refine/
    for idx in $(seq 0 9)
    do
        $(${origPwd}/get_lowest_energy_structure.sh ${idx})
        while read -r line;
        do
            cp ${origPwd}/${dir}/${line}.pdb ${origPwd}/${dir}/rescore_refine_selected/$(printf '%05d' $counter).pdb;
            counter=$((counter+1))
        done < best_filenames${idx}.txt
    done
    cd ${origPwd}
done

