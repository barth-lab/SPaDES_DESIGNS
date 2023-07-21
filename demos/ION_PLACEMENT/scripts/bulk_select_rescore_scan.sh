#!/bin/bash

read -p "WARNING: overwrites files [CTRL+C to cancel] " null
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

origPwd=$(pwd)
directory='rescore_scan'
top=100

for pdb in active_des active_wt inactive_wt inactive_des;
do
    echo $pdb
    cd $origPwd
    cd ${pdb}/${directory}
    touch ../${directory}_top${top}
    rm -r ../${directory}_top${top}
    mkdir -p ../${directory}_top${top}/
    grep 'total_score' *.pdb > tmp.txt
    cat tmp.txt | sort -r -n -k2 | tail -${top} | sed 's/:.*$//g' > top${top}.txt
    rm tmp.txt

    count=0
    for i in $(tac top${top}.txt)
    do
        count=$((count+1))
        cp $i ../${directory}_top${top}/$(printf '%04d' $count)_$(basename $(pwd))_$i
    done
    cd ../${directory}_top${top}
    touch top_sodium.txt
    rm top_sodium.txt
    for i in 0*.pdb
    do
        sed -n '/ WAT /!p' $i > sodiumonly_${i}
        echo sodiumonly_${i} >> top_sodium.txt
    done
    cd ../${directory}_top${top}
    touch cluster_test
    rm cluster*
    ${SCRIPT_DIR}/Cluster_ligands_position_scan.py -l top_sodium.txt -t cluster -s ${pdb}
    # remove total from last line in the wc and the first line, which is just cluster.cluster
    # cluster.cluster file contains name of pdb : cluster_id, clusternum contains the PDBs belonging to a cluster 
    wc -l cluster.cluster* | head -n -1 | tail -n+2 | sort -k1 -nr | head -12 | tail -10 | awk '{print $2}' > top_clusters
    touch ../${directory}_selected
    rm -r ../${directory}_selected
    mkdir -p ../${directory}_selected
    for i in $(cat top_clusters)
    do
        cp $(cat $i | head -1).pdb ../${directory}_selected/${i}_$(cat $i | head -1).pdb
    done
    cd $origPwd
done

