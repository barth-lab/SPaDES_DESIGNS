#!/bin/bash
set -e
if [[ -z $3 ]]
then
    echo "Usage: $0 pdb_without_ion pdb_single_selected_ion index_label"
    echo "  pdb_without_ion = structure with NO ion in pdb"
    echo "  pdb_single_selected_ion = pdb file of protein with ion containing the specific selected ion line" 
    echo "  index_label = suffix label for grid pdbs"
    exit
fi

pdb=$1
pdb_ion=$2
label=$3

firstline_pdb=$(cat $pdb | head -1)
lastline_pdb_ion=$(cat $pdb_ion | grep ' NA    NA')
atomNum=$(printf '% 5s' $(echo "${firstline_pdb:7:5} - 1" | bc -l))
resNum=$(printf '% 4s' $(echo "${firstline_pdb:22:4} - 1" | bc -l))
x1=${lastline_pdb_ion:30:8}
y1=${lastline_pdb_ion:38:8}
z1=${lastline_pdb_ion:46:8}

mkdir -p grid_refine
count=0
for xshift in -0.20 -0.15 -0.10 -0.05 0.00 0.05 0.10 0.15 0.20
do
    for yshift in -0.20 -0.15 -0.10 -0.05 0.00 0.05 0.10 0.15 0.20
    do
        for zshift in -0.20 -0.15 -0.10 -0.05 0.00 0.05 0.10 0.15 0.20
        do
            count=$((count+1))
            idx=$(printf "%05d" ${count})

            xnew1=$(printf '% 8s' $(printf '%.03f' $(echo "$x1 + $xshift" | bc -l)))
            ynew1=$(printf '% 8s' $(printf '%.03f' $(echo "$y1 + $yshift" | bc -l)))
            znew1=$(printf '% 8s' $(printf '%.03f' $(echo "$z1 + $zshift" | bc -l)))
            indexed_pdb=grid_refine/$(basename ${pdb} .pdb)_${label}_${idx}.pdb
            cp ${pdb} ${indexed_pdb}
            #echo "HETATM${atomNum} NA   NA  A${resNum}    ${xnew1}${ynew1}${znew1}  1.00 20.00          NA" >> ${indexed_pdb}
            echo -e "HETATM   1   NA  NA  A${resNum}    ${xnew1}${ynew1}${znew1}  1.00 20.00          NA\n$(cat ${indexed_pdb})" > ${indexed_pdb}
        done
    done
done
