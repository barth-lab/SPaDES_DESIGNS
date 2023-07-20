#!/bin/bash

if [[ -z $2 ]]
then
    echo "Usage: $0 pdb_without_ion pdb_only_ion"
    echo "  pdb_without_ion = structure with NO ion in pdb"
    echo "  pdb_only_ion = pdb file containing only the relevant ion line" 
    exit
fi

pdb=$1
pdb_ion=$2

lastline_pdb=$(cat $pdb | tail -1)
lastline_pdb_ion=$(cat $pdb_ion | tail -1)
atomNum=$(printf '% 5s' $(echo "${lastline_pdb:7:5} + 1" | bc -l))
resNum=$(printf '% 4s' $(echo "${lastline_pdb:22:4} + 1" | bc -l))
x1=${lastline_pdb_ion:30:8}
y1=${lastline_pdb_ion:38:8}
z1=${lastline_pdb_ion:46:8}

mkdir -p grid_scan

count=0
for xshift in -2.50 -2.00 -1.50 -1.00 -0.50 0.00 0.50 1.00 1.50 2.00 2.50
do
    for yshift in -2.50 -2.00 -1.50 -1.00 -0.50 0.00 0.50 1.00 1.50 2.00 2.50
    do
        for zshift in -2.50 -2.00 -1.50 -1.00 -0.50 0.00 0.50 1.00 1.50 2.00 2.50
        do
            count=$((count+1))
            idx=$(printf "%05d" ${count})

            xnew1=$(printf '% 8s' $(printf '%.03f' $(echo "$x1 + $xshift" | bc -l)))
            ynew1=$(printf '% 8s' $(printf '%.03f' $(echo "$y1 + $yshift" | bc -l)))
            znew1=$(printf '% 8s' $(printf '%.03f' $(echo "$z1 + $zshift" | bc -l)))
            indexed_pdb=grid_scan/$(basename ${pdb} .pdb)_${idx}.pdb
            cp ${pdb} ${indexed_pdb}
            echo "HETATM${atomNum} NA   NA  A${resNum}    ${xnew1}${ynew1}${znew1}  1.00 20.00          NA" >> ${indexed_pdb}
        done
    done
done
