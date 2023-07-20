#!/bin/bash

read -p "WARNING: About to start multiple jobs on cluster [CTRL+C to cancel] " null
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for i in active_des active_wt inactive_des inactive_wt
do
    pdb=$(basename $i)
    for idx in 0 1 2 3 4 5 6 7 8 9
    do
        jobname=${pdb}_${idx}
        sbatch -J ${jobname} -D $(pwd) --export=pdb=${pdb},idx=${idx} ${SCRIPT_DIR}/rescore_grid_scan.slurm
    done
done

