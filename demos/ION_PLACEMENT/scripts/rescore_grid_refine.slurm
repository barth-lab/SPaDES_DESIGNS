#!/bin/bash
#SBATCH --time=72:00:00 -n 1 
#SBATCH --cpus-per-task=1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --partition=serial
#SBATCH --job-name=hyd_sodium

rand=$$
origPwd=$(pwd)
inputs="${origPwd}/inputs/"

rand=$$
origPwd=$(pwd)

if [[ -z $pdb ]] || [[ -z $idx1 ]] || [[ -z $idx2 ]]
then
    echo "Usage: $0 pdb idx1 idx2"
    echo "  script meant to be run via bulk_rescore*.sh"
    exit
fi

hydrate="/path/to/rosetta/main/source/bin/hydrate.linuxgccrelease"
export LD_LIBRARY_PATH=/path/to/rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/
DATABASE="/path/to/rosetta/main/database"

cd ${pdb}
mkdir -p rescore_refine
$hydrate @${inputs}/rehydrate_flags \
     -in:file:s grid_refine/${pdb}*_???{${idx1},$((idx1+5))}${idx2}.pdb \
     -hydrate:hyfile ${inputs}/${pdb}.hyfile \
     -packing:resfile ${inputs}/${pdb}.resfile \
     -hydrate:ignore_fa_sol_at_positions $(cat ${inputs}/${pdb}.nofasol) \
     -in:file:spanfile ${inputs}/STRUCTURE.span \
     -in:file:extra_res_fa ${inputs}/${pdb}_ligand.fa.params \
     -in:file:extra_res_cen ${inputs}/${pdb}_ligand.cen.params \
     -in:auto_setup_metals \
     -overwrite \
     -constant_seed=1 \
     -seed_offset ${rand} \
     -out:prefix rescore_refine/part_${idx1}${idx2}_${rand}_ >> rescore_refine/part_${idx1}${idx2}_${rand}.log 2>&1

