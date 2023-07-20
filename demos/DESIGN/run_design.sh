#!/bin/bash

# hydratable location
hydrate="/path/to/rosetta/main/source/bin/hydrate.linuxgccrelease"
# Include libraries (your specific compiler may mean this path doesn't exist, change as needed)
export LD_LIBRARY_PATH=/path/to/rosetta/main/source/build/external/release/linux/5.13/64/x86/gcc/9/default/
# Rosetta database location
DATABASE="/path/to/rosetta/main/database/"

rand=$$ # seed for rosetta
currentpwd=$(pwd)

# state = input pdb structure (active or inactive)
state=${1}
# wildtype (wt) or design (des) ?
ptype=${2}

RUN=${3}
touch RUN${RUN}
mkdir RUN${RUN}
cd RUN${RUN}

# assign nofasol value (read file with commas)
nofasol="$(cat ${currentpwd}/inputs/${state}_${ptype}.nofasol)"

# Run hydrate!
$hydrate @${currentpwd}/inputs/general_flags.flags \
	-hydrate:hyfile ${currentpwd}/inputs/${state}_${ptype}.hyfile \
	-packing:resfile ${currentpwd}/inputs/${state}_${ptype}_high10.resfile \
	-in:file:spanfile ${currentpwd}/inputs/STRUCTURE.span \
	-hydrate:ignore_fa_sol_at_positions ${nofasol} \
        -database ${DATABASE} \
        -in:file:s ${currentpwd}/inputs/${state}_wt.pdb  \
	-in:file:extra_res_cen ${currentpwd}/inputs/${state}_ligand.cen.params \
	-in:file:extra_res_fa ${currentpwd}/inputs/${state}_ligand.fa.params \
	-seed_offset ${rand} \
	-out:prefix ./${state}_${ptype}/  > ${state}_${ptype}.log

