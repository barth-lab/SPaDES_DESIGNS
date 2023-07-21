# DESIGN

We take HIGH10 as the example for this demo, which includes the following mutations from WT:

L48A, S91T, L194M, I238Y, V239L, A243L

For the other designs, one only needs to modify appropriately the resfile & potentially the hyfile and repack.

# Inputs

Inside this folder you'll find the various inputs needed to run the four sets of SPaDES runs we run:
 
(1) Active state, high10 design

(2) Active state, WT

(3) Inactive state, high10 design

(4) Inactive state, WT

Feel free to chop and change these. The PDBs are also provided for each state, but of course in theory you can take any PDB and use these files as a basis to setup your design protocol. 

We include the following key files:

(1) X.resfile

Standard resfile to make active design decisions on specific sites with PIKAA. NATAA is used to repack all residues within 5 Ang. of the designed residues and the Na+ ion, while NATRO is used to keep other residues static (and avoid unnecessary computation). It is possible to use the other design keywords.

(2) X.hyfile

All residues specified in this file (based on pose numbering), will be considered hydratable, and SPaDES will attempt to place waters in contact with them, if energetically favourable.

(3) X.nofasol

Residues to avoid implicit solvent contributions to the scoring function. This should correspond to the hyfile residues, as these will have explicit solvent contributions present.

(4) X.pdb

PDB file for the state

(5) X.cen.params / X.fa.params

Rosetta forcefield parameters for included ligand. In a pure protein system, or one in which the ligand parameters are pre-existing within Rosetta, these files wouldn't be necessary. 

(6) X_ligand.cen.pdb / X_ligand.fa.pdb

PDB files for the ligand to correspond to the above parameter files. Again, not necessary outside of this system.

(7) opte_mb_elec.wts

Optimal weights file for SPaDES using pre_talaris_2013 behaviour. Note that general_flags (below) needs to be updated with the direction this this file

(8) STRUCTURE.span

Span file to inform Rosetta of the membrane embedded regions

(9) general_flags.flags

Flag file for general terms applicable to every system, such as repacking information, the scoring function, membrane initialisation, the maximum number of water rotamers to consider (set to 500 default, only change if dealing with a very large hydratable system), hbond threshold to limit the number of new water rotamers that can be added based on the hbond needing to be of this REU strength. Note the folder directory that needs to be updated for opte_mb_elec.wts

# Running the SPaDES designs

run_design.sh can be used to run Rosetta SPaDES

The important things to note in this file is to change the location of Rosetta in the top of the file (i.e. set hydrate, LD_LIBRARY_PATH and DATABASE variables) after installation of Rosetta.

You can run with:

./run_design.sh active des 1

The first argument corresponds to the state (active or inactive), the second whether it is the design or WT (des or wt), while the third simply creates a new folder called RUN1 (RUN2 etc.), within which it will create the active_des folder to place the output PDBs and logfiles.

The important things to play with this file is the actual $hydrate command at the end, which you can change the individual arguments for for your own projects.

# Validation / Screening

After SPaDES has completed it's run, you should end up with 200 individual PDB files (based on the out:nstruct flag in general_flags.flags). Feel free to change this, though note it will take longer. To parallelise, you could make additional runs with different integer numbers for the RUN argument above.

get_best_results.sh can be run in the folder just to simply report on the best structure. It gives the line in the score file corresponding to the best overall rosetta energy.

get_top10_best_result.sh can be run in the folder to get the best 10% of models and copy the output into a new file, top10.sc. Structures notified in file can then be further analysed. 
