# ION_PLACEMENT

We will take HIGH10 as the example for this demo, which includes the following mutations:

L48A, S91T, L194M, I238Y, V239L, A243L

For the other designs, one only needs to modify appropriately the resfile & potentially the hyfile and repack (see all_input for more details).

The majority of needed files and scripts are similar to the DESIGN demo, with some minor differences centered around the need to sample different sodium positions in the systems. When it comes to runtime, the main Rosetta SPaDES command is identical to that shown in DESIGN. The differences, therefore, are in the pre- and post-processing.

# Inputs

Inside this folder you'll find the various inputs needed to run the four sets of SPaDES runs we run: 

(1) Active state, high10 design

(2) Active state, WT

(3) Inactive state, high10 design

(4) Inactive state, WT

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

(10) sodium_X.pdb

Sodium PDB file to insert into our PDB structures during the grid sampling

# Running the Ion Sampling routine

The order of input is given below, as well as a brief description of what each script is doing. Note that some scripts are run on a supercomputing cluster in order to parallelise the operation to cover the 1000+ structures that need to be modelled.

As a word of warning, keep an eye on your hard drive disk space as you run through this process. The bulk generation of PDBs needed for the ION SAMPLING can be quite demanding (several GB), and if you are running this on many systems this can become problematic.

Run from the ION_PLACEMENT base folder

(1) Prep input for initial ion grid scan placement

./scripts/bulk_cubic_grid_scan.sh

This loops through each system (the 4 listed above) using ./scripts/cubic_grid_scan.sh and places the NA in a grid pattern around its current position in the host PDB file. Note this first phase is very coarse, and thus we move in units of 0.5 Ang. between -2.50 Ang. and 2.50 Ang. around the current site. The NA is always inserted as the first atom in the produced PDBs, and is set to hydratable in the hyfiles. The output PDB files are placed in a new folder called grid_scan in each active_des, active_wt etc. system folder.

This script can take ~5 minutes to run.

(2) Copy this folder over to your cluster/folder of choice.

The scripts used to run SPaDES for the ION SAMPLING are slurm scripts designed for our host cluster. Please change the parameters as needed. 

(3) Run SPaDES hydrate on your cluster

./scripts/bulk_rescore_grid_scan.sh

This again loops through each system, and then subsequently splits the files for each into 10 batches and runs hydrate on each batch independently, i.e. it submits 40 jobs total. Specifically, it submits the ./scripts/rescore_grid_scan.slurm file, which you will need to change as necessary. This file structure follws the general input seen in the DESIGN demo, so please refer to that for more information. Specifically, you will likely need to update the SBATCH headings at the top, and you will certainly need to update the hydrate, LD_LIBRARY_PATH and DATABASE variables as necessary depending on where your Rosetta is installed.

This run normally takes ~24-48 hours to complete

(4) Copy the files back to your home computer

(5) Cluster the resultant PDBs

bulk_select_rescore_scan.sh

This step geometrically clusters the Na positions and sidechain conformations to maximise the diversity of samples for the refinement stage while ensuring the energies are reasonable. You should end up with ~30 clustered results after this run.

Specifically, it runs the Cluster_ligands_position_scan.py which relies on a series of clustering scripts given in the clustering_scripts folder. It is currently setup to be applied to the active and inactive states discussed in this demo. While the hyperparameters on the RMSD bounds (ClusterContact.Iterate_Cluster_vector function in script) are general, you may need to play around with them to get that ~30 sample diversity - though of course this depends on how much diversity you want in your system.

It relies on the following additional modules to be installed:

- numpy
- BioPython

An alternative to the above strategy is to take the top 30 Na scoring structures from the pool of solutions given in rescore_scan and take those forward into step 6. We found this gave similar results to the clustering overall (but may not be generalisable long term).

(6) Rename files

The above clustering should create a new folder called rescore_scan_selected. If you use the second approach mentioned based on the top 30 scoring results, these should be placed in a folder of the same name. Now we run:

bulk_rename_rescore_scan_selected.sh

This just renames the files into a nicer format for the following stages.

(7) Generate a finer set of gridpoints around each NA site

./bulk_cubic_grid_refine.sh

Similar to bulk_cubic_grid_scan, this generates the grid points of Na positions about the sites identified in part 5. This relies on the script

cubic_grid_refine.sh

to run. Again, the NA atom is placed as the first line. This time, NA atoms are placed at 0.05 Ang. intervals in the -0.2 -> 0.2 Ang. grid about the central atom.

This script takes about ~15 minutes to run, largely because we have many more files now that need to be generated.

(8) Copy files back onto your cluster

(9) Generate new hydrated SPaDES models with hydrate

./bulk_rescore_grid_refine.sh

Similar to step 3 where we ran the original coarse hydrate, we will once again be running hydratable design in a manner similar to the DESIGN demo. It launches a slurm script:

./rescore_grid_refine.slurm

which you will need to edit depending on your cluster / rosetta installation specifications. 

Please be aware this will likely launch ~1000 jobs on your cluster partitian, and will take roughly ~24-48 hours, so if you pay for cluster use be aware of this.

(10) Score output

./bulk_select_rescore_refine.sh

Finally, you can score and get the top 10% of results (moved into a folder rescore_refine_selected) using the get_lowest_energy_structure.sh script. You can run this on the cluster or locally depending on hard disk requirements. Alternatively, you can just get the best scored models (see DESIGNS README) if need be. 
