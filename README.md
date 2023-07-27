# SPaDES Designs

Within this repo you will find information pertaining to running Rosetta SPaDES based on our recent paper TBC. 

We include specific installation instructions, as well as detailed demos - specifically relating to the ion sampling and designs discussed in the article. In particular, the design based demo is transferable to any SPaDES hydrate design problem.

**IMPORTANT TO NOTE** The hydrate movers, that is those compatabile with XML scripts and those needed to interface SPaDES with phenix, are available on the other repo (https://github.com/barth-lab/phenix_with_spades). These movers and the subsequent changes to protocols files etc. within Rosetta are entirely compatabile with the functions provided here, as well as the most up to date version of Rosetta. In other words, feel free to compile a version of Rosetta with the changes listed here or with everything as given in the other repo.

## Installation

With the exception of phenix based structure refinement, which is described in a separate repo (https://github.com/barth-lab/phenix_with_spades), the work described in the article was entirely obtained with Rosetta 2022. 

The specific hydrate functions and updated options file are entirely compatabile with the up-to-date version of Rosetta (as of July 2023), however they are yet to be pushed into the main branch of Rosetta. We offer here those functions to copy into your own version of Rosetta.

First, you will need to download the latest Rosetta build (https://www.rosettacommons.org/software/license-and-download), or create a new version of your current Rosetta installation for this install.

If needed, move into the Rosetta directory as you prefer and decompress:
```
tar -xvf rosetta.tar.gz
```
The hydrate files, i.e. those files in SPaDES_DESIGNS/hydrate_files need to be moved into:
```
/path/to/rosetta/source/src/protocols/hydrate
```
While the options file needs to be moved into:
```
/path/to/rosetta/source/src/basic/options
```
After this, you can compile as normal. For example:
```
cd /path/to/rosetta/source
./scons.py bin mode=release -j 10
```
Where j=10 refers to the number of CPUs you want to compile with.

## Running the demos

We offer two demos to recreate some of the data we produced for the article, specifically based on the designed High10 example. We split these into the DESIGN and ION_SAMPLING demos, where the DESIGN example is specific to using SPaDES for design, i.e. facilitating water network optimisation during Rosetta based design. The ION_SAMPLING is based on the DESIGN approach, but more expansive as it involves sampling ion placement within the NA binding site of A2AR and optimisation of residues/water around said site.

The DESIGN example in particular is designed to be transferable, i.e. one can take the scripts used within and apply them to their own design problems. 

The ION_SAMPLING realistically requires running on a supercomputing cluster owing to the sheer number of ion sites and possible rotamer/water network sampling. We offer general slurm scripts both to run the ION SAMPLING demo, and to generally run bulk design.

There are detailed README.md files in each of the demo folders containing instructions to run the demos.

## All Inputs

The general inputs used for the rest of the design/ion sampling results in the paper can be found in the all_input folder.

The files are set out in such a way that they can be inserted into the demos and run in place of the HIGH10 example. The relationship between the LOW/HIGH numbering and the actual mutations from WT are given in a README file.

## Contact

If you have any issues with the installation, demo or otherwise, please either submit a ticket above or contact me at lucas.rudden@epfl.ch
