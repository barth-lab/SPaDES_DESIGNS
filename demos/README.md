# INFORMATION

Within this folder you will find two demos to run the ion sampling with Rosetta SPaDES and a design with SPaDES. 

The design is more Rosetta "native", that is, like other Rosetta functions, it is effectively a single command line with various flags of your choosing to be included. I will outline each of the files included and it's purpose within the design, for more general information, please see https://new.rosettacommons.org/docs/latest/full-options-list#hydrate. 

The ion sampling requires some intermediate steps as it involves both an initial coarse search, followed by geometric clustering, then a more refined search to get the final results. The ion sampling is embarrassingly parallelisable, that is each ion position can be run individually, therefore we structure the SPaDES runs as being run on a large cluster (specifically a slurm based cluster) to maximise generation of results in the shortest timeframe possible.

There are additional README.md files in each of the demo folders containing detailed instructions to run the demos.

The general inputs used for the rest of the results discussed in the paper can be found in the previous ../all_input folder.
