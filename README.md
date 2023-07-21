# Introduction

Within this repo you will find information pertaining to running Rosetta SPaDES based on our recent paper TBC. 

We include specific installation instructions, as well as detailed demos - specifically relating to the ion sampling and designs discussed in the article. In particular, the design based demo is transferable to any SPaDES hydrate design problem.

# Installation


# Demo descriptions

The design is more Rosetta "native", that is, like other Rosetta functions, it is effectively a single command line with various flags of your choosing to be included. I will outline each of the files included and it's purpose within the design, for more general information, please see https://new.rosettacommons.org/docs/latest/full-options-list#hydrate.

The ion sampling requires some intermediate steps as it involves both an initial coarse search, followed by geometric clustering, then a more refined search to get the final results. The ion sampling is embarrassingly parallelisable, that is each ion position can be run individually, therefore we structure the SPaDES runs as being run on a large cluster (specifically a slurm based cluster) to maximise generation of results in the shortest timeframe possible.

There are detailed README.md files in each of the demo folders containing instructions to run the demos.

# All Inputs

The general inputs used for the rest of the results discussed in the paper can be found in the previous all_input folder.

All structural inputs used within the context of the manuscript
