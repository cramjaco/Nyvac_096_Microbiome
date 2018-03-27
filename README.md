Nyvac096Microbiome

Jacob A. Cram
jcram@fhcrc.org
cramjaco@gmail.com

This directory accompanies the (in prep as of this writing) manuscript "The human gut microbiota associated with baselin eimmune status and response to HIV vaccines". It contains materials to run both the upstream processing of the microbiome data (demultiplexing with Qiime, sequence variant assignment with DADA2, phylogenetic tree with phangorn, SV taxonomic assignment with DADA2) and downstream analysis (statistics in a jupyter notebook file in R).

The downstream analysis can be run without re-doing the upstream portion. We default to using the files generated by the upstream analysis from our initial run, for consistancy between runs. There appears to be some variability in the results that one gets between upstream runs.


# Dependencies:
## For upstream analysis.
 * A running version of Anacondas.
 * To run the demultiplexing, you will need to install qiime1. I recommend using anacondas to set up the following environment
`# conda create -n qiime1 numpy=1.10 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda`

 * To run dada2, you will need R (I am usign version 3.4.1).
 ** You should be able to do this with `conda install -c r r-essentials r-irkernel`

 ** Within R you will need to install, with install.packages()
 *** dada2 (version 1.4.0) and
 *** ggplot2 (I am using 2.2.1.9000, you can probably get away with using a near by normal version thereof, I ended up installing the special version trying to work around some other problem and never bother to revert.)'

 * For the phylogenetic tree you will need to install the following, within R with install.packages
 ** tidyverse
 ** phangorn
 ** DECIPHER

## For downstream analyis:
 * jupyter notebook or jupyter lab running in condas
 *  a working R kernel for jupyter notebook (we used version 3.4.1).
 ** You should be able to do this with `conda install -c r r-essentials r-irkernel` (If you haven't already)
 * View the notebook file and install every package in the __Library__'s section by opening an R terminal and running install.packages() for each one not already present on your system. Some of them require bioconductor.

# How to run analyses
## Upstream analysis

The order for this is:
* demultiplex
* call SVs
* make tree
* generate taxonomic information.

These scripts can be called in order by calling, from the `scripts\` directory
`all_upstream.sh`

Individual pieces can be run as follows:

 * To demultiplex, run `sh scripts/demultBothPlates.sh`

 * To remake dada2 sequence varients run `Rscript dada2work-March2018Run.R`. One can also open the r script and run it in any R interpreter. (This is true of all of the subsequent R steps. Such a process makes for substantially easier debugging.

 * To make the phylogenetic tree `Rscript makeTree.R`

 * To generate taxonomic information first acquire necessary training data by running `sh pull_training.sh`. Then run `Rscript dada2taxonomy-March2018Run.R`.

## Downstream analysis.

This can be run independently of the downstream analysis. It defaults to using data from the `data\` directory. In theory, all one should need to do is open the Mar2018_096.ipynb file in jupyter notebook or jupyter lab and run all of the cells.