#!/bin/bash
ENV=${1:-"nyvaclab"}
conda create -y -n $ENV
conda install -y -n $ENV -c conda-forge r r-essentials r-irkernel r-ape=5.0
conda install -y -n $ENV -c bioconda bioconductor-dada2
conda install -y -n $ENV -c bioconda bioconductor-phyloseq
conda install -y -n $ENV -c cramjaco r-mirkat
# source activate $ENV
# Rscript R_Setup.R
# source deactivate
