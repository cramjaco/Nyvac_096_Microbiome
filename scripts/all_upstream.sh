#!/bin/bash

bash demultBothPlates.sh

## Install nyvac-lab-2 conda environment before running this
## conda env create cramjaco/nyvac-lab-2

source activate nyvac-lab-2

Rscript dada2work-March2018Run.R

Rscript dada2taxonomy-March2018Run.R

Rscript makeTree.R

source deactivate
