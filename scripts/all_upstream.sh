#!/bin/bash

bash demultBothPlates.sh

Rscript dada2work-March2018Run.R

Rscript dada2taxonomy-March2018Run.R

Rscript makeTree.R
