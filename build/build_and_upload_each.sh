#!/bin/bash

conda config --set anaconda_upload yes
# conda build r-bayesm
conda build r-compositions
conda build r-decipher
conda build r-imputemissings
conda build r-purrrlyr
conda build r-tensora
