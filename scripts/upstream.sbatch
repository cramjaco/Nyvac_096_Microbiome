#!/bin/bash
#SBATCH --partition largenode
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --exclusive
#SBATCH --time 24:00:00
#SBATCH --mem-per-cpu 16000
#SBATCH --job-name nyvac-upstream
#SBATCH --output nyvac-upstream-log-%J.txt

bash all_upstream.sh