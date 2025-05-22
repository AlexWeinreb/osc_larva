#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=runR_velocyto_load_loom
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --time=20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "Start $(date)"

ml R

Rscript R/velocyto_load_loom.R



echo "End $(date)"



