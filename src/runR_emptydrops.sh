#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=emptydrops
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "Loading R script, $(date)"
ml R

Rscript R/emptydrops.R

echo "~~~~ done ~~~~"

