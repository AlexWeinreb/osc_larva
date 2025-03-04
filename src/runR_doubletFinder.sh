#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=doubletFinder
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --time=5:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "~~~~ Loading R script, $(date) ~~~~"
ml R

Rscript R/doubletFinder.R

echo "~~~~ done ~~~~"

