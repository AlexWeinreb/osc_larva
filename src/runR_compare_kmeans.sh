#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=compare_kmeans
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --time=15:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "~~~~ Loading R script, $(date) ~~~~"
ml R

Rscript R/compare_kmeans.R

echo "~~~~ done $(date) ~~~~"

