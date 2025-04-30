#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=comp_clust
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --time=23:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "~~~~ Loading R script, $(date) ~~~~"
ml R

Rscript R/compare_clustering_methods.R

echo "~~~~ done ~~~~"

