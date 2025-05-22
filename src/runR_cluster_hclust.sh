#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=compare_kmeans
#SBATCH -c 1
#SBATCH --mem=130G
#SBATCH --time=30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "~~~~ Loading R script, $(date) ~~~~"
ml R

Rscript R/cluster_cellgenes.R

echo "~~~~ done $(date) ~~~~"

