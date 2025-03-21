#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=preproc_ptDE
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH --time=20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "Start $(date)"

ml R

Rscript R/step_3_curve_shape.R



echo "End $(date)"



