#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=cr_count
#SBATCH -c 24
#SBATCH --mem=80G
#SBATCH --time=6-05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


set -ue

# Define variables
WS="WS295"


dir="/gpfs/gibbs/project/hammarlund/aw853/osc_larva/"



# List the samples to process ----


mapfile -t sample_paths < $dir/data/raw_paths.txt
mapfile -t sample_ids < $dir/data/raw_sample_ids.txt




echo "Processing ${#sample_paths[@]} samples..."


module load CellRanger/8.0.1

for (( index=0; index<${#sample_paths[@]}; index++ ))
do 
	cur_id="${sample_ids[$index]}"
	cur_path="${sample_paths[$index]}"
	
	echo "-------   Starting: $cur_id at $(date)"
	echo "$cur_path"
	
	
	cellranger count \
		--id=$cur_id \
		--fastqs=$cur_path \
		--transcriptome=$dir/intermediates/c_el_${WS}_GFP_genome \
		--create-bam=true \
		--localcores=$SLURM_CPUS_PER_TASK \
		--localmem=$(( $SLURM_MEM_PER_NODE / 1024 ))
	
	
	echo
	echo
	echo "    > alignment done."
	
	


done


echo
echo
echo " ---- Finished at $(date)"



