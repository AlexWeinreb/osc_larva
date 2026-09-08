module load SAMtools
module load miniconda
conda activate velocyto

echo "Start $(date)"

dir="/vast/palmer/scratch/hammarlund/aw853/250331_align"
gtf="/gpfs/ycga/work/hammarlund/aw853/references/WS295/c_elegans.PRJNA13758.WS295.canonical_geneset_GFP.gtf"

sample=$1

echo "-------   Starting: $sample at $(date)"


if [ ! -d "$dir/$sample/" ]
then
  echo "Sample not found: $sample"
  echo "at path $dir/$sample/"
  exit 1
fi





velocyto run10x -@ $SLURM_CPUS_PER_TASK \
	-vvv \
	$dir/$sample/ \
	$gtf
  



echo "End $(date)"

