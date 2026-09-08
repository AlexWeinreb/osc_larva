#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=prep_cellranger
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


set -ue



# Prepare annotation with GFP




ref="/gpfs/ycga/work/hammarlund/aw853/references"
WS="WS295"





# Checks

echo "-------   Checks   -------"
if [[ ! -f ${ref}/"nsIs198/nsIs198.gtf" ]]; then
	echo "Error, transgene GTF file not found."
	exit
fi

if [[ ! -f ${ref}/"nsIs198/nsIs198.fa" ]]; then
	echo "Error, transgene FASTA file not found."
	exit
fi

# From Wormbase

if [[ ! -d $ref/$WS ]]; then
	echo "Error: Create the references directory first!"
	exit
fi
echo "    > done."
echo


# define variables

base_url="https://downloads.wormbase.org/releases/"${WS}/"species/c_elegans/PRJNA13758"

ref_dir=${ref}/${WS}

gtf="c_elegans.PRJNA13758."${WS}".canonical_geneset.gtf"
fa="c_elegans.PRJNA13758."${WS}".genomic.fa"

gtf2="c_elegans.PRJNA13758."${WS}".canonical_geneset_GFP.gtf"
fa2="c_elegans.PRJNA13758."${WS}".genomic_GFP.fa"



echo "-------   Download references   -------"
echo
if [[ ! -f $ref_dir/$gtf ]]; then
	wget -P $ref_dir $base_url/$gtf.gz
	gunzip $ref_dir/$gtf.gz
fi

if [[ ! -f $ref_dir/$fa ]]; then
	wget -P $ref_dir $base_url/$fa.gz
	gunzip $ref_dir/$fa.gz
fi

echo "    > done."
echo




echo "-------   Add transgene to annotation   -------"

if [[ ! -f $ref_dir/$gtf2 ]]; then
	cp $ref_dir/$gtf $ref_dir/$gtf2

  cat $ref/"nsIs198/nsIs198.gtf" \
    >> $ref_dir/$gtf2
fi


if [[ ! -f $ref_dir/$fa2 ]]; then
	cp $ref_dir/$fa $ref_dir/$fa2
  cat $ref/"nsIs198/nsIs198.fa"  >> $ref_dir/$fa2
fi

echo "    > done."
echo




# Filter GTF
module load CellRanger/8.0.1

echo
echo "-------   Filter GTF   -------"
echo

cellranger mkgtf \
	$ref_dir/$gtf2 \
	$ref_dir/"c_el_GFP_filt.gtf" \
	--attribute=gene_biotype:protein_coding

echo
echo "-------   Make index   -------"
echo

cellranger mkref \
	--genome=c_el_${WS}_GFP_genome \
	--fasta=$ref_dir/$fa2 \
	--genes=$ref_dir/"c_el_GFP_filt.gtf"

mv c_el_${WS}_GFP_genome intermediates/index_$WS

echo
echo "-------   Reference created.   -------"
echo




echo
echo "-------   Done.   -------"




