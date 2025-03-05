
# Data

Samples from "glia_10x", raw:
```
ycga_work/glia_10x_old/raw/*
├── 2020-07-30_batch1
│   └── CHB3840b
├── 2020-10-13_batch2
│   ├── CHB3840b_CEG_fqs
│   └── CHB3841_CEG_fqs
├── 2021-04-13_batch3
│   ├── CHB3840b
│   └── CHB3841
├── 2021-04-20_batch4
│   ├── CHB3841
│   └── hmn
├── 2021-04-27_batch5
│   ├── CHB3840b
│   └── CHB3841
└── 2022-02-10_OH17400
```

Samples in grl-18 sorts, raw: 
```
project/10x_grl18/data/raw/*
├── 230414
├── 230421_AM
├── 230421_PM
├── 230505_AM
├── 230505_PM
├── 231019GRL
└── 240111
```

Recorded in `data/raw_paths.txt`, with the absolute path to raw and `data/raw_sample_ids.txt` as second file (keep separate files to easily `mapfile`).


# alignment

Create transcriptome file with `src/cr_prep_annot.sh`.

Save index in `intermediates`


# count

In `src/cr_count.sh`, load lists of samples from `data/raw_xx.txt` and run on each sample. Save CellRanger alignments in `data/250225_aligned/` (note, the script just saved them in the current directory; results were transfered after running the script).


# emptydrops

Following previous scripts in `larval_devt` repo, first run emptydrops as dsq, then load it and filter the actual cells.

File `joblists/emptydrops.dsq.txt` created with:
```
writeClipboard(
  
  paste0(
    r'(module load R; R -e 'library(DropletUtils) |> suppressPackageStartupMessages(); data_dir <- "data/250225_aligned"; samples <- list.files(data_dir); sample <- samples[[)',
    1:17,
    r'(]]; dat <- read10xCounts(file.path(data_dir,sample,"outs","raw_feature_bc_matrix")); e.out <- emptyDrops(dat, lower = 100, niter = 1e6); qs::qsave(e.out, file.path("intermediates/2502/emptydrops", paste0(sample,"_e_out.qs"))) ')'
  ))
```

Then the dsq file is created interactively with:
```
ml dSQ; dsq --job-file joblists/emptydrops.dsq.txt  --cpus-per-task 1 --mem 7G --time 4:40:00 --partition day
```

Note sample #5 (job 04) failed after 4h40', rerun with 23h.

Second phase: running `src/runR_emptydrops.sh` which calls `R/emptydrops.R`. This is modified from `larval_devt` but removing the iterative cleaning (not required for this dataset).

Inputs:
* raw matrix from "data/250225_aligned"
* precomputed emptydrops from "intermediates/2502/emptydrops"
Then load the precomputed emptydrops object and subset the matrix
Output:
* filtered matrix in "intermediates/2502/empytdrops_filt/{sample}filt_feature_bc_matrix"

(note the typo at "empytdrops", note the lack of separator between sample name and "filt")



## SoupX


With `soupX.R`, for each sample, `for` loop on cluster
Can be run interactively for individual samples, for log run `src/soupX.sh`
* inputs:
  * the filtered matrix from emptyDrops
  * the raw droplets matrix aligned
* basic Seurat clustering, run SoupX using list of tissue-specific marker genes
* output: `soupx_feature_bc_matrix` for each sample


Note: when rerunning soupX, we can't overwrite existing files. Clean up with:
```
cd data/cellranger/
rm L*/soupx_feature_bc_matrix/*.tsv.gz
rm L*/soupx_feature_bc_matrix/*.mtx.gz
rmdir L*/soupx_fea*
```





### DoubletFinder


With `doubletFinder.R` called from `src/runR_doubletFinder.sh`
* input: `soupx_feature_bc_matrix` for each sample
* basic filtering, pca, clustering, DoubletFinder (discarded, we only keep the annotation)
* output: Seurat object from `soupx_feature_bc_matrix` with an additional metadata column "doubletFinder", saved in `2502/250304_doubletFinder/250304_{sample}_dblts.qs`

Note the log file can be filtered with `grep -e "-----  "`.

Important: examining plots, it appears the pK parameter selection may be incorrect for several samples: maximum at either first of last value, instead of a peak.



### QC per sample


With `qc_per_sample.R`, interactive on cluster, for each sample:
* inputs: `2502/250304_doubletFinder/250304_{sample}_dblts.qs`
* for each sample, remove cells with high mito or low count
* cluster (SCT), basic annotation of clusters, eliminate bad clusters
* output: `250304_filt_ds/{sample}.qs` for each sample


(note, script adapted from `larval_devt/seurat_qc_fiilter_clusters`)




















