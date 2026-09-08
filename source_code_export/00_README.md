
# Data

Samples from 2 series of experiments.

From "glia_10x":
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

From grl-18: 
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

Recorded in `data/raw_paths.txt`, with the absolute path to raw and `data/raw_sample_ids.txt`
as second file (keep separate files to easily `mapfile`).

Description in Table EV1



# Reads/cells pre-processing

## 01/ alignment


Create transcriptome file with `01_cellranger_prep_annotation.sh`.

Save index in `intermediates`

Everything with WS295.



## 02/ count

In `02_cellranger_align_count.sh`, load lists of samples from `data/raw_xx.txt` and run on each sample. 
Save CellRanger alignments in `data/250225_aligned/`.

Note (see below), we will need the BAM files for velocyto.



## 03-04/ emptydrops



Run emptydrops as a Dead Simple Queue (dsq) job using a computing cluster, then load it and filter the actual cells.

File `03_emptydrops.dsq.txt` created with:
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


Running `04_emptydrops.R` (called from a SLURM wrapper on computing cluster).



* Inputs:
 * raw matrix from "data/250225_aligned"
 * precomputed emptydrops from "intermediates/2502/emptydrops"
Then load the precomputed emptydrops object and subset the matrix
* Output:
 * filtered matrix in "intermediates/2502/empytdrops_filt/{sample}filt_feature_bc_matrix"

(note the typo at "empytdrops", note the lack of separator between sample name and "filt")



## 05/ SoupX


With `05_soupX.R`, for each sample, `for` loop on cluster
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





## 06/ DoubletFinder


With `06_doubletFinder.R` (called from wrapper)
* input: `soupx_feature_bc_matrix` for each sample
  * basic filtering, pca, clustering, DoubletFinder (discarded, we only keep the annotation)
* output: Seurat object from `soupx_feature_bc_matrix` with an additional metadata column "doubletFinder",
        saved in `2502/250304_doubletFinder/250304_{sample}_dblts.qs`

Note the log file can be filtered with `grep -e "-----  "`.

Important: examining plots, it appears the pK parameter selection may be incorrect for several
samples: maximum at either first of last value, instead of a peak.

Folloing visual examination, the results of this doublet detection were not explicitly used in subsequent steps.



## 07/ QC per sample


With `07_qc_per_sample.R`, run interactively, for each sample:
* inputs: `2502/250304_doubletFinder/250304_{sample}_dblts.qs`
  * for each sample, remove cells with high mito or low count
  * cluster (SCT), basic annotation of clusters, eliminate bad clusters
* output: `250304_filt_ds/{sample}.qs` for each sample





# Annotation

The samples went through several rounds of annotation.

First, the samples were gathered by condition, all samples from the same condition merged. 
The Seurat object for each condition was clustered, the clusters annotated (first pass). 
Then the objects were split by tissue (getting a separate object for each tissue and condition), 
and further clustered and annotated at the cell-type level (second pass). Then a third 
pass specifically looked at the clusters that still had subclusters (notably neurons) to
ensure final annotation only contained atomic clusters.

Each round of annotation uses "forward" and "reverse" method. Forward: cluster markers obtained 
from data (e.g. Seurat's FindMarkers), pasted into CengenApp's heatmap to see if these markers
are enriched for an annotated tissue. Reverse: Starting with a precompiled list of tissue-
specific markers (e.g. myo-3 for muscles), look at their expression in data. The annotation 
is then manually compiled from combining this information.


## 08/ Assemble per condition

4 groups:
* group 1, L2 herma
* group 2, L4 herma
* group 3, L4 male/feminized
* group 4, adult herma

note we don't group by promoter.

| sample_name                    | sex           | promoter | stage | group |
| ------------------------------ | ------------- | -------- | ----- | ----- |
| 200730_batch1_CHB3840b         | hermaphrodite | mir-228  | L2    | 1     |
| 201013_batch2_CHB3840b_CEG_fqs | hermaphrodite | mir-228  | L2    | 1     |
| 201013_batch2_CHB3841_CEG_fqs  | male          | mir-228  | L4    | 3     |
| 210413_batch3_CHB3840b         | hermaphrodite | mir-228  | L4    | 2     |
| 210413_batch3_CHB3841          | male          | mir-228  | L4    | 3     |
| 210420_batch4_CHB3841          | male          | mir-228  | L4    | 3     |
| 210420_batch4_hmn              | other         | mir-228  | L4    | 3     |
| 210427_batch5_CHB3840b         | hermaphrodite | mir-228  | L4    | 2     |
| 210427_batch5_CHB3841          | male          | mir-228  | L4    | 3     |
| 220210_OH17400                 | hermaphrodite | mam-5    | L4    | 2     |
| 230414                         | hermaphrodite | grl-18   | L4    | 2     |
| 230421_AM                      | hermaphrodite | grl-18   | L4    | 2     |
| 230421_PM                      | hermaphrodite | grl-18   | L4    | 2     |
| 230505_AM                      | hermaphrodite | grl-18   | L4    | 2     |
| 230505_PM                      | hermaphrodite | grl-18   | L4    | 2     |
| 231019GRL                      | hermaphrodite | grl-18   | adult | 4     |
| 240111                         | hermaphrodite | grl-18   | adult | 4     |
  

+++ This data is also provided in Table EV1


For each group, load the corresponding files, assemble them (no integration), and annotate (next step).




With `08_assemble_per_condition.R`, interactively:
* inputs: `250304_filt_ds/{sample}.qs` grouped based on `data/samples_table.tsv`
  * Seurat merge per condition
  * SCT, clustering, find markers
* outputs:
  * `250305_per_condition/250305_seu_merged_{group}.qs`
  * `250305_per_condition/250305_marks_merged_{group}.qs`






## 09/ First pass (annotate tissues per condition)

Two scripts (forward and reverse, see above).

On PC, run `09-1_first_pass_annotation_fwd.R`:
* input: `250305_per_condition/250305_marks_merged_{group}.qs`
  * go through clusters, compare cengenapp


On cluster, run `09-2_first_pass_annotation_rev.R`:
* input: `250305_per_condition/250305_seu_merged_{group}.qs`
  * look at known tissue markers, compare with fwd results

Save annotation in csv files as `250305_per_condition/250305_annot_{group}.csv`.




## 10-11/ Second pass per tissue and condition

With `10_second_pass_subset.R`, on cluster
* inputs:
  * `250305_per_condition/250305_seu_merged_{group}.qs`
  * annotation in `250305_per_condition/250305_annot_{group}.csv`
* For each stage, separate clusters by tissue, recluster, find markers
* outputs:
  * `250306_per_condition_tissue/250306_seu_{group}_{tissue}.qs`
  * `250306_per_condition_tissue/250306_marks_{group}_{tissue}.qs`


Note: discarding "unclear" clusters.



With `11-1_second_pass_fwd.R`, on PC
* input: `250306_per_condition_tissue/250306_marks_{group}_{tissue}.qs`
* forward marker selection

With `11-2_second_pass_rev.R` on cluster


Save annotations in `250306_per_condition_tissue/250306_annot_tissues.xlsx` (piling the different groups and tissues in the same file).



## 12-13/ Third pass

Separate/assemble by second-pass annotation, reannotate to ensure no subclusters. In `12_third_pass_subset.R`.

Saved in `intermediates/2502/250311_third`:
  * `250311_third_subsets.xlsx`: description of the subsets
  * `250311_(marks|seu)_{group}_{tissue}.qs`: single tissue or cell type subset

Further reprocessing those that are not an atomic cell type (single cluster), in
`13-1_third_pass_rev.R`. Note this script covers both reverse and forward approach (no separate script to look at markers).

For neurons, reprocessed separately in `13-2_third_pass_neurs.R`.

The individually processed files are finally saved in `intermediates/2502/250529_third_processed`.

The list of cell types, their putative identities, and the gene markers, are recorded in `data/cell_annotations.xlsx`.



## 14 Assemble and plot

In `14_assemble_main_seurat_object.R` assemble all cell types from L2 and L4 (herma only).

* inputs: `250313_third_processed/230318_g(1|2)_{cell_type}.qs` for each cell type and each stage
  * assemble into single big Seurat object, re-umap
  * Compute cell phases, permutation tests, add as columns in Seurat object
  * produces several plots (copied out), several intermediates (saved in `250509_assembled`)
* outputs:
  * Main Seurat object saved as `250509_assembled/250509_seu_all_herma.qs`

+++ Fig. 1B, C
+++ Fig. 2A, B, C

+++ Fig. 3A, B, D

+++ Fig. EV 1A, B
+++ Fig. EV 2A, B, C, D




### 15/ impute, preprocess for fit


In `15_impute_preprocess_cell_type.R`, (called from wrapper):
* inputs: Seurat object `250605_assembled/250606_seu_all_herma.qs` from `14_assemble...`
  * impute or not, run SCT, run PCA; save plots and objects
* outputs in `250609_step1` for each cell type







### 16/ Pseudotime, GAM fit and smooth curve processing

Compute pseudotime (from smoothed PCA step 1) and fit GAM (using unsmoothed counts).

We fit GAM twice in a row: First one "uncentered", that we use for representations and timings.
We use the peak position of the uncentered to position the peak in the middle, 
and run a second "centered" fit. From this, we keep
 * coefficients of the GAM
 * descriptors of the smooth curve (amplitude, auc, ...)
 * similarity of smooth curve to thin peak (DTW distance)

These are used for clustering.





`16_gam_smooth_curve_processing.R` (run as dSQ job `16-1_gam.dsq.txt`). Contents:
  * load intermediates from step 1 `250609_step1`
  * prefilter, run ElPiGraph and a NB-GAM (with size factors), save the models, descriptors, smooth fits
  * save intermediates in "intermediates/2502/250624_step2"

Outputs:
  * `{cell_type}_descriptors.qs` used for clustering (step 3)
  * `{cell_type}_preds.qs` (uncentered smooth) and `{cell_type}_preds_cent_clipped.qs` (centered, clipped smooth), for plotting
  * the model objects from `mgcv::gam` (shouldn't be needed)



Jobfile created with:
```r
paste("module load mlq; ml R; Rscript R/step_2a_gam_nb.R",
       "--dir_step1 'intermediates/2502/250609_step1'",
      "--out_dir 'intermediates/2502/250624_step2'",
      "--i", seq_along(list.files('intermediates/2502/250609_step1', pattern = "_seu\\.qs$")),
      "--prop_thres 0.05 --cnt_thres 20") |>
  writeLines("joblists/step_2a_gam.dsq.txt")
```


Job prepared with:
```
ml dSQ; dsq --job-file joblists/step_2a_gam.dsq.txt  --cpus-per-task 1 --mem 15G --time 00:20:00 --partition day; ml unload dSQ
```

For plotting: `16-x1_plot_elpigraph.R`
+++ Fig. EV3B






# Velocyto

Overview:
  * rerun CellRanger if needed,
  * replace filtered matrix,
  * Velocyto to quantify unspliced reads,
  * reorganize by cell type (based on step 1a),
  * plot each cell type


## CellRanger



Note: the bams are often 10-20 GB, not easily stored on long term. If needed, 
rerun `02_cellranger_align_count.sh` with saving bam (in scratch dir: `250331_align`).


## 17/ Prepare filtered matrix as Velocyto input

Problem: velocyto will always read `filtered_feature_bc_matrix`, usually provided by CellRanger. 
But here I use a manual EmptyDrops reannotation keeping more cells.

In `17_save_filtered_matrices_for_velocyto.R`, for each sample we rename the 
CellRanger `filtered_feature_bc_matrix` and replace it with an export of the Seurat object.
* inputs: `scratch/250331_align` bams and `250328_assembled/250329_seu_all_herma.qs`
* for each sample,
 * rename `filtered_...` to `cr_filtered_feature_bc_matrix`
 * take the count matrix from `assembled` (subset sample), rename cell bc if needed
 * save this subset of "assembled" in `filtered_feature_bc_matrix`



## 18/ Velocyto quantification and plot

In `18_velocyto_sample.sh` (dsq array, cf 18-1)
* input: `250331_align/{sample}/outs/` for each sample (uses bam and filtered matrix)
* run `velocyto run10x`
* output in `250331_align/{sample}/velocyto/sample.loom`


dsq prepared with:
  ```
ml dSQ; dsq --job-file joblists/velocyto_samples.dsq.txt  --cpus-per-task 6 --mem 20G --time 1:50:00 --partition day
```



## Save copy of loom files


Copy files out of scratch, work directly from them later.

```
cp -v /vast/palmer/scratch/hammarlund/aw853/250331_align/*/velocyto/*.loom intermediates/2502/250409_loom/
```


## 19/ Reorganize by cell type

Consistently with other approach, we split by cell type and process each cell type separately.

In `19_velocyto_load_loom.R` (called from wrapper)
* inputs: `250409_loom/{sample}.loom`, `250605_assembled/250606_seu_all_herma.qs`
* Process:
 * read all loom files using velocyto.R
 * combine into big "spliced" and "unspliced" matrices
* output: matrices in `250825_velocyto/emat_tot.qs` and `nmat_tot.qs`



## 20/ Plot velocyto

For each cell type, run `20_velocyto_cell_type.R`. Inputs:
 * from 15_impute..., the seu_unsmoothed object to reuse its PCA and average phase precomputed
 * from dir_velocyto, the total matrices emat and nmat, to compute velocity.
Outputs:
 * plots in pdf and png of PCA with/without colors and arrows
 * preprocessed objects to replot (qs format).


Using joblist (cf 20-1):

```r
paste(
  "module load R; Rscript R/velocyto_cell_type.R",
  "--dir_step1 'intermediates/2502/250609_step1'",
  "--dir_velocyto 'intermediates/2502/250825_velocyto'",
  "--i", seq_along(list.files('intermediates/2502/250609_step1', pattern = "_seu_unsmoothed\\.qs$"))
) |>
  writeLines("joblists/velocyto_cell_type.dsq.txt")
```

Job run with
```
ml dSQ; dsq --job-file joblists/velocyto_cell_type.dsq.txt  --cpus-per-task 1 --mem 35G --time 00:40:00 --partition day; ml unload dSQ
```



#### 20.2 replot nicely

The previous velocyto plots may not look good. But since we save the preprocessed object, 
we can easily replot cells of interest with custom parameters.

In `20-2_velocyto_celltypes_replot.R`, code copied from `velocyto_cell_type.R` but for more interactive use.

+++ Fig. 2D
+++ Fig. 3C



# Calling pulsatile genes

## 21/ Clustering of gene expression profiles

In `21_hclust_identify_pulsatile_genes.R` (called from wrapper):
  * load descriptors from `16_gam_smooth_curve_processing`
  * transformations (exp(-a*x)), normalize (box-cox), scale; cluster with fastclust::hclust
  * save


## 22/ Identification of pulsatile genes

In `22_analyze_plot_hclust_results.R`, load 21's clustering, cutree and cluster identification.
All results saved in `250624_cluster`.

+++ Fig. 4C, D, E

+++ Fig. EV3C, D

+++ Table EV3
+++ Table EV4





## 23/ process cell types, heatmaps


In `23_analyze_oscillatory_cell_types.R`, look at each cell type's pulsatile genes. 
Categorize cell types as oscillatory or not based on entropy of peaks.


+++ Fig. 5A, B, C

+++ Table EV5


### 24/ Analysis of pulsatile genes


In `24_analyze_pulsatile_genes.R` look at pulsatile genes in oscillatory cell types.

See also `x01_panther_build_dict.R` for the PANTHER.db dictionary that is used in that script.

+++ Fig. 4A, F
+++ Fig. 5D, E, F, G

+++ Fig. EV3A
+++ Fig. EV4B, C
+++ Fig. EV5A, B

+++ Table EV6 (manually compiled, used as input)



# Additional analyzes


## 25/ Transcription factors and timing

In `25_CelEst_transcription_factors.R` looks at TFs.

Load the CelEst GRN 

+++ Fig. 6A, B
+++ Fig. EV6A, B
+++ Fig. EV7
+++ Table EV7


## 26/ Peak width in the Meeuse dataset


`26_meeuse_peak_width.R` loads data provided from Meeuse et al (2020) to analyze peak width.
 
+++ Fig. EV4A


## 27/ Subsample perplexity

`27_subsample_perplexity.R` 

+++ Fig. EV5C





# Figures vs scripts


| Script | Figures |
|:-------|:--------|
| 14     | Fig 1B, 1C, 2A, 2B, 2C, 3A, 3B, 3D; EV 1A, 1B, 2A, 2B, 2C, 2D |
| 16     | EV 3B |
| 20.2   | Fig 2D, 3C |
| 22     | Fig 4D, 4E; EV 3C, 3D |
| 23     | Fig 5A, 5B, 5C |
| 24     | Fig 4A, 4F, 5D, 5E, 5F, 5G; EV 3A, 4B, 4C, 5A, 5B |
| 25     | Fig 6A, 6B; EV 6A, 6B, 7 |
| 26     | EV 4A |
| 27     | EV 5C |








