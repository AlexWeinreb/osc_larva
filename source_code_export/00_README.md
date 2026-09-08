
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
+++ Fig. 2A

+++ Fig. 3B, D

+++ Fig. EV 1A, B
+++ Fig. EV 2A




### 15 Step 1a: impute, preprocess for fit


In `R/step_1a_preproc_ct.R`, called with `src/step_1a_preproc_ct.sh`:
  * inputs: Seurat object `250605_assembled/250606_seu_all_herma.qs.qs` from assembled
* impute or not, run SCT, run PCA; save plots and objects
* outputs in `250609_step1` for each cell type


note older versions: `250330_step1`: only kept cell types with cells from L2 and L4. `250502`: process all cell types with > 20 cells.



### 16 Step 2a: GAM fit and smooth curve processing

Compute pseudotime (from smoothed PCA step 1) and fit GAM (using unsmoothed counts). We fit GAM twice in a row: once without centering, that we use for representations and timings.

We use the peak of the uncentered to run a second fit, on a pre-centered curve. From this, we keep
* coefficients of the GAM
* descriptors of the smooth curve (amplitude, auc, ...)
* similarity of smooth curve to thin peak (DTW distance)

These are used for clustering.





`step_2a_gam.R` is to be run on cluster as dsq jobarray, called from dSQ. Contents:
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

Notes older versions:
* previously used pseudotimeDE at this step, along with filtering on curve shape as step 3. No longer useful: most/all genes appear DE with pseudotime, replace with simple GAM and curve shape filtering. Keeping state of repo at that point in branch `pseudotimede`.
* used binomial fit in some versions, later used NB-GAM on raw counts without offset, and Gaussian GAM on SCT data
* bootstraps on dtw in previous version: CI wasn't obviously a better predictor than the dtw distance itself








### Step 1b: Velocyto

Overview:
  * rerun CellRanger if needed,
* replace filtered matrix,
* Velocyto to quantify unspliced reads,
* reorganize by cell type (based on step 1a),
* plot each cell type


#### CellRanger

Reran `cr_count.sh` with saving bam (in scratch dir: `250331_align`).

Notes the bams are often 10-20 GB, not easily stored on long term.


#### 17 Replace filtered matrix

Problem: velocyto will always read `filtered_feature_bc_matrix`. But here I use a manual reannotation keeping more cells.

In `R/save_filtered_matrices_for_velocyto.R`, for each sample we rename the CellRanger `filtered_feature_bc_matrix` and replace it with an export of the Seurat object.
* inputs: `scratch/250331_align` bams and `250328_assembled/250329_seu_all_herma.qs`
* for each sample,
* rename `filtered_...` to `cr_filtered_feature_bc_matrix`
* take the count matrix from `assembled` (subset sample), rename cell bc if needed
* save this subset of "assembled" in `filtered_feature_bc_matrix`


#### 18 Velocyto quantification and plot

In `src/velocyto_sample.sh`
* input: `250331_align/{sample}/outs/` for each sample (uses bam and filtered matrix)
* run `velocyto run10x`
* output in `250331_align/{sample}/velocyto/sample.loom`

Run as dsq array; joblist file in `joblists/velocyto_samples.dsq.txt`

Contents (one row per sample):
  ```
bash ./src/velocyto_sample.sh "200730_batch1_CHB3840b"
bash ./src/velocyto_sample.sh "201013_batch2_CHB3840b_CEG_fqs"
...
```

dsq prepared with:
  ```
ml dSQ; dsq --job-file joblists/velocyto_samples.dsq.txt  --cpus-per-task 6 --mem 20G --time 1:50:00 --partition day
```



#### Save copy of loom files


Copy files out of scratch, work directly from them later.
```
cp -v /vast/palmer/scratch/hammarlund/aw853/250331_align/*/velocyto/*.loom intermediates/2502/250409_loom/
  ```


#### 19 Reorganize by cell type

Consistently with other approach, we split by cell type and process each cell type separately.

In `R/velocyto_load_loom.R`, called from `src/runR_velocyto_load_loom.sh`
* inputs: `250409_loom/{sample}.loom`, `250605_assembled/250606_seu_all_herma.qs`
* Process:
  * read all loom files using velocyto.R
* combine into big "spliced" and "unspliced" matrices
* output: matrices in `250825_velocyto/emat_tot.qs` and `nmat_tot.qs`


#### 20 Plot velocyto

For each cell type, run `velocyto_cell_type.R`. Inputs:
  * from step 1, the seu_unsmoothed object to reuse its PCA and average phase precomputed
* from dir_velocyto, the total matrices emat and nmat, to compute velocity.
Outputs:
  * plots in pdf and png of PCA with/without colors and arrows
* preprocessed objects to replot (qs format).


Using joblist:
  
  
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

The previous velocyto plots may not look good. But since we save the preprocessed object, we can easily replot cells of interest with custom parameters.

In `velocyto_celltypes_replot.R`, code copied from `velocyto_cell_type.R` but for more interactive use.

Save in `presentatons/figures/250825/velocyto`.




###### older

Note: older version
instead, interactively ran `R/velocyto.R`
* for each cell type annotated in "assembled", subset the corresponding cells, create AnnData object
* output: `250409_anndata/{sample}.h5ad`


At the end of `R/velocyto.R`, additional code for velocity estimate and plotting with velocyto.R (not used in practice).


Note: previous attempts to use VeloCycle, in folder `ipynb/`. While it seems to work, it's not answering the questions I have here.


then step 2b: scVelo



With `step_2b_scvelo_cell_type.R` run with dSQ:
* input: `250409_anndata/{sample}.h5ad`
* minimal filtering, recover_dynamics (not actually used), velocity with *stochastic* model, extract genes by fit_likelihood
* outputs:
  * `250501_scvelo/{sample}_velocity.png`
  * `250501_scvelo/{sample}_scvelo_fit.qs` table of genes by fit_likelihood
  * `250501_scvelo/{sample}_adata.pkl` with processed object (not used)

Using joblist:


```r
paste(
"module load R; Rscript R/step2b_scvelo_cell_type.R",
"--dir_in_anndata 'intermediates/2502/250409_anndata'",
"--dir_out_scvelo 'intermediates/2502/250501_scvelo'",
"--i", seq_along(list.files(params$dir_in_anndata, pattern = "\\.h5ad$") |> str_subset("scvel", negate = TRUE))
) |>
  writeLines("joblists/step_2b_scvelo_ct.dsq.txt")
```

Job run with
```
ml dSQ; dsq --job-file joblists/step_2b_scvelo_ct.dsq.txt  --cpus-per-task 1 --mem 20G --time 00:40:00 --partition day; ml unload dSQ
```

Tests and manual version in `test_scVelo.R` (not used).









### 21-22 Clustering

In `R/cluster_cellgenes.R`, called from `runR_cluster_hclust.sh`:
* load descriptors from step 2a `250624_step2`
* transformations (exp(-a*x)), normalize (box-cox), scale; cluster with fastclust::hclust
* save

In `R/hclust_results.R`, load this clustering, cutree and cluster identification.


All results saved in `250624_cluster`.


In `explore_step3_timeseries_distances.R`, temporary explorations, to delete later.

In `R/step_3_process_celltypes.R`, temporary explorations (manually annotate some genes to compare to gene clustering).




### 23 Step 3: process cell types, curve shape, heatmaps


In `step3_heatmap_from_clust.R`, look at each cell type's pulsatile genes. Categorize cell types as oscillatory or not based on entropy of peaks.



### 24-25 Step 4: Analysis


In `R/step_4_analysis.R` look at genes.

In `R/CelEst.R` looks at TFs.










