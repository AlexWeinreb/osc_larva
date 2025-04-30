
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

Everything with WS295.

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



## Assemble per condition

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

This table is also saved in `data/samples_table.tsv`


For each group, load the corresponding files, assemble them (no integration), and annotate (next step).




With `assemble_per_condition.R`, interactive on cluster:
* inputs: `250304_filt_ds/{sample}.qs` grouped based on `data/samples_table.tsv`
* Seurat merge per condition
* SCT, clustering, find markers
* outputs:
    * `250305_per_condition/250305_seu_merged_{group}.qs`
    * `250305_per_condition/250305_marks_merged_{group}.qs`


(note: script adapted from `larval_devt/bseu_assemble.R`)




## First pass per condition

Two scripts

On PC, run `first_pass_fwd.R`:
* input: `250305_per_condition/250305_marks_merged_{group}.qs`
* go through clusters, compare cengenapp


On cluster, run `first_pass_rev.R`:
* input: `250305_per_condition/250305_seu_merged_{group}.qs`
* look at known tissue markers, compare with fwd results

Save annotation in csv files as `250305_per_condition/250305_annot_{group}.csv`.

Keep notes and UMAPs in `presentations/250305_first_pass_per_condition.pptx`.


(note: scripts adapted from `larval_devt/bseu_first_pass_forward` and `bseu_first_pass_reverse`)




### Second pass per tissue and condition

With `second_pass_subset.R`, on cluster
* inputs:
    * `250305_per_condition/250305_seu_merged_{group}.qs`
    * annotation in `250305_per_condition/250305_annot_{group}.csv`
* For each stage, separate clusters by tissue, recluster, find markers
* outputs:
    * `250306_per_condition_tissue/250306_seu_{group}_{tissue}.qs`
    * `250306_per_condition_tissue/250306_marks_{group}_{tissue}.qs`


Note: discarding "unclear" clusters.



With `second_pass_fwd.R`, on PC
* input: `250306_per_condition_tissue/250306_marks_{group}_{tissue}.qs`
* forward marker selection

With `second_pass_rev.R` on cluster


Save annotations in `250306_per_condition_tissue/250306_annot_tissues.xlsx` (piling the different groups and tissues in the same file).



### Third pass

Separate/assemble by second-pass annotation, reannotate to ensure no subclusters. In `R/third_pass_subset.R`.

Saved in `intermediates/2502/250311_third`:
* `250311_third_subsets.xlsx`: description of the subsets
* `250311_(marks|seu)_{group}_{tissue}.qs`: single tissue or cell type subset

Further reprocessing those that are not an atomic cell type (single cluster), in `R/third_pass_rev.R`. Note it includes the forward approach (no separate script to look at markers).

For neurons, reprocessed separately in `third_pass_neurs.R`.

The individually processed files are finally saved in `intermediates/2502/250313_third_processed`



### Assemble and plot

In `assemble_osc.R` (interactive on cluster) assemble all cell types from L2 and L4 (herma only).

* inputs: `250313_third_processed/230318_g(1|2)_{cell_type}.qs` for each cell type and each stage
* assemble into single big Seurat object, re-umap
* Compute cell phases, permutation tests, add as columns in Seurat object
  * produces several plots (copied out), several intermediates (saved in `250328_assembled`)
* outputs:
  * Main Seurat object saved as `250328_assembled/250329_seu_all_herma.qs`
  * Re-exporting individual cell types in `250328_assembled/250330_cell_types/{cell_type.qs}`



# Velocyto (test)

Reran `cr_count.sh` with saving bam (in scratch dir).

Problem: velocyto will always read `filtered_feature_bc_matrix`. But my annotation has more cells. In `save_filtered_matrices_for_velocyto`, for each sample we rename the CellRanger `filtered_feature_bc_matrix` and replace it with an export of the Seurat object.

Then we run velocyto on each samples (as dsq array). Joblist file in `joblists/velocyto.dsq.txt`

Contents (one row per sample):
```
bash ./src/velocyto_sample.sh "200730_batch1_CHB3840b"
bash ./src/velocyto_sample.sh "201013_batch2_CHB3840b_CEG_fqs"
```

dsq prepared with:
```
ml dSQ; dsq --job-file joblists/velocyto_samples.dsq.txt  --cpus-per-task 6 --mem 20G --time 1:50:00 --partition day
```





### Step 1: impute, preprocess


In `R/step_1_preproc_ct.R`, called with `src/step_1_preproc_ct.sh`:
* inputs: Seurat object `250328_assembled/250329_seu_all_herma.qs` from assembled
* impute, run SCT, run PCA
* outputs in `250330_step1` for each cell type




### Step 2


`step_2_gam_binom.R` is to be run on cluster as dsq jobarray, called from dSQ. Contents:
* load intermediates from step 1, prefilter, run ElPiGraph and a GAM with binomial family, save the model along with curve amplitude and dev explained
* save intermediates in "intermediates/2502/250330_step2"




Jobfile created with:
```r
paste("module load R; Rscript R/step_2_gam_binom.R",
       "--batch_rmed_dir 'intermediates/2502/250330_step1'",
      "--out_dir 'intermediates/2502/250330_step2'",
      "--i", seq_along(list.files('intermediates/2502/250330_step1', pattern = "_seu\\.qs$")),
      "--model 'auto' --prop_thres 0.05 --cnt_thres 20") |>
  writeLines("joblists/step_2_gam_binom.dsq.txt")
```


Job prepared with:
```
ml dSQ; dsq --job-file joblists/step_2_gam_binom.dsq.txt  --cpus-per-task 1 --mem 5G --time 00:10:00 --partition day
```

Note: previously used pseudotimeDE at this step, along with filtering on curve shape as step 3. No longer useful: most/all genes appear DE with pseudotime, replace with simple GAM and curve shape filtering. Keeping state of repo at that point in branch `pseudotimede`.



Testing with bootstraps:
```
paste("module load R; Rscript R/step_2_gene_expr_nb_bootstrap_dtw.R",
       "--batch_rmed_dir 'intermediates/2502/250330_step1'",
      "--out_dir 'intermediates/2502/250424_step2_boot_nb'",
      "--i", seq_along(list.files('intermediates/2502/250330_step1', pattern = "_seu\\.qs$")),
      "--prop_thres 0.05 --cnt_thres 20") |>
  writeLines("joblists/step_2_gam_boot_nb.dsq.txt")
```

Run with:
```
ml dSQ; dsq --job-file joblists/step_2_gam_boot_nb.dsq.txt  --cpus-per-task 1 --mem 15G --time 20:00:00 --partition day; ml unload dSQ
```

Fail for ILso, pha epith, hyp. Rerun:
```
ml dSQ; dsq --job-file joblists/step_2_gam_boot_nb.dsq.txt  --cpus-per-task 1 --mem 15G --time 3-20:00:00 --partition week; ml unload dSQ
```


Third alternative: fit with NB,
```
dsq --job-file joblists/step_2_gam_nb.dsq.txt --cpus-per-task 1 --mem 5G --time 00:10:00 --partition day
```




### Step 3: process cell types, curve shape, heatmaps

On cluster, interactively, run `src/runR_step_3_process_celltypes`:
* inputs:
  * result of step 2 `2502/250325_step2_binom/{celltype}_res_gam.qs` for each cell type
  * also use step1 to examine individual genes
* determine what genes are peaky based on curve shape
* plot heatmaps per cell type, determine if heatmap diagonal
* save results in `2502/250330_step3_genes_by_celltype/`













