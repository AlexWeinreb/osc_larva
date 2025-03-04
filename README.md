
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
























