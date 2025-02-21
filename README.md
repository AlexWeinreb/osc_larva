
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












