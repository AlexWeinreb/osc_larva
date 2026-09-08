# Inits ----


library(tidyverse)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

dir_prefilt <- "intermediates/2502/250304_filt_ds"
dir_out <- "intermediates/2502/250305_per_condition"



sample_table <- read_tsv("data/samples_table.tsv",
                         col_types = "cffff")

i <- 0


# condition ----
i <- i+1

cond_here <- unique(sample_table$group)[[i]]

cond_here

samples_here <- sample_table |>
  filter(group == cond_here) |>
  pull(sample_name)

datasets <- file.path(dir_prefilt,
                      samples_here |>
                        paste0(".qs") ) |>
  set_names(samples_here) |>
  map(qs::qread) |>
  map(~ {DefaultAssay(.x) <- "RNA"; .x}) |>
  imap(~ {.x[["orig.ident"]] <- .y; .x})

datasets

if(length(datasets) > 1){
  seu <- merge(datasets[[1]],
               datasets[-1])
} else{
  seu <- datasets[[1]]
}




#~ Clustering ----
seu <- SCTransform(seu)


seu <- RunPCA(seu, npcs = 200, verbose = FALSE)

ElbowPlot(seu, ndims = 200)

ndims <- 30

seu <- FindNeighbors(seu,
                     dims = 1:ndims,
                     verbose = FALSE)



seu <- RunUMAP(seu,
               dims = 1:ndims)

DimPlot(seu,
        reduction = "umap",
        group.by = "orig.ident",
        label = TRUE,
        pt.size = 2,
        alpha = .08) +
  NoLegend()

seu <- FindClusters(seu,
                    resolution = .5)
DimPlot(seu,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend()



seu <- PrepSCTFindMarkers(seu)



all_markers <- FindAllMarkers(seu,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)


file.exists(file.path( dir_out, paste0("250305_seu_merged_",cond_here,".qs") ))
file.exists(file.path( dir_out, paste0("250305_marks_merged_",cond_here,".qs") ))


qs::qsave(seu,
          file.path( dir_out, paste0("250305_seu_merged_",cond_here,".qs") ))
qs::qsave(all_markers,
          file.path( dir_out, paste0("250305_marks_merged_",cond_here,".qs") ))





