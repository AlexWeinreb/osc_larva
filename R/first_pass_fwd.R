# first pass, annotate each stage

# Inits ----

library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

dir_save <- "intermediates/2502/250305_per_condition"



marker_files <- list.files(dir_save)

i <- 0




# Each stage ----
i <- i+1
cond_here <- marker_files[[i]]

cond_here


markers_sub <- qs::qread( file.path(dir_save, cond_here) )


unique(markers_sub$cluster)




cl <- 0

cl <- cl +1
markers_sub |>
  filter(cluster == cl) |>
  filter(p_val_adj < 0.1,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0) |>
  pull(gene_name) |>
  clipr::write_clip()
cl

markers_sub |>
  filter(cluster == cl) |>
  filter(p_val_adj < 0.1,
         pct.1 > 0.1,
         pct.1 > pct.2,
         avg_log2FC > 0) |>
  pull(gene_name) |>
  clipr::write_clip()

markers_sub |>
  filter(cluster == cl) |>
  filter(pct.1 > pct.2) |>
  pull(gene_name) |>
  clipr::write_clip()





