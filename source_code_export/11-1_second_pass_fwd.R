
# second pass, annotate each tissue in each stage



library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

dir_save <- "intermediates/2502/250306_per_condition_tissue"



marker_files <- list.files(dir_save, pattern = "^250306_marks_.+qs")

i <- 0




# Each file ----
i <- i+1
file_here <- marker_files[[i]]

file_here


markers_sub <- qs::qread( file.path(dir_save, file_here) )


unique(markers_sub$cluster)




cl <- 0

cl <- cl +1
markers_sub |>
  filter(cluster == cl) |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0) |>
  pull(gene_name) |>
  clipr::write_clip()
cl

markers_sub |>
  filter(cluster == cl) |>
  filter(p_val_adj < 0.1,
         pct.1 > 0.1,
         avg_log2FC > 0) |>
  pull(gene_name) |>
  clipr::write_clip()

markers_sub |>
  filter(cluster == cl) |>
  filter(pct.1 > 0.2) |>
  pull(gene_name) |>
  clipr::write_clip()






