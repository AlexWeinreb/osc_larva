
# Inits ----
Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")
Sys.unsetenv("RETICULATE_PYTHON")

library(reticulate)

use_virtualenv("scvelo", required = TRUE)

scv <- import("scvelo")
sc <- import("scanpy")
plt <- import("matplotlib.pyplot")



library(tidyverse)
library(wbData)


gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = "a",
                  gene_id = "nsIs198",
                  symbol = "GFP",
                  sequence = "GFP",
                  status = "Live",
                  biotype = "protein_coding_gene",
                  name = "GFP"
  )


library(getopt)


if(! interactive()){
  spec <- matrix(c(
    'dir_in_anndata', 'f', 1, 'character',
    'dir_out_scvelo', 'o', 1, 'character',
    'i', 'i', 1, 'integer'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    dir_in_anndata = "intermediates/2502/250409_anndata",
    dir_out_scvelo = "intermediates/2502/250501_scvelo",
    i = 1
  )
}





cell_types <- list.files(params$dir_in_anndata, pattern = "\\.h5ad$") |>
  str_subset("scvel", negate = TRUE) |>
  str_remove("\\.h5ad$")

length(cell_types)




# cell type ----

ct <- cell_types[[params$i]]

message("----  ",params$i, ": ", ct," ----")

base_in <- file.path( params$dir_in_anndata, ct )
base_out <- file.path( params$dir_out_scvelo, ct )

adata <- sc$read(paste0(base_in, ".h5ad"))

#~ number of cells ----
message("Number of cells: ", nrow(adata$X))

if(nrow(adata$X) < 50){
  message(" --> fewer than 50 cells, skip")
  quit()
}


#~ filter ----
scv$pp$filter_genes(adata)
scv$pp$normalize_per_cell(adata)

sc$pp$log1p(adata)

#~ preprocess ----
sc$tl$pca(adata)
sc$pp$neighbors(adata, n_neighbors=30L, n_pcs=30L)

scv$pp$moments(adata)

#~ velocity ----
scv$tl$recover_dynamics(adata, show_progress_bar = FALSE)

scv$tl$velocity(adata, mode = "stochastic")
scv$tl$velocity_graph(adata)

#~ gene likelihood ----
scv$tl$latent_time(adata)


scv$pl$velocity_embedding_stream(adata, basis='pca', color='latent_time', show = FALSE)
plt$savefig(paste0(base_out, "_velocity.png"), dpi = 300)
plt$close()

py_save_object(adata, filename = paste0(base_out, "_adata.pkl"))



scvelo_fit <- adata$var['fit_likelihood'] |>
  rownames_to_column("gene_id") |>
  mutate(cell_type = ct,
         gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1) |>
  as_tibble() |>
  arrange(desc(fit_likelihood))

qs::qsave(scvelo_fit, paste0(base_out, "_scvelo_fit.qs"))


message("======  Package versions  ======")

sessioninfo::session_info()

py_list_packages()[,1:2]



