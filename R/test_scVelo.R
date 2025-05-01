# Inits ----
Sys.getenv("RETICULATE_PYTHON_FALLBACK")
Sys.getenv("RETICULATE_PYTHON")
Sys.getenv("RETICULATE_PYTHON_ENV")

Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")

library(reticulate)
virtualenv_create("scvelo", packages = NULL)
use_virtualenv("scvelo", required = TRUE)

virtualenv_install("scvelo",
                   packages = c("pandas==1.3.5", "scvelo", "igraph","louvain", "scanpy"))

use_virtualenv("scvelo", required = TRUE)

py_config()

# py_require("scvelo")
# # py_require("scvelo@git+https://github.com/theislab/scvelo@main")
# py_require("tqdm")
# py_require("igraph")
# py_require("louvain")
# # py_require("pandas==1.3.5")

scv <- import("scvelo")
sc <- import("scanpy")
# np <- import("numpy", convert = FALSE)


dir_out_anndata <- "intermediates/2502/250409_anndata"

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


cell_types <- list.files(dir_out_anndata, pattern = "\\.h5ad$") |>
  str_subset("scvel", negate = TRUE) |>
  str_remove("\\.h5ad$")

length(cell_types)

for(ct in cell_types){
  
  message("----  ", ct," ----")
  
  basename <- file.path( dir_out_anndata, ct )
  
  data <- sc$read(paste0(basename, ".h5ad"))
  
  #~ filter ----
  scv$pp$filter_genes(adata)
  scv$pp$normalize_per_cell(adata)
  
  sc$pp$log1p(adata)
  
  #~ preprocess ----
  sc$tl$pca(adata)
  sc$pp$neighbors(adata, n_neighbors=30L, n_pcs=30L)
  
  scv$pp$moments(adata)
  
  #~ velocity ----
  scv$tl$recover_dynamics(adata)
  
  scv$tl$velocity(adata, mode = "stochastic")
  scv$tl$velocity_graph(adata)
  
  #~ gene likelihood ----
  scv$tl$latent_time(adata)
  
  scv$pl$velocity_embedding_stream(adata, basis='pca', color='latent_time', show = FALSE)
  plt$savefig(paste0(basename, "_velocity.png"), dpi = 300)
  plt$close()
  
  py_save_object(adata, filename = paste0(basename, "_adata.pkl"))
  
  scvelo_fit <- adata$var['fit_likelihood'] |>
    rownames_to_column("gene_id") |>
    mutate(cell_type = ct,
           gene_name = i2s(gene_id, gids, warn_missing = TRUE),
           .before = 1) |>
    as_tibble() |>
    arrange(desc(fit_likelihood))
  
  qs::qsave(scvelo_fit, paste0(basename, "_scvelo_fit.qs"))
}


scv$pl$velocity_embedding_stream(adata, basis='pca', color='latent_time', show = FALSE)
plt <- import("matplotlib.pyplot")
plt$savefig("velocuty.png", dpi = 300)
plt$close()

# Manual ----
adata <- sc$read(file.path( dir_out_anndata, paste0(ct, ".h5ad") ))


scv$pl$proportions(adata)


# keep minimal filtering
scv$pp$filter_genes(adata)
scv$pp$normalize_per_cell(adata)

sc$pp$log1p(adata)


sc$tl$pca(adata)
sc$pp$neighbors(adata, n_neighbors=30L, n_pcs=30L)

scv$pp$moments(adata)

scv$tl$recover_dynamics(adata)

# adata$write_h5ad(file.path(dir_out_anndata, paste0(ct, "scvelo.h5ad") ))
# adata <- sc$read(file.path(dir_out_anndata, paste0(ct, "scvelo.h5ad") ))



scv$tl$velocity(adata, mode = "stochastic")
scv$tl$velocity_graph(adata)

scv$pl$velocity_embedding_stream(adata, basis='pca')




# Check individual genes
ex_genes <- c('grl-18',  'col-53', 'lgg-1', 'cutl-16') |> s2i(gids)
scv$pl$velocity(adata, ex_genes, ncols=2)


scv$pl$scatter(adata, ex_genes[[1]], color=c('velocity'))

#~ rank genes ----
scv$tl$rank_velocity_genes(adata, min_corr=.3)


py_run_string("import scvelo as scv")
py_run_string("import pandas as pd")

py$adata <- adata
py_run_string("ranked_genes_df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])")
df = py$ranked_genes_df
dim(df)
df |>
  mutate(rank = row_number()) |>
  pivot_longer(-rank,
               names_to = "cluster",
               values_to = "gene_id") |>
  mutate(gene_name = i2s(gene_id, gids))

i <- 0

i <- i+1
scv$pl$scatter(adata, df[1:5,i], color=c('velocity'))
i2s(df[1:5,i], gids)

df[1:4,1:5]


#~ kinetic rate parameters ----

# with the uv version
# df <- lapply(adata$var, \(.np) .np$tolist() ) |>
#   bind_cols() |>
#   add_column(gene_id = adata$var_names,
#              .before = 1) |>
#   mutate(gene_name = i2s(gene_id, gids),
#          .before = 1)

# with the venv version

df <- adata$var |>
  rownames_to_column("gene_id") |>
  mutate(gene_name = i2s(gene_id, gids),
         .before = 1) |>
  as_tibble()



df2 <- df |>
  filter(fit_likelihood > .1,
         velocity_genes)

patchwork::wrap_plots(
  df2 |>
    ggplot() +
    theme_classic() +
    scale_x_log10() +
    geom_histogram(aes(x = fit_alpha)),
  df2 |>
    ggplot() +
    theme_classic() +
    scale_x_log10() +
    geom_histogram(aes(x = fit_beta * fit_scaling)),
  df2 |>
    ggplot() +
    theme_classic() +
    scale_x_log10() +
    geom_histogram(aes(x = fit_gamma))
)


#~ latent time ----

scv$tl$latent_time(adata)

scv$pl$scatter(adata, color='latent_time', color_map='gnuplot', size=80)


scv$pl$velocity_embedding_stream(adata, basis='pca', color='latent_time')

plot(sort(adata$var['fit_likelihood']$fit_likelihood, decreasing = TRUE))

scvelo_fit <- adata$var['fit_likelihood'] |>
  rownames_to_column("gene_id") |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1) |>
  as_tibble()

top_genes <- adata$var['fit_likelihood'] |>
  arrange(desc(fit_likelihood)) |>
  slice_head(n = 200) |>
  rownames()

i2s(top_genes, gids)

scv$pl$heatmap(adata, var_names=top_genes, sortby='latent_time', n_convolve=100L)


scv$pl$scatter(adata, basis=top_genes[1:15], ncols=5L, frameon=FALSE)



# Compare clustering ----

# From clustering
dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"

cluster_results <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv"))

cluster_results |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(prop_puls)

both <- cluster_results |>
  filter(cell_type == ct) |>
  full_join(scvelo_fit,
            by = join_by(gene_name))

table(is.na(both$cell_type), is.na(both$fit_likelihood))


both |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = as.factor(cluster), y = fit_likelihood))


both |>
  filter(cluster == 5,
         fit_likelihood >= .3)


both |>
  filter(cluster == 6,
         fit_likelihood <= .1)


both |> filter(is.na(cell_type), is.na(fit_likelihood))












