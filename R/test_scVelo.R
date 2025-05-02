

# Inits ----

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






# load scvelo ----

dir_out_scvelo <- "intermediates/2502/250501_scvelo"


res_scvelo <- list.files(dir_out_scvelo, pattern = "_scvelo_fit\\.qs$") |>
  str_remove("_scvelo_fit\\.qs$") |>
  map_dfr(\(ct){
    qs::qread(file.path(dir_out_scvelo, paste0(ct, "_scvelo_fit.qs")))
  })


table(is.nan( res_scvelo$fit_likelihood ))


#~ Compare clustering ----

dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"
cluster_results <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv"))




both <- cluster_results |>
  full_join(res_scvelo,
            by = join_by(gene_name, cell_type))

table(clust = is.na(both$shape), scvelo = is.na(both$gene_id))
table(clust = is.na(both$shape), scvelo = is.nan(both$fit_likelihood))


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




manual_amphso <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

manual_several <- readxl::read_excel("intermediates/2502/250330_step2/manual_annotation.xlsx") |>
  dplyr::select(-amplitude, -dev_expl)

inner_join(manual_several, manual_amphso, by = join_by(gene_name, cell_type))

manual <- bind_rows(
  manual_several |> filter(manual %in% c("yes","no")),
  manual_amphso |> filter(manual %in% c("yes","no")) |> select(-descr)
) |> distinct()

stopifnot(anyDuplicated(paste0(manual$cell_type, manual$gene_name)) == 0L)



both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  left_join(manual,
            by = join_by(gene_name, cell_type)) |>
  mutate(rank = row_number(),
         manual = if_else(is.na(manual),
                          "undef",
                          manual)) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(yes = 1, no = 1, undef = 0.3)) +
  scale_size_manual(values = c(yes = 2, no = 2, undef = 0.8)) +
  scale_color_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  facet_wrap(~cluster) +
  geom_point(aes(x = rank, y = fit_likelihood,
                 alpha = manual,
                 color = manual,
                 size = manual)) +
  geom_hline(aes(yintercept = .25),
             linetype = 'dotted')



# all cells
both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  left_join(manual,
            by = join_by(gene_name, cell_type)) |>
  mutate(rank = row_number(),
         manual = if_else(is.na(manual),
                          "undef",
                          manual)) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(yes = 1, no = 1, undef = 0.1)) +
  scale_size_manual(values = c(yes = 2, no = 2, undef = 0.8)) +
  scale_color_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  geom_point(aes(x = rank, y = fit_likelihood,
                 alpha = manual,
                 color = manual,
                 size = manual)) +
  geom_hline(aes(yintercept = .25),
             linetype = 'dotted')


both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  mutate(rank = row_number(),
         cluster = case_match(shape,
                              "pulsatile" ~ "yes",
                              "nonpulsatile"~"no",
                              .default = "undef")) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(yes = .2, no = .2, undef = 0.1)) +
  scale_size_manual(values = c(yes = 1.5, no = 1.5, undef = 0.8)) +
  scale_color_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  scale_y_sqrt() + scale_x_sqrt() +
  geom_point(aes(x = rank, y = fit_likelihood,
                 alpha = cluster,
                 color = cluster,
                 size = cluster)) +
  geom_hline(aes(yintercept = .25),
             linetype = 'dotted')



both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  mutate(rank = row_number(),
         cluster = case_match(shape,
                              "pulsatile" ~ "yes",
                              "nonpulsatile"~"no",
                              .default = "undef")) |>
  ggplot() +
  theme_classic() +
  scale_fill_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  scale_x_sqrt() +
  geom_density(aes(x = fit_likelihood,
                 fill = cluster),
               alpha = .5) +
  geom_vline(aes(xintercept = .25),
             linetype = 'dotted')




both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  left_join(manual,
            by = join_by(gene_name, cell_type)) |>
  mutate(rank = row_number(),
         manual = if_else(is.na(manual),
                          "undef",
                          manual)) |>
  filter(cell_type == "AM_PHso") |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(yes = 1, no = 1, undef = 0.01)) +
  scale_size_manual(values = c(yes = 2, no = 2, undef = 0.8)) +
  scale_color_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  geom_point(aes(x = rank, y = fit_likelihood,
                 alpha = manual,
                 color = manual,
                 size = manual)) +
  geom_hline(aes(yintercept = .25),
             linetype = 'dotted')



both |>
  arrange(desc(fit_likelihood)) |>
  left_join(manual,
            by = join_by(gene_name, cell_type)) |>
  mutate(rank = row_number(),
         manual = if_else(is.na(manual),
                          "undef",
                          manual)) |>
  filter(cell_type == "AM_PHso",
         manual != "undef")




both |>
  filter(!is.na(fit_likelihood)) |>
  arrange(desc(fit_likelihood)) |>
  mutate(rank = row_number()) |>
  filter(cell_type == "ADE_PDEso") |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "yes",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "yes",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "yes",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "no",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "no",
    .default = "undef"
  )) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(yes = 1, no = 1, undef = 0.3)) +
  scale_size_manual(values = c(yes = 2, no = 2, undef = 0.8)) +
  scale_color_manual(values = c(yes = "green4", no = "red4", undef = "grey")) +
  geom_point(aes(x = rank, y = fit_likelihood,
                 alpha = guess,
                 color = guess,
                 size = guess)) +
  geom_hline(aes(yintercept = .25),
             linetype = 'dotted')


cell_types_with_clust <- both |>
  filter(shape == "pulsatile") |>
  pull(cell_type) |>
  unique()

both |>
  filter(cell_type %in% cell_types_with_clust) |>
  mutate(pulsatile_cluster = shape == "pulsatile",
         pulsatile_scvelo = fit_likelihood > .25) |>
  select(cell_type, gene_name, starts_with("pulsatile_")) |>
  nest(.by = cell_type) |>
  deframe() |>
  imap(\(.dat, .nm){
    list(cluster = .dat[["gene_name"]][ which(.dat[["pulsatile_cluster"]]) ],
         scvelo = .dat[["gene_name"]][ which(.dat[["pulsatile_scvelo"]]) ]) |>
      eulerr::euler() |>
      plot(main = list(label = .nm, cex = .8),
           quantities = list(cex = .5),
           labels = list(cex = .5))
  }) |>
  patchwork::wrap_plots()






# _____ ----
stop()
# Manual ----

# the gist of the below was run as dSQ in `step2b_scvelo_cell_type`, but here more explorations of the adata objects

#~ inits ----
# Sys.getenv("RETICULATE_PYTHON_FALLBACK")
# Sys.getenv("RETICULATE_PYTHON")
# Sys.getenv("RETICULATE_PYTHON_ENV")
# 
# Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")

library(reticulate)
# virtualenv_create("scvelo", packages = NULL)
# use_virtualenv("scvelo", required = TRUE)
# 
# virtualenv_install("scvelo",
#                    packages = c("pandas==1.3.5", "scvelo", "igraph","louvain", "scanpy"))

use_virtualenv("scvelo", required = TRUE)



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
ct <- "ILso"



#~ load 1 ct ----
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

# adata$write_h5ad(file.path(dir_out_anndata, paste0(ct, "_scvelo.h5ad") ))
# adata <- sc$read(file.path(dir_out_anndata, paste0(ct, "_scvelo.h5ad") ))



scv$tl$velocity(adata, mode = "dynamical")
scv$tl$velocity_graph(adata, show_progress_bar = FALSE)

scv$pl$velocity_embedding_stream(adata, basis='pca')




# Check individual genes
ex_genes <- c('grl-18',  'col-53', 'lgg-1', 'cutl-16') |> s2i(gids)
scv$pl$velocity(adata, ex_genes, ncols=2)


scv$pl$scatter(adata, ex_genes[[1]], color=c('velocity'))



#~ latent time ----

scv$tl$latent_time(adata)

scv$pl$scatter(adata, color='latent_time', color_map='gnuplot', size=80)


scv$pl$velocity_embedding_stream(adata, basis='pca', color='latent_time')

plot(sort(adata$var['fit_likelihood']$fit_likelihood, decreasing = TRUE))

scvelo_fit <- adata$var['fit_likelihood'] |>
  rownames_to_column("gene_id") |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         .before = 1) |>
  as_tibble() |>
  arrange(desc(fit_likelihood))

top_genes <- adata$var['fit_likelihood'] |>
  arrange(desc(fit_likelihood)) |>
  slice_head(n = 200) |>
  rownames()

i2s(top_genes, gids)

scv$pl$heatmap(adata, var_names=top_genes, sortby='latent_time', n_convolve=100L)


scv$pl$scatter(adata, basis=top_genes[1:10], ncols = 5L, color = 'latent_time')


top_genes[1:10] |> i2s(gids)


bad_genes <- scvelo_fit |> filter(fit_likelihood < .2) |> pull(gene_name) |> sample(10)

scv$pl$scatter(adata, basis=s2i(bad_genes, gids), ncols = 5L, color = 'latent_time')


#~ rank genes ----
# cluster-specific velocity
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
table(clust = is.na(both$shape), scvelo = is.na(both$gene_id))
table(clust = is.na(both$shape), scvelo = is.nan(both$fit_likelihood))




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












