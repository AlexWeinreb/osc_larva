# Inits ----


library(tidyverse)
library(Seurat)
library(wbData)
library(RAPToR)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

dir_prefilt <- "intermediates/2502/250304_filt_ds"
# dir_out <- "intermediates/2502/250305_per_condition"



### Load ----
samples_here <- tibble(
  sample_name = read_tsv("data/samples_table.tsv",
                         col_types = "cffff") |>
    pull(sample_name),
  file_path = file.path(dir_prefilt,
                        sample_name |>
                          paste0(".qs") )
)



mat_list <- samples_here$file_path |>
  map(~{
    qs::qread(.x) |>
      GetAssayData(assay = "RNA", layer = "counts") |>
      rowSums() |>
      as.matrix()
  })


X <- do.call(cbind, mat_list)
colnames(X) <- samples_here$sample_name
rownames(X) <- str_replace(rownames(X), "E-BE45912.2", "E_BE45912.2")
rownames(X) <- s2i(rownames(X), gids, warn_missing = TRUE)


ref <- prepare_refdata("Cel_larv_YA", "wormRef", 600)

ae_X <- ae(X, ref)

# check output
summary(ae_X)
plot(ae_X) # plot all sample estimates
plot_cor(ae_X) # plot individual correlation profiles of samples

# get results
ae_X$age.estimates


comp <- ae_X$age.estimates |>
  as.data.frame() |>
  add_column(time_on_food = c(
    26,26,26,51,51,NA,NA,44,44,48,
    49.5,46.2,42.5,48.5,46.5,61,63
  )) |>
  add_column(promoter = rep(c("mir-228","mam-5", "grl-18"), times = c(9,1,7)))

comp |>
  ggplot() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  aes(x = time_on_food, y = age.estimate, color = promoter) +
  geom_smooth(method = "lm") +
  geom_point() +
  geom_errorbar(aes(ymin = lb, ymax = ub))















