# Inits ----

library(tidyverse)

dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"
dir_clust <- "intermediates/2502/250624_cluster"
dir_out   <- "intermediates/2502/260527_decomposition"

# Load ----

cell_types_osc <- qs::qread(file.path(dir_step3, "cell_types.qs")) |>
  filter(p_coherence_adj < .05,
         perplexity > 30) |>
  pull(cell_type)

clust <- read_csv(file.path(dir_clust, "250624_cluster_results.csv"),
                  show_col_types = FALSE) |>
  mutate(cell_type = if_else(cell_type == "coelomyocyte", "coelomocyte", cell_type))

stopifnot(all(
  cell_types_osc %in% unique(clust$cell_type)
))



clust_osc <- clust |>
  filter(cell_type %in% cell_types_osc)



n_puls_per_gene <- tibble(
  level = fct_inorder(c("C7", "C7+8", "all expressed")),
  puls_clusters = list(7, c(7, 8), 1:9)
) |>
  mutate(
    per_gene = map(puls_clusters, \(.cl)
                   clust_osc |>
                     filter(cluster %in% .cl) |>
                     distinct(gene_name, cell_type) |>
                     count(gene_name, name = "n_ct_pulsatile")
    )
  ) |>
  select(level, per_gene) |>
  unnest(per_gene)



# sharing statistics across the three definitions of "pulsatile"
sweep_summary <- n_puls_per_gene |>
  group_by(level) |>
  summarize(
    total_pulsatile = n(),
    frac_exactly_1  = mean(n_ct_pulsatile == 1),
    frac_5_or_more  = mean(n_ct_pulsatile >= 5),
    n_in_all_ct     = sum(n_ct_pulsatile == length(cell_types_osc))
  )


# Fig 5E-style histogram with the three thresholds overlaid
n_puls_per_gene |>
  mutate(n_ct_pulsatile = as.factor(n_ct_pulsatile)) |>
  count(level, n_ct_pulsatile,
        name = "nb_genes",
        .drop = FALSE) |>
  ggplot() +
  theme_classic() +
  labs(x = "# cell types in which gene is pulsatile",
       y = "# genes",
       fill = "Pulsatile\nclusters") +
  scale_fill_manual(values = c(C7 = "#4ddb3a", `C7+8` = "#fbc589", `all expressed` = "black")) +
  geom_col(aes(x = n_ct_pulsatile, y = nb_genes, fill = level),
           position = "dodge")


