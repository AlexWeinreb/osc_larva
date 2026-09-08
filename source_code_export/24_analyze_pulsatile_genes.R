

# Inits ----
library(tidyverse)

library(wbData)


gids <- wb_load_gene_ids(295) |>
  add_row(X = "a",
          gene_id = "nsIs198",
          symbol = "GFP",
          sequence = "GFP",
          status = "Live",
          biotype = "protein_coding_gene",
          name = "GFP"
  )


opar <- par(no.readonly = TRUE)



out_dir <- "presentations/figures/250624_puls_genes"


dir_assembled <- "intermediates/2502/250605_assembled/"

dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"
dir_clust <- "intermediates/2502/250624_cluster"


# Load ----


# Clusters and whether they are oscillatory according to bulk
# bulk_coherence <- qs::qread(file.path(dir_assembled, "250610_coherence_by_ct.qs"))


osc_raw <- readxl::read_excel("data/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID, refresh = Inf),
         gene_name = i2s(gene_id, gids) )


#> There are 21 genes that are "suppressed" or "dead" (e.g. F41B5.6, F44E2.2, ...)
#>  and end up with gene_name == NA. These can be removed.
#>  
#>  
#> In addition, there are duplicated gene names, removed with slice_sample();
#> this seems mostly due to "dead" genes in the annotation used by Meeuse 2020 compared to now
#> for ex, ZK1151.1 and T02G6.1 are "Osc" and T02G6.2 nonOsc, but since 2019 all 3 are vab-10.
#> to look at the genes:
#> 
# duplicated_genes <- osc_raw$gene_name[duplicated(osc_raw$gene_name)] |> setdiff(NA_character_)
# map(duplicated_genes,
#     ~ filter(osc_raw, gene_name == .x))
#> 
#> How to handle them? One approach is to take one randomly, the other to take the highest
#> There are only a handful of genes (6) where that might matter:
# duplicated_genes_amp_diff <- map(duplicated_genes,
#     ~ filter(osc_raw, gene_name == .x)) |>
#   map_dfr(\(.df){
#     amp <- .df$OscAmplitude
#     amp[is.na(amp)] <- 0
#     
#     data.frame(gene_name = .df$gene_name[[1]],
#                min_amp = min(amp),
#                max_amp = max(amp))
#   })
# 
# duplicated_genes_amp_diff |>
#   ggplot() + theme_classic() +
#   geom_point(aes(x = min_amp, y = max_amp),
#              alpha = .3, size = 2) +
#   geom_hline(aes(yintercept = 1.5)) + geom_vline(aes(xintercept = 1.5))
# 
# duplicated_genes_amp_diff |>
#   filter(min_amp < 1.5, max_amp > 1.5)




osc_table <- osc_raw |>
  filter(! is.na(gene_name)) |>
  mutate(osc_amplitude = if_else(is.na(OscAmplitude), 0, OscAmplitude)) |>
  select(gene_name, bulk_class = Class, osc_amplitude, bulk_peak = PeakPhase) |>
  group_by(gene_name) |>
  slice_max(osc_amplitude,
            with_ties = FALSE) |>
  ungroup()


stopifnot(!any(is.na(osc_table$gene_name)))
stopifnot(anyDuplicated(osc_table$gene_name) == 0L)




#~ Load results ----

cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- read_csv(file.path(dir_clust, "250624_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(cell_type %in% cell_types_info$cell_type)



stopifnot(all.equal(
  cell_types_info$cell_type |> unique() |> sort(),
  all_genes$cell_type |> unique() |> sort()
))

# check which cell types we're filtering out
# setdiff(
#   all_genes$cell_type |> unique(),
#   cell_types_info$cell_type |> unique()
# )



#~ PANTHER ----
dict_panther <- qs::qread("data/gene_families/panther_dict.qs")
nb_genes_in_fam <- dict_panther[,-1] |> colSums()
dict_panther <- dict_panther[ , c(TRUE, nb_genes_in_fam > 5) ]






# General plots ----




#~ Plot nb of puls ----

cell_types_info |>
  mutate(cell_type_noneur = if_else(tissue == "neuron", "", cell_type),
         bulk = if_else(p_coherence_adj < .05, "oscillatory", "nonoscillatory"),
         prop_puls = prop_puls/100) |>
  ggplot() +
  theme_classic() +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_continuous(labels = scales::label_comma(),
                     limits = c(0,NA)) +
  scale_shape_manual(values = c(19,8)) +
  scale_size_manual(values = c(1, 2.5)) +
  aes(x = n_tot,
      y = prop_puls,
      color = tissue) +
  xlab("Number of genes tested") +
  ylab("Proportion pulsatile") +
  geom_point(aes(shape = bulk, size = bulk)) +
  ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)


# ggsave("nb_puls_genes_by_celltype.png",
#        path = out_dir,
#        width = 20, height = 15, units = "cm")






ordered_cells <- cell_types_info |>
  arrange(tissue, desc(prop_puls)) |>
  pull(cell_type) |> fct_inorder() |> levels()





cell_types_info |>
  mutate(n_nonpuls = n_tot - n_puls) |>
  select(- n_tot) |>
  pivot_longer(cols = starts_with("n_"),
               names_to = "type",
               values_to = "count",
               names_prefix = "n_") |>
  mutate(cell_type = factor(cell_type, levels = ordered_cells)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c(puls = "red4", nonpuls = "grey70")) +
  xlab(NULL) + ylab("Number of genes") +
  geom_col(aes(x = cell_type, y = count, fill = type),
           show.legend = FALSE) +
  geom_text(aes(x = cell_type,
                y = nb_puls,
                label = prop_puls),
            data = all_genes |>
              summarize(nb_puls = sum(shape == "pulsatile"),
                        nb_tested = n(),
                        .by = cell_type) |>
              mutate(prop_puls = round( 100 * nb_puls / nb_tested) |> paste0("%") ),
            nudge_y = 100
  )









# Compare clusters to bulk ----


# Note: this reproduces the plot generated in `22_analyze_plot_hclust_results.R`, 
# that one was used as figure, this is just a check.

cell_types_info |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence") +
  ylab("Perplexity") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = mean_coherence, y = perplexity, color = tissue,
                 shape = p_coherence_adj < .05),
             size = 3)



cell_types_osc <- cell_types_info |>
  filter(p_coherence_adj < .05,
         perplexity > 30) |>
  pull(cell_type)






# Compare genes to known bulk ----

#~ overlap peaky/Osc ----

osc_genes_compare <- all_genes |>
  filter(cell_type %in% cell_types_osc) |>
  select(gene_name, shape) |>
  summarize(is_puls = any(shape == "pulsatile"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(`single-cell` = if_else(is_puls, "pulsatile", "nonpulsatile"))

osc_genes_compare |>
  (\(df) table(bulk = df$bulk_class, `single-cell` = df$`single-cell`,
               useNA = "ifany"))()

osc_genes_compare |>
  (\(df) table(bulk = if_else(!is.na(df$bulk_class) & df$bulk_class == "Osc",
                              "Osc",
                              "nonOsc"),
               `single-cell` = df$`single-cell`))()


# pdf(file.path(out_dir, "osc_vs_pulsatile_euler.pdf"),
#     width = 2, height = 2)
list(
  bulk = osc_genes_compare$gene_name[which(osc_genes_compare$bulk_class == "Osc")],
  sc = osc_genes_compare$gene_name[which(osc_genes_compare$`single-cell` == "pulsatile")]
) |>
  eulerr::euler() |>
  plot(quantities = TRUE,
       edges = FALSE)

# dev.off()

# > not used in figures



# Check genes Osc but non-pulsatile
genelist_bulk_nonpuls <- osc_genes_compare |>
  filter(bulk_class == "Osc",
         `single-cell` == "nonpulsatile") |>
  pull(gene_name)

length(genelist_bulk_nonpuls)

dict <- wormbaseEnrich::fetch_dictionary("tissue", cache = Inf)
dict <- wormbaseEnrich::fetch_dictionary("go", cache = Inf)

enrres <- wormbaseEnrich::enrichment_analysis(
  genelist_bulk_nonpuls |> s2i(gids),
  dict
)

wormbaseEnrich::plot_enrichment_results(enrres)

enrres |>
  mutate(is_devt = str_detect(term_name,
                              "^(C|AB)[aplrvd]+$")) |>
  ggplot() +
  theme_classic() +
  scale_alpha_manual(values = c(`TRUE` = .03, `FALSE` = 1)) +
  geom_point(aes(x = observed, y = -log10(FDR),
                 size = enrichment_fc,
                 alpha = is_devt,
                 color = is_devt))





## Compare OscAmplitude (cf also hclust_results)

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  ggplot() +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.8,.7),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  xlab("Bulk-annotated amplitude") +
  geom_density(aes(x = osc_amplitude, fill = `single-cell`),
               alpha = .5)

# ggsave("osc_vs_pulsatile_density.pdf",
#        path = out_dir,
#        width = 75, height = 45, units = "mm")



# > +++ Fig. 4F +++ ----

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  ggplot() +
  theme_classic() +
  theme(
    legend.position = "none",
    # legend.position.inside = c(.8,.7),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  ylab("Bulk-annotated amplitude") +# xlab(NULL) +
  scale_fill_manual(values = c(nonpulsatile = "#b4a2ce", pulsatile = "#bb7e76")) +
  geom_boxplot(aes(x = `single-cell`, y = osc_amplitude,
                   fill = `single-cell`))
  # ggbeeswarm::geom_quasirandom(aes(x = `single-cell`, y = osc_amplitude,color = `single-cell`),
  #              alpha = .5)

# ggsave("osc_vs_pulsatile_boxplot.pdf",
#        path = out_dir,
#        width = 50, height = 45, units = "mm")

## stat test

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  wilcox.test(osc_amplitude ~ `single-cell`, data = _)

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  summarize(median_amplitude = median(osc_amplitude),
            mad = mad(osc_amplitude),
            .by = `single-cell`)














# intersections of genes ----



#~ proportion of expressed genes ----

in_osc_ct <- all_genes |>
  filter(cell_type %in% cell_types_osc) |>
  mutate(is_puls = (shape == "pulsatile")) |>
  select(cell_type, gene_name, is_puls)



# proportion of testable genes (high genes) that are peaky in both
intersections <- expand_grid(
  cell_1 = unique(in_osc_ct$cell_type),
  cell_2 = unique(in_osc_ct$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(
        in_osc_ct$gene_name[in_osc_ct$cell_type == .x &
                              in_osc_ct$is_puls],
        in_osc_ct$gene_name[in_osc_ct$cell_type == .y &
                              in_osc_ct$is_puls]
      ))
    ),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(
        in_osc_ct$gene_name[in_osc_ct$cell_type == .x],
        in_osc_ct$gene_name[in_osc_ct$cell_type == .y]
      ))
    )
  )

mat_intersections <- intersections |>
  mutate(prop_inter = 100 * intersection / union) |>
  select(cell_1, cell_2, prop_inter) |>
  pivot_wider(id_cols = cell_1,
              names_from = "cell_2",
              values_from = "prop_inter") |>
  column_to_rownames("cell_1") |>
  as.matrix()

dist_mat_prop <- mat_intersections/max(mat_intersections, na.rm = TRUE)



hc <- as.dist(1-dist_mat_prop) |>
  hclust(method = "ward.D2") |>
  as.dendrogram() |>
  reorder(wts = diag(mat_intersections)) |>
  rev() |>
  as.hclust()


pheatmap::pheatmap(mat_intersections,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # width = 3.15, height = 2.5,
                   # filename = file.path(out_dir, "intersections_expr.png"),
                   annotation_row = cell_types_info |> filter(cell_type %in% rownames(mat_intersections)) |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_col = cell_types_info |> filter(cell_type %in% rownames(mat_intersections)) |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_colors = list(tissue = c(
                     glia = scales::hue_pal()(8)[[1]],
                     # muscle = scales::hue_pal()(8)[[2]],
                     # neuron = scales::hue_pal()(8)[[3]],
                     other = scales::hue_pal()(8)[[4]],
                     pharynx = scales::hue_pal()(8)[[5]],
                     # reproductive = scales::hue_pal()(8)[[6]],
                     skin = scales::hue_pal()(8)[[7]]
                   )),
                   main = "Proportion of pulsatile among expressed genes")


# dev.off()




#~ Prop puls in both among puls in one ----
# non-symmetric


intersections <- expand_grid(
  cell_1 = unique(in_osc_ct$cell_type),
  cell_2 = unique(in_osc_ct$cell_type)
) |>
  mutate(
    in_1 = map_int(cell_1,
                   \(.ct1) sum(in_osc_ct$is_puls[in_osc_ct$cell_type == .ct1]) ),
    intersection = map2_int(
      cell_1, cell_2,
      \(.ct1, .ct2){
        intersect(
          in_osc_ct$gene_name[in_osc_ct$cell_type == .ct1 &
                                in_osc_ct$is_puls],
          in_osc_ct$gene_name[in_osc_ct$cell_type == .ct2 &
                                in_osc_ct$is_puls]
        ) |> length()
      })
  )

mat_intersections <- intersections |>
  mutate(prop_inter = 100 * intersection / in_1) |>
  select(cell_1, cell_2, prop_inter) |>
  pivot_wider(id_cols = cell_1,
              names_from = "cell_2",
              values_from = "prop_inter") |>
  column_to_rownames("cell_1") |>
  as.matrix()

stopifnot(!any( is.na(mat_intersections) ))

diag(mat_intersections) <- NA_real_



dist_mat_prop <- mat_intersections/max(mat_intersections, na.rm = TRUE)

hc <- as.dist(1-dist_mat_prop) |>
  hclust(method = "ward.D2")


pheatmap::pheatmap(mat_intersections,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # width = 3.15, height = 2.5,
                   # filename = file.path(out_dir, "intersections_pulsFirst.png"),
                   annotation_row = cell_types_info |> filter(cell_type %in% rownames(mat_intersections)) |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_col = cell_types_info |> filter(cell_type %in% rownames(mat_intersections)) |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_colors = list(tissue = c(
                     glia = scales::hue_pal()(8)[[1]],
                     # muscle = scales::hue_pal()(8)[[2]],
                     # neuron = scales::hue_pal()(8)[[3]],
                     other = scales::hue_pal()(8)[[4]],
                     pharynx = scales::hue_pal()(8)[[5]],
                     # reproductive = scales::hue_pal()(8)[[6]],
                     skin = scales::hue_pal()(8)[[7]]
                   )),
                   main = "Prop pulsatile in both (among puls in first)")


# dev.off()





#~ proportion of puls genes in both (out of puls genes in each) ----

intersections <- expand_grid(
  cell_1 = unique(in_osc_ct$cell_type),
  cell_2 = unique(in_osc_ct$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(
        in_osc_ct$gene_name[in_osc_ct$cell_type == .x &
                              in_osc_ct$is_puls],
        in_osc_ct$gene_name[in_osc_ct$cell_type == .y &
                              in_osc_ct$is_puls]
      ))),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(
        in_osc_ct$gene_name[in_osc_ct$cell_type == .x &
                              in_osc_ct$is_puls],
        in_osc_ct$gene_name[in_osc_ct$cell_type == .y &
                              in_osc_ct$is_puls]
      )))
  )

mat_intersections <- intersections |>
  mutate(prop_inter = 100 * intersection / union) |>
  select(cell_1, cell_2, prop_inter) |>
  pivot_wider(id_cols = cell_1,
              names_from = "cell_2",
              values_from = "prop_inter") |>
  column_to_rownames("cell_1") |>
  as.matrix()

# remove empty cell types, remove diagonal
# mat_intersections <- mat_intersections[!is.na(rowSums(mat_intersections)),
#                                        !is.na(colSums(mat_intersections))]
diag(mat_intersections) <- NA_real_



dist_mat_prop <- mat_intersections/max(mat_intersections, na.rm = TRUE)

hc <- as.dist(1-dist_mat_prop) |>
  hclust(method = "ward.D2")


# > +++ Fig. 5D +++ ----

pheatmap::pheatmap(mat_intersections,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # width = 6, height = 5,
                   # fontsize = 7,
                   # filename = file.path(out_dir, "intersections_pulsOne.pdf"),
                   annotation_row = cell_types_info |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_col = cell_types_info |> column_to_rownames("cell_type") |> select(tissue),
                   annotation_colors = list(tissue = c(
                     glia = scales::hue_pal()(8)[[1]],
                     muscle = scales::hue_pal()(8)[[2]],
                     neuron = scales::hue_pal()(8)[[3]],
                     other = scales::hue_pal()(8)[[4]],
                     pharynx = scales::hue_pal()(8)[[5]],
                     reproductive = scales::hue_pal()(8)[[6]],
                     skin = scales::hue_pal()(8)[[7]]
                   )),
                   main = "Prop pulsatile in both (among puls in one)")


# dev.off()




#~ Intersection sizes ----
# > +++ Fig. 5E +++ ----

in_osc_ct |>
  filter(is_puls) |>
  summarize(nb_cell_types = n(),
            .by = gene_name) |>
  summarize(nb_genes = n(),
            .by = nb_cell_types) |>
  mutate(nb_cell_types = factor(as.character(nb_cell_types),
                                levels = rev(1:13))) |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 90,
                               vjust = .5,
                               hjust = 1),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  xlab("Number of genes") +
  ylab("Number of cell types in which pulsatile") +
  scale_x_continuous(n.breaks = 7,
                     labels = scales::label_comma()) +
  scale_y_discrete(breaks = seq(1,14,by = 2)) +
  geom_col(aes(y = nb_cell_types, x = nb_genes))


# ggsave("puls_per_ct_intersections.pdf",
#        path = out_dir,
#        width = 30, height = 70, units = "mm")



in_osc_ct |>
  filter(is_puls) |>
  summarize(nb_cell_types = n(),
            .by = gene_name) |>
  summarize(nb_genes = n(),
            .by = nb_cell_types) |>
  # filter(nb_cell_types <= 1) |>
  filter(nb_cell_types >= 5) |>
  pull(nb_genes) |>
  sum()

in_osc_ct |>
  summarize(nb_cell_types = sum(is_puls),
            .by = gene_name) |>
  arrange(desc(nb_cell_types)) |>
  View()




# Gene families ----


in_osc_ct <- all_genes |>
  filter(cell_type %in% cell_types_osc)


#~ GO ----

dict <- wormbaseEnrich::fetch_dictionary("go", cache = Inf)

all_go <- cell_types_osc |>
  set_names() |>
  map(\(ct){
    ct_genes <- all_genes |>
      filter(cell_type == ct) |>
      mutate(gene_id = s2i(gene_name, gids))
    
    gene_list <- ct_genes$gene_id[ct_genes$shape == "pulsatile"] |> intersect(dict$wbid)
    background_genes <- ct_genes$gene_id |> intersect(dict$wbid)
    
    
    enr_go <- wormbaseEnrich::enrichment_analysis(
      gene_list = gene_list,
      dictionary = dict,
      background_genes = background_genes
    )
    enr_go |>
      filter(observed > 5, FDR < .1) |>
      pull(term_name)
  })

imap_dfr(all_go,
         ~enframe(.x, name = NULL) |> add_column(ct = .y)) |> as.data.frame()








#~ PANTHER families ----


all_panther <- map_dfr(
  cell_types_osc,
  \(ct){
    
    ct_genes <- all_genes |>
      filter(cell_type == ct) |>
      mutate(gene_id = s2i(gene_name, gids))
    
    gene_list <- ct_genes$gene_id[ct_genes$shape == "pulsatile"] |> intersect(dict_panther$wbid)
    background_genes <- ct_genes$gene_id |> intersect(dict_panther$wbid)
    
    
    enr_res <- wormbaseEnrich::enrichment_analysis(
      gene_list = gene_list,
      dictionary = dict_panther,
      background_genes = background_genes,
      filter = FALSE
    )
    
    enr_res |>
      select(-term_id) |>
      add_column(cell_type = ct, .before = 1) |>
      mutate(term_name = if_else(startsWith(term_name, " PTHR"),
                                 paste0("Unnamed", term_name),
                                 term_name)) |>
      separate_wider_regex(term_name,
                           patterns = c(family_terms = "^[[:print:] ]+", " ",
                                        family_id = "PTHR[[:digit:]]+(?:\\:SF[[:digit:]]+)?$"))
  }) |>
  filter(expected > 0) |>
  mutate(FDR = p.adjust(p_value, method = "BH"))




all_panther |>
  filter(observed >= 5) |>
  mutate(family_terms = if_else(family_terms == "Unnamed", family_id, family_terms)) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~cell_type) +
  aes(x = observed, y = -log10(FDR), label = family_terms) +
  geom_point() +
  ggrepel::geom_text_repel()



panther_filt <- all_panther |>
  filter(observed >= 5,
         FDR < .1)


# # manual annotation based on the genes in each family
# i <- 0
# 
# i <- i+1
# panther_filt[i,]
# col_nb <- str_detect(colnames(dict_panther), paste0(panther_filt$family_id[i], " WBbt:0000000"))
# colnames(dict_panther)[col_nb]
# 
# genelist <- dict_panther$wbid[dict_panther[col_nb] == 1L]
# genelist |> i2s(gids)



# > +++ Table EV6 +++ ----

manual_annot_panther_families <- readxl::read_excel("data/gene_families/manual_annot_panther_families.xlsx")


all(panther_filt$family_id %in% manual_annot_panther_families$family)
# panther_filt$family_id[! panther_filt$family_id %in% manual_annot_panther_families$family]

panther_filt <- left_join(panther_filt,
                          manual_annot_panther_families,
                          by = c(family_id = "family")) |>
  select(cell_type, description = short_description, family_id, expected, observed, enrichment_fc, FDR)


# when both family and subfamily, only keep main family
subfams <- str_match(panther_filt$family_id, "^(PTHR[0-9]+)\\:SF[0-9]+$")
stopifnot(all(
  na.omit(subfams[,2]) %in% panther_filt$family_id
))

panther_filt <- panther_filt |>
  filter(! family_id %in% subfams[,1])




#~~ panther res ----

# Remove underscores in cell types

panther_filt <- panther_filt |>
  mutate(cell_type = str_replace_all(cell_type, "_", " "))

panther_filt |>
  # filter(family_id %in% families_in_multiple) |>
  mutate(family = paste0(family_id,": ", description),
         logFDR = -log10(FDR + 1e-16)) |>
  pivot_wider(id_cols = family,
              names_from = cell_type,
              values_from = logFDR,
              values_fill = 1) |>
  column_to_rownames("family") |>
  as.matrix() |>
  pheatmap::pheatmap()

mat <- panther_filt |>
  # filter(family_id %in% families_in_multiple) |>
  mutate(logFDR = -log10(FDR + 1e-16)) |>
  pivot_wider(id_cols = family_id,
              names_from = cell_type,
              values_from = logFDR,
              values_fill = 1) |>
  column_to_rownames("family_id") |>
  as.matrix()

hc_fams <- dist(mat) |> hclust()
hc_cts <- dist(t(mat)) |> hclust()

pheatmap::pheatmap(mat)

panther_filt_plot <- panther_filt |>
  # filter(family_id %in% families_in_multiple) |>
  mutate(family_id = factor(family_id,
                         levels = rev(hc_fams$labels[hc_fams$order])),
         cell_type = factor(cell_type,
                            levels = hc_cts$labels[hc_cts$order]))


# > +++ Fig. 5G +++ ----

panther_filt_plot |>
  ggplot() +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "top",
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, "mm"),
    legend.margin = margin(b = 0),
    legend.box.margin = margin()
  ) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(
    low = "grey75",
    high = "orange2",
    limits = c(0, NA)
  ) +
  scale_alpha_continuous(transform = c("log10", "reverse"),
                         limits = c(1, 1e-8),
                         range = c(.2,1)) +
  geom_point(aes(
    y = cell_type, x = description,
    color = -log10(FDR + 1e-16),
    size = log2(enrichment_fc),
    alpha = FDR + 1e-16
  ),
  shape = 16)

# ggsave("panther_terms_enrichement_dotplot.pdf", path = out_dir,
#        width = 220, height = 95, units = "mm")


# With description on the right
# https://github.com/tidyverse/ggplot2/issues/3171#issuecomment-1340221650
# note, sec.axis should become possible with next version of ggplot2 (3.5.1 doesn't have it)







#~ manually defined families ----

genelist_hedgehog <- readxl::read_excel("data/gene_families/hedgehog_hao2006_table1.xlsx", skip = 1) |>
  pull(gene_name) |>
  wb_clean_gene_names()




genelist_zp <- readxl::read_excel("data/gene_families/zp_cohen2019.xlsx") |>
  pull(gene_id) |>
  unique()



genelist_appg <- readxl::read_excel("data/gene_families/David-Raizen_2014_bio20147500-sup-table_s5.xlsx",
                                 sheet = 1,
                                 range = readxl::cell_cols(2),
                                 col_names = "gene_name") |>
  pull(gene_name) |>
  wb_clean_gene_names()



genelist_col <- readxl::read_excel("data/gene_families/teuscher_2019_table_S1.xlsx",
                                    sheet = "Ce Matrisome",
                                    skip = 1L) |>
  filter(`Matrisome Category` == "Cuticular Collagens") |>
  pull("WormBase ID") |>
  wb_clean_gene_names()







#~| Plot proportions ----


bind_rows(
  tibble(family = "collagens",
         gene_id = genelist_col),
  tibble(family = "Hh",
         gene_id = genelist_hedgehog),
  tibble(family = "APPG",
         gene_id = genelist_appg),
  tibble(family = "ZP",
         gene_id = genelist_zp)
) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE)) |>
  inner_join(all_genes |> filter(cell_type %in% cell_types_osc),
            by = "gene_name") |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_non_puls = sum(shape != "pulsatile"),
            .by = c(family, cell_type)) |>
  arrange(family, n_puls) |>
  mutate(cell_type = fct_inorder(cell_type)) |>
  pivot_longer(-c(family, cell_type),
               names_to = "category",
               values_to = "count",
               names_prefix = "n_") |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  xlab(NULL) + ylab("Number of genes") +
  scale_fill_manual(values = c("puls" = scales::muted("red"), "non_puls" = "grey40")) +
  facet_wrap(~family) +
  geom_col(aes(x = cell_type, y = count, fill = category),
           show.legend = FALSE)

# ggsave(paste0(man_fam,"_proportions.pdf"), path = out_dir,
#        width = 40, height = 70, units = "mm")




man_fam <- "collagens"

for(man_fam in c("collagens", "hedgehog","appg","zp")){
  
  genelist <- switch (man_fam,
                      collagens = genelist_col,
                      hedgehog = genelist_hedgehog,
                      appg = genelist_appg,
                      zp = genelist_zp
  )
  
  # length(genelist)
  # table( i2s(genelist, gids) %in% all_genes$gene_name)
  # table( i2s(genelist, gids) %in% all_genes$gene_name[all_genes$shape == "pulsatile"])
  
  dat <- all_genes |>
    filter(gene_name %in% i2s(genelist, gids),
           cell_type %in% cell_types_osc) |>
    summarize(n_puls = sum(shape == "pulsatile"),
              n_non_puls = sum(shape != "pulsatile"),
              .by = cell_type) |>
    arrange(n_puls) |>
    mutate(cell_type = str_replace_all(cell_type, "_", " "),
           cell_type = fct_inorder(cell_type)) |>
    pivot_longer(-cell_type,
                 names_to = "category",
                 values_to = "count",
                 names_prefix = "n_")
  
  n_bars <- length(unique(dat$cell_type))
  
  tot_height <- 3.4*n_bars + 8.7
  
  # > +++ Fig. 5F +++ ----
  
  dat |>
    ggplot() +
    theme_classic() +
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 10),
      legend.position = "none",
      plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    coord_flip() +
    xlab(NULL) + ylab("Number of genes") +
    scale_fill_manual(values = c("puls" = scales::muted("red"), "non_puls" = "grey40")) +
    geom_col(aes(x = cell_type, y = count, fill = category))
  
  # ggsave(paste0(man_fam,"_proportions.pdf"), path = out_dir,
  #        width = 40, height = tot_height, units = "mm")
  
}











# Illustrate GAM ----

library(Seurat)
library(tidyverse)
# library(patchwork)

dir_step2 <- "intermediates/2502/250624_step2/"
dir_step1 <- "intermediates/2502/250609_step1/"

# if working from external HDD
# dir_step2 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250624_step2/"
# dir_step1 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step1/"

# dir_step2 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"
# dir_step1 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250609_step1/"

# dir_step2 <- "D:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"
# dir_step1 <- "D:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250609_step1/250609_step1/"


dir_fig_gam <- "presentations/figures/260907_gam_illustrations"
# dir.create(dir_fig_gam)




# clustering results

dir_clust <- "intermediates/2502/250624_cluster"
dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"


cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- read_csv(file.path(dir_clust, "250624_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(cell_type %in% cell_types_info$cell_type)



stopifnot(all.equal(
  cell_types_info$cell_type |> unique() |> sort(),
  all_genes$cell_type |> unique() |> sort()
))




#~ ILso ----

ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu_unsmoothed.qs")) )

mods_uncentered <- qs::qread(file.path(dir_step2, paste0("ILso", "_mods_uncentered.qs")))



smooth_centered <- list.files(dir_step2,
                              pattern = "_preds_cent_clipped\\.qs$") |>
  enframe(value = "filename",
          name = NULL) |>
  separate_wider_regex(filename,
                       patterns = c(
                         cell_type = "^.+",
                         "_preds_cent_clipped\\.qs"
                       ),
                       cols_remove = FALSE) |>
  pmap(\(cell_type, filename){
    mat_preds <- qs::qread(file.path(dir_step2,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)


# computed same for all genes
mean_sf <- lapply(mods_uncentered,
                  \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
  unlist() |>
  log() |>
  mean() |>
  exp()








#~~ genes ----

# > +++ Fig. 4A +++ ----
# > +++ Fig. EV3A +++ ----
# > +++ Fig. EV4C +++ ---- 

# 4A and EV4C use same code: 4A use selected genes and ggsave individual plots
# EV4C prepare a list of genes to plot, and assemble them before ggsave


# goi <- "grl-18"
# goi <- "nhr-23"
# goi <- "col-109"
# goi <- "pugs-11"
# goi <- "mam-5"
# 
# 
# goi <- "rps-27A"
# goi <- "dnj-1"


# goi <- rownames(ilso_subseu) |> sample(1)

tab <- inner_join(
  wormOsc::table_osc_genes |>
    select(gene_name, osc_amplitude, peak_phase_deg) |>
    mutate(oscillating = !is.na(osc_amplitude),
           high_osc_amplitude = osc_amplitude > 2),
  all_genes |>
    filter(cell_type == "ILso") |>
    select(gene_name, sc_shape = shape, sc_cluster = cluster),
  by = "gene_name"
)


# > +++ Fig. EV4B +++ ----
tab |>
  count(high_osc_amplitude, sc_shape)



set.seed(456)
genes_to_plot <- bind_rows(
  tab |>
    filter(
      high_osc_amplitude,
      sc_shape == "pulsatile"
      )
  ,
  tab |>
    filter(
      high_osc_amplitude,
      sc_shape == "nonpulsatile"
    )
  ,
  tab |>
    filter(
      !oscillating,
      sc_shape == "pulsatile"
    )
  ,
  tab |>
    filter(
      !oscillating,
      sc_shape == "nonpulsatile"
    )
) |>
  slice_sample(n = 4,
               by = c(high_osc_amplitude, oscillating, sc_shape)) |>
  pull(gene_name)


goi <- genes_to_plot[[1]]

dir_fig_gam <- "presentations/figures/260623_random_gam_plots"
# dir.create(dir_fig_gam)


# for(goi in c("grl-18", "nhr-23", "col-109", "pugs-11", "rps-27A", "dnj-1")){
for(goi in genes_to_plot){

  
  row_nb <- which(wormOsc::table_osc_genes$gene_name == goi)
  stopifnot( length(row_nb) == 1L )
  goi_ampl <- wormOsc::table_osc_genes[["osc_amplitude"]][[row_nb]]
  goi_phase <- wormOsc::table_osc_genes[["peak_phase_deg"]][[row_nb]]
  
  
  row_nb_clust <- which(all_genes$cell_type == "ILso" & all_genes$gene_name == goi)
  stopifnot( length(row_nb_clust) == 1L )
  goi_shape <- all_genes[["shape"]][[row_nb_clust]]
  goi_cluster <- all_genes[["cluster"]][[row_nb_clust]]
  
  if(is.na(goi_ampl)){
    title_bulk <- ("not oscillating")
  } else{
    title_bulk <- bquote(
        "amplitude " * .(round(goi_ampl, 1)) * ", " *
        "phase " * .(round(goi_phase)) * "°; "
    )
  }
  
  title_clust <- bquote( .(goi_shape) * " (cluster " * .(goi_cluster) * ")" )
  
  
#~| cells ----

dat <- FetchData(ilso_subseu, vars = c("PC_1","PC_2",goi))

# for ILso, invert axes for easier interpretation
dat$PC_1 <- -dat$PC_1



gg_pca <- dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 10),
    legend.position = "top",
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  labs(x = "PC 1", y = "PC 2") +
  scale_color_gradient(low = alpha("grey", .3), high = alpha("blue2", .8)) +
  ggtitle(bquote(
    italic(.(goi)) * ": " * .(title_bulk) * "; " * .(title_clust)
    )) +
  ggrastr::geom_point_rast(aes(x = PC_1, y = PC_2,
                               color = .data[[goi]]),
                           shape = 16,
                           size = 1,
                           raster.dpi = 500)



# ggsave(paste0(goi, "_expr.pdf"),
#        path = dir_fig_gam,
#        width = 55, height = 55, units = "mm",
#        scale = 1)
# 
# ggsave(paste0(goi, "_expr.png"),
#        path = dir_fig_gam,
#        width = 55, height = 55, units = "mm",
#        scale = 2)




#~| GAM ----

# clipped
mod <- mods_uncentered[[goi]]

dat <- data.frame(
  pseudotime = mod$model$pseudotime,
  count = log10( 1 + mean_sf * mod$model$expr / exp(mod$model$`offset(log(size_factors))`) ),
  fit = log10( 1 + mean_sf * mod$fitted.values / exp(mod$model$`offset(log(size_factors))`) )
)

clip <- max(
  dat$count |> quantile(probs = .95),
  1.1 * max(dat$fit)
)

gg_gam <- dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  ylab("Expression") +
  # ggtitle(goi) +
  scale_x_continuous(breaks = 0:1) +
  scale_y_continuous(limits = c(0, clip),
                     oob = scales::squish) +
  ggrastr::geom_point_rast(aes(x = pseudotime,
                               y = count),
                           alpha = .3,
                           size = 1,
                           shape = 16,
                           raster.dpi = 500) +
  geom_line(aes(x = pseudotime,
                y = fit),
            color = 'orange2',
            linewidth = 1)

# ggsave(paste0(goi, "_devexpl.pdf"),
#        path = dir_fig_gam,
#        width = 50, height = 45, units = "mm",
#        scale = 1)
# 
# ggsave(paste0(goi, "_devexpl.png"),
#        path = dir_fig_gam,
#        width = 50, height = 45, units = "mm",
#        scale = 2)


title_grob <- grid::textGrob(
  bquote(italic(.(goi)) * ": " * .(title_bulk) * "; " * .(title_clust)),
  gp = grid::gpar(fontsize = 10),
  x = 0, hjust = 0  # left-align; drop these for centered
)

plots_aligned <- cowplot::plot_grid(
  gg_pca + ggtitle(NULL), gg_gam,
  nrow = 1,
  align = "v",   # aligns vertical extents
  axis = "tb"    # matches top and bottom axes
)

gg_assembled <- patchwork::wrap_elements(plots_aligned) +
  patchwork::plot_annotation(
    title = bquote(italic(.(goi)) * ": " * .(title_bulk) * "; " * .(title_clust)),
    theme = theme(plot.title = element_text(size = 10))
  )


ggsave(paste0(goi, "_assembled.pdf"),
       plot = gg_assembled,
       path = dir_fig_gam,
       width = 70, height = 40, units = "mm",
       scale = 1.5)

}






#~ BWM ----


subseu <- qs::qread( file.path(dir_step1,
                                        paste0("BWM", "_seu_unsmoothed.qs")) )

mods_uncentered <- qs::qread(file.path(dir_step2, paste0("BWM", "_mods_uncentered.qs")))



# computed same for all genes
mean_sf <- lapply(mods_uncentered,
                  \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
  unlist() |>
  log() |>
  mean() |>
  exp()


# > +++ Fig. EV5A +++ ----
# > +++ Fig. EV5B +++ ----


#~ gene ----

goi <- "lgc-34"
goi <- "vha-11"



#~| cells ----

dat <- FetchData(subseu, vars = c("PC_1","PC_2",goi))



dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.title = element_text(face = "italic",
                              size = 10),
    legend.position = "top",
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  labs(x = "PC 1", y = "PC 2") +
  scale_color_gradient(low = alpha("grey", .3), high = alpha("blue2", .8)) +
  # ggtitle(goi) +
  ggrastr::geom_point_rast(aes(x = PC_1, y = PC_2,
                               color = .data[[goi]]),
                           shape = 16,
                           size = 1,
                           raster.dpi = 500)



# ggsave(paste0(goi, "_bwm_expr.pdf"),
#        path = dir_fig_gam,
#        width = 55, height = 55, units = "mm")




#~| GAM ----


# clipped
mod <- mods_uncentered[[goi]]

dat <- data.frame(
  pseudotime = mod$model$pseudotime,
  count = log10( 1 + mean_sf * mod$model$expr / exp(mod$model$`offset(log(size_factors))`) ),
  fit = log10( 1 + mean_sf * mod$fitted.values / exp(mod$model$`offset(log(size_factors))`) )
)

clip <- max(
  dat$count |> quantile(probs = .95),
  1.1 * max(dat$fit)
)

dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  ylab("Expression") +
  # ggtitle(goi) +
  scale_x_continuous(breaks = 0:1) +
  scale_y_continuous(limits = c(0, clip),
                     oob = scales::squish) +
  ggrastr::geom_point_rast(aes(x = pseudotime,
                               y = count),
                           alpha = .3,
                           size = 1,
                           shape = 16,
                           raster.dpi = 500) +
  geom_line(aes(x = pseudotime,
                y = fit),
            color = 'orange2',
            linewidth = 1)


# ggsave(paste0(goi, "_bwm_devexpl.pdf"),
#        path = dir_fig_gam,
#        width = 50, height = 45, units = "mm")








