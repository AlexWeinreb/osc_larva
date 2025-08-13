

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
  mutate(gene_id = wb_clean_gene_names(WB_ID),
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


# Note: the version in paper is done in `hclust_results.R`, before running step3.

cell_types_info |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence") +
  ylab("Perplexity") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 19)) +
  geom_point(aes(x = mean_coherence, y = perplexity, color = tissue,
                 shape = p_coherence_adj < .05),
             size = 3)
# note same plot was saved in step 3 (use that one for figures)


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

dict <- wormbaseEnrich::fetch_dictionary("tissue")
dict <- wormbaseEnrich::fetch_dictionary("go")

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

dict <- wormbaseEnrich::fetch_dictionary("go")

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



#~~ res ----

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


genelist_zp <- read_tsv("data/gene_families/IPR001507_ZP_proteins.tsv",
                        skip = 1L) |>
  pull(`WormBase Gene ID`) |>
  unique()

genelist_col <- read_tsv("data/gene_families/IPR008160_col.tsv",
                         skip = 1L) |>
  pull(`WormBase Gene ID`) |>
  unique()

genelist_appg <- readxl::read_excel("data/gene_families/David-Raizen_2014_bio20147500-sup-table_s5.xlsx",
                                 sheet = 1,
                                 range = readxl::cell_cols(2),
                                 col_names = "gene_name") |>
  pull(gene_name) |>
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














# Export ILso genes in families ----

lookup_famid <- dict_panther |>
  pivot_longer(-wbid,
               names_to = "family_long",
               values_to = "present") |>
  filter(present == 1L) |>
  separate_wider_regex(family_long,
                       patterns = c("^[[:print:]]* ",
                                    family_id = "PTHR[0-9]+(?:\\:SF[0-9]+)?",
                                    " WBbt:0000000$")) |>
  inner_join(panther_filt |>
               filter(cell_type == "ILso") |>
               select(family_id) |>
               mutate(family_terms = case_match(
                 family_id,
                 "PTHR47327" ~ "cuticlins/apple domain",
                 "PTHR22907" ~ "ZP domain cuticlins",
                 "PTHR10127" ~ "NAS Zn metalloproteinases",
                 "PTHR24637" ~ "collagens",
                 "PTHR10796" ~ "patched-related")),
             by = "family_id") |>
  mutate(gene_name = i2s(wbid, gids)) |>
  select(gene_name, family_id, family_description = family_terms) |>
  mutate(source = "PANTHER") |>
  bind_rows(tibble(
    gene_name = genelist_hedgehog,
    family_id = "N.A.",
    family_description = "hedgehog-related",
    source = "Hao 2006 Table 1"
  )) |>
  bind_rows(tibble(
    gene_name = genelist_cutl |> setdiff(genelist_cutl_extended) |> i2s(gids),
    family_id = "N.A.",
    family_description = "cuticlin",
    source = "Additional cuticlin by name"
  )) |>
  bind_rows(tibble(
    gene_name = genelist_mam |> i2s(gids),
    family_id = "N.A.",
    family_description = "mam",
    source = "By name and direct paralogy"
  )) |>
  bind_rows(tibble(
    gene_name = genelist_sundaram |> i2s(gids),
    family_id = "N.A.",
    family_description = "Sundaram",
    source = "Sundaram and Pujol 2024 Supp table"
  ))


osc_table_nodup <- osc_table |>
  summarize(osc_amplitude = mean(osc_amplitude),
            .by = c(gene_name, bulk_class))

all_genes |>
  filter(cell_type == "ILso") |>
  inner_join(lookup_famid,
             by = "gene_name") |>
  relocate(family_description, .after = gene_name) |>
  left_join(osc_table_nodup) |> View()
  # writexl::write_xlsx("data/gene_families/250625_ILso_osc_genes_from_defined_families.xlsx")

all_genes |>
  filter(cell_type == "ILso") |>
  inner_join(lookup_famid,
             by = "gene_name") |>
  relocate(family_description, .after = gene_name) |>
  left_join(osc_table_nodup) |>
  filter(shape == "pulsatile", bulk_class == "nonOsc")


# block A
genes_block_A <- c("C10B5.3","dhs-16","C26B9.3","grl-10","let-653","noah-1","fkb-5","atf-8","ets-4",
                 "F32B4.8","fasn-1","phy-2","ugt-58","F46C3.6","F47B7.2","txdc-12.2","noah-2","F52C9.5",
                 "bus-18","ptr-10","dad-1","F59B10.5","H10E21.5","strm-1","egg-6","atf-2","calu-1",
                 "spi-1","ugt-59","acs-3","T10B5.10","tgn-38","peb-1","cutl-6","mlt-11","lpr-6","lpr-5",
                 "lpr-4","lpr-3","Y105E8A.13","nhr-91","Y38H8A.1","Y43F4A.1","cutl-25","rml-5","nex-1",
                 "cuti-1","wrt-10","glam-5","fbn-1")
# block B
genes_block_B <- c("B0361.9","nas-4","col-39","gpx-5","pana-1","cal-5","C34D1.4","gmap-1","C54D10.9",
                 "F10D7.11","epic-1","lips-10","F25B5.3","zipt-2.2","tag-290","spig-16","clec-139",
                 "col-34","sox-3","F41D9.2","cpn-1","F47G4.4","F49E2.5","F53F4.2","F59E11.7","famk-1",
                 "ceeh-1","K04G2.7","vem-1","lag-1","epic-2","tyr-2","M05D6.9","lin-46","cutl-11",
                 "srap-1","madf-4","T23F2.5","kel-8","ztf-6","Y43E12A.2","Y54G2A.11","Y57G11A.4","ztf-29",
                 "sto-4","ZC416.2","ZC449.4","mlt-7")

all_genes |>
  filter(cell_type == "ILso") |>
  inner_join(lookup_famid,
             by = "gene_name") |>
  relocate(family_description, .after = gene_name) |>
  left_join(osc_table_nodup) |>
  filter(gene_name %in% c(genes_block_A, genes_block_B))



genes_in_blocks_and_ilso <- all_genes |>
  filter(cell_type == "ILso") |>
  inner_join(lookup_famid,
             by = "gene_name") |>
  filter(gene_name %in% c(genes_block_A, genes_block_B)) |>
  distinct()

all_genes |>
  filter(gene_name %in% c(genes_block_A, genes_block_B)) |>
  summarize(nb_expressed = n(),
            nb_pulsatile = sum(shape == "pulsatile"),
            nb_non_pulsatile = sum(shape == "nonpulsatile"),
            .by = gene_name) |>
  inner_join(lookup_famid,
             by = "gene_name") |>
  select(gene_name, family_description, nb_pulsatile, nb_non_pulsatile, nb_expressed)
  # flextable::flextable() |> print(preview = "docx")



# Compare sets ----
all_genes |>
  filter(gene_name %in% i2s(genelist, gids),
         cell_type %in% cell_types_osc_both) |>
  summarize(n_peaky = sum(shape == "pulsatile"),
            n_tot = n(),
            `%` = round(100*n_peaky / n_tot),
            .by = cell_type) |>
  arrange(desc(n_peaky))


list(ILso = all_genes |>
       filter(gene_name %in% i2s(genelist, gids),
              cell_type == "ILso") |>
       pull(gene_name),
     hypodermis = all_genes |>
       filter(gene_name %in% i2s(genelist, gids),
              cell_type == "hypodermis") |>
       pull(gene_name)) |>
  eulerr::euler() |>
  plot(quantities = TRUE,
       main = "expressed")


list(ILso = all_genes |>
       filter(gene_name %in% i2s(genelist, gids),
              cell_type == "ILso",
              peaky == "peak") |>
       pull(gene_name),
     hypodermis = all_genes |>
       filter(gene_name %in% i2s(genelist, gids),
              cell_type == "hypodermis",
              peaky == "peak") |>
       pull(gene_name)) |>
  eulerr::euler() |>
  plot(quantities = TRUE,
       main = "pulsatile")


list(
  ILso_expr = all_genes |>
    filter(gene_name %in% i2s(genelist, gids),
           cell_type == "ILso") |>
    pull(gene_name),
  hypodermis_expr = all_genes |>
    filter(gene_name %in% i2s(genelist, gids),
           cell_type == "hypodermis") |>
    pull(gene_name),
  ILso_spiky = all_genes |>
    filter(gene_name %in% i2s(genelist, gids),
           cell_type == "ILso",
           peaky == "peak") |>
    pull(gene_name),
  hypodermis_spiky = all_genes |>
    filter(gene_name %in% i2s(genelist, gids),
           cell_type == "hypodermis",
           peaky == "peak") |>
    pull(gene_name)
)|>
  UpSetR::fromList() |>
  UpSetR::upset()




#~ plot btw cell types ----

# look for genes that are puls in one cell type, nonpuls in another
# this approach based on puls/nonpuls, below with gene expression filtering

genes_differ_btw_ct <- all_genes |>
  filter(cell_type %in% setdiff(cell_types_osc_both, "PHsh")) |>
  summarize(nb_ct_where_puls = sum(shape == "pulsatile"),
            nb_ct_where_expr = n(),
            .by = gene_name) |>
  filter(nb_ct_where_puls > 0,
         nb_ct_where_puls < nb_ct_where_expr)

all_genes |>
  filter(cell_type %in% setdiff(cell_types_osc_both, "PHsh")) |>
  filter(gene_name %in% genes_differ_btw_ct$gene_name) |>
  arrange(gene_name) |>
  filter(str_detect(gene_name, "\\-"))



#~ Genes osc in one, expr in another ----

dir_step1 <- "intermediates/2502/250330_step1/"


gene_expressions <- map_dfr(cell_types_osc_both,
                            ~ {
                              qs::qread(file.path(dir_step1, paste0(.x, "_gene_expressions.qs"))) |>
                                add_column(cell_type = .x,
                                           .before = 1) |>
                                as_tibble()
                            })
stopifnot(all.equal(
  gene_expressions |>
    filter(cell_type %in% cell_types_osc_both) |>
    mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
    filter(prop_cells >= .05) |>
    pull(cellgene) |>
    sort(),
  all_genes |>
    filter(cell_type %in% cell_types_osc_both) |>
    pull(cellgene) |>
    sort()
))
# note: everything in all_genes is above threshold of prop_cells >.05




all_genes |>
  filter(cell_type %in% c("ILso", "seam")) |>
  left_join(gene_expressions,
            by = c("gene_name", "cell_type")) |>
  filter(prop_cells > .2, nb_cells > 50) |>
  summarize(prop_puls = mean(shape == "pulsatile"),
            n_expressed = n(),
            .by = gene_name) |>
  filter(n_expressed >= 2,
         prop_puls >= .2,
         prop_puls <= .8) |>
  filter(str_detect(gene_name, "\\-"))



# check individual genes
ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu.qs")) )


seam_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("seam", "_seu.qs")) )

amphso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("AM_PHso", "_seu.qs")) )

hyp_subseu <- qs::qread( file.path(dir_step1,
                                      paste0("hypodermis", "_seu.qs")) )

goi <- "pana-1"

all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  left_join(gene_expressions,
            by = c("gene_name", "cell_type")) |>
  filter(gene_name == goi) |>
  select(gene_name, cell_type, shape, nb_cells, prop_cells)




Seurat::FeaturePlot(ilso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "ILso")


Seurat::FeaturePlot(seam_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "seam")

Seurat::FeaturePlot(hyp_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "hypodermis")




Seurat::FeaturePlot(amphso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "AM/PHso")










# Families in time ----

source("R/utils_heatmap_processing.R")

# dir_step2 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"
# dir_step4 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250801_step4_analysis/"

all(paste0(cell_types_osc, "_mods_uncentered.qs") %in% list.files(dir_step2))




puls_genes_timed <- map_dfr(cell_types_osc,
        \(ct) {
          
          message(ct)
          
          mods_uncentered <- qs::qread(file.path(dir_step2, paste0(ct, "_mods_uncentered.qs")))
          
          # computed same for all genes
          mean_sf <- lapply(mods_uncentered,
                            \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
            unlist() |>
            log() |>
            mean() |>
            exp()
          
          len <- 128
          
          preds_uncentered <- vapply(mods_uncentered,
                                     \(.mod) predict(.mod,
                                                     type = "response",
                                                     newdata = data.frame(
                                                       pseudotime = (0:(len-1))/len ,
                                                       size_factors = rep(mean_sf, len))
                                     ),
                                     FUN.VALUE = double(len))
          
          
          time_max <- apply(preds_uncentered, 2, which.max) / len
          
          # preds_uncentered1 <- preds_uncentered[,1:20]
          # time_max1 <- time_max[1:20]
          # ngenes <- ncol(preds_uncentered1)
          # printMat::matimage(log1p(preds_uncentered1))
          # points((1:ngenes - 1)/(ngenes - 1), (1 - time_max1), pch = "-", cex = 3.5, col = 'purple')
          
          
          
          # align to bulk phase
          
          peak_times <- left_join(
            enframe(time_max,
                    name = "gene_name",
                    value = "peak_pseudotime"),
            osc_table,
            by = "gene_name"
          ) |>
            filter(bulk_class == "Osc")
          
          alignment <- align_circular(peak_times$bulk_peak,
                                      peak_times$peak_pseudotime*360)
          
          time_max_deg <- if (alignment$invert) {
            ((360 - time_max*360) - alignment$shift) %% 360
          } else {
            (time_max*360 - alignment$shift) %% 360
          }
          
          png(file.path(dir_step4, paste0(ct, "_alignment.png")))
          plot(peak_times$bulk_peak,
               time_max_deg[peak_times$gene_name])
          dev.off()
          
          
          
          
          all_genes |>
            filter(cell_type == ct,
                   shape == "pulsatile") |>
            mutate(time_peak_deg = time_max_deg[gene_name])
          
          
        })


# qs::qsave(puls_genes_timed,
#           file.path(dir_step4, "puls_genes_timed.qs"))


# ordered by bulk (ugly because many pulsatile genes have NA phase in bulk)
puls_genes_timed |>
  mutate(
    gene_name = factor(gene_name,
                       levels = rev(
                         osc_table |>
                           arrange(bulk_peak) |>
                           pull(gene_name) |>
                           intersect(puls_genes_timed$gene_name)
                       ))
  ) |>
  ggplot() +
  theme_classic() +
  scale_x_discrete(breaks = rev(puls_genes_ordered)[(1:20)*5200/20]) +
  coord_flip() +
  geom_point(aes(x = gene_name, y = time_peak_deg),
             alpha = .2, shape = 16)




# pre-order genes
puls_genes_averaged <- puls_genes_timed |>
  mutate(phase_peak = circular::circular(time_peak_deg, units = "degree")) |>
  summarize(mean_phase = circular::mean.circular(phase_peak),
            .by = gene_name) |>
  mutate(mean_phase_deg = (as.numeric(mean_phase) + 360) %% 360 ) |>
  arrange(desc(mean_phase_deg))

puls_genes_ordered <- puls_genes_averaged |>
  pull(gene_name) |>
  fct_inorder() |>
  levels()



puls_genes_timed |>
  mutate(gene_name = factor(gene_name, levels = rev(puls_genes_ordered))) |>
  ggplot() +
  theme_classic() +
  scale_x_discrete(breaks = rev(puls_genes_ordered)[(1:20)*5200/20]) +
  coord_flip() +
  geom_point(aes(x = gene_name, y = time_peak_deg),
             alpha = .2, shape = 16)


left_join(puls_genes_averaged,
          osc_table,
          by = "gene_name") |>
  ggplot() +
  theme_classic() +
  xlab("Phase from bulk") + ylab("Average phase from sc") +
  geom_point(aes(x = bulk_peak, y = mean_phase_deg),
             alpha = .3)


puls_genes_timed |>
  pivot_wider(id_cols = gene_name,
              values_from = time_peak_deg,
              names_from = cell_type) |>
  ggplot() + theme_classic() +
  geom_point(aes(x = ILso, y = AM_PHso))




#By family
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
  mutate(gene_name = i2s(gene_id, gids)) |>
  left_join(puls_genes_averaged,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  ylab(NULL) + xlab("Phase of peak (avg btw cell types)") +
  scale_x_continuous(limits = c(0,360)) +
  geom_point(aes(x = mean_phase_deg, y = family, color = family),
             alpha = .6, show.legend = FALSE)
  # geom_density(aes(x = mean_phase_deg, fill = family),
  #              alpha = .3)


# Timings from Sundaram and Pujol

bind_rows(
  tibble(family = "precuticule",
         gene_name = c("lpr-3", "noah-1", "sym-1", "noah-2", "fbn-1")),
  tibble(family = "furrow col",
         gene_name = c("dpy-2", "dpy-3", "dpy-7", "dpy-8", "dpy-9", "dpy-10" )),
  tibble(family = "annuli col",
         gene_name = c("sqt-3", "dpy-4", "dpy-5", "dpy-13"))
) |>
  left_join(puls_genes_averaged,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  ylab(NULL) + xlab("Phase of peak (avg btw cell types)") +
  scale_x_continuous(limits = c(0,360)) +
  geom_point(aes(x = mean_phase_deg, y = family, color = family),
             alpha = .6, show.legend = FALSE)



bind_rows(
  tibble(family = "precuticule",
         gene_name = c("lpr-3", "noah-1", "sym-1", "noah-2", "fbn-1")),
  tibble(family = "furrow col",
         gene_name = c("dpy-2", "dpy-3", "dpy-7", "dpy-8", "dpy-9", "dpy-10" )),
  tibble(family = "annuli col",
         gene_name = c("sqt-3", "dpy-4", "dpy-5", "dpy-13"))
) |>
  left_join(puls_genes_averaged,
            by = "gene_name") |>
  ggplot() +
  theme_classic() +
  ylab(NULL) + xlab("Phase of peak (avg btw cell types)") +
  scale_x_continuous(limits = c(0,360)) +
  coord_polar() +
  geom_segment(aes(x = mean_phase_deg, y = 0, yend = 1, color = family),
             alpha = .6, show.legend = FALSE)





#~ heatmap ----

ct
all_tests_bin <- qs::qread(file.path(dir_step4, "bins", paste0(ct, ".qs")))

signif_fams <- all_tests_bin |>
  filter(FDR < 0.05, observed > 0) |>
  summarize(family = list(family_id),
            .by = time_bin) |>
  deframe() |> unlist() |> unique()

if(length(signif_fams) < 2) next

fam_by_time <- all_tests_bin |>
  filter(family_id %in% signif_fams) |>
  arrange(time_bin) |>
  mutate(signif = -log10(FDR)) |>
  # mutate(signif = odds_ratio) |>
  pivot_wider(id_cols = family_id,
              names_from = time_bin,
              values_from = signif) |>
  column_to_rownames("family_id") |>
  as.matrix()


pt_colnames <- round(bins_start + (bin_width/2), 1)
pt_colnames[2 * (1:(length(pt_colnames)/2))] <- ""
colnames(fam_by_time) <- pt_colnames


pheatmap::pheatmap(fam_by_time,
                   cluster_rows = TRUE,
                   clustering_distance_rows = "correlation",
                   cluster_cols = FALSE,
                   scale = "none")

n_tfs <- nrow(tf_by_time)

# height 2.5 for ILso/main figure, for supp,  2.7 mm (.1 in) per gene + 6 mm (.24 in) for legend
# divide by 2 in figure
pheatmap::pheatmap(
  tf_by_time,
  color = colorRampPalette(c("white", "#C994C7", "#DD1C77"))(100),
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  filename = paste0(dir_out, "/heatmaps_timebins/", ct,".pdf"),
  fontsize = 10,
  width = 9, height = .5+n_tfs*.15
)















# Illustrate GAM ----



dir_step2 <- "intermediates/2502/250624_step2/"
dir_step1 <- "intermediates/2502/250609_step1/"

# if working from external HDD
# dir_step2 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250624_step2/"
# dir_step1 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step1/"

dir_fig_gam <- "presentations/figures/250625_gam_illustrations"
# dir.create(dir_fig_gam)



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

# goi <- "grl-18"
# goi <- "nhr-23"
# goi <- "col-109"
# goi <- "pugs-11"
# 
# 
# goi <- "rps-27A"
# goi <- "dnj-1"

# for(goi in c("grl-18", "nhr-23", "col-109", "pugs-11", "rps-27A", "dnj-1")){



#~| cells ----

dat <- FetchData(ilso_subseu, vars = c("PC_1","PC_2",goi))

# for ILso, invert axes for easier interpretation
dat$PC_1 <- -dat$PC_1



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
                           size = 1.5,
                           raster.dpi = 500)



ggsave(paste0(goi, "_expr.pdf"),
       path = dir_fig_gam,
       width = 55, height = 55, units = "mm",
       scale = 1)






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

ggsave(paste0(goi, "_devexpl.pdf"),
       path = dir_fig_gam,
       width = 50, height = 45, units = "mm",
       scale = 1)

# }




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





#~ gene ----

goi <- "lgc-34"


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
                           size = 1.5,
                           raster.dpi = 500)



ggsave(paste0(goi, "_bwm_expr.pdf"),
       path = dir_fig_gam,
       width = 55, height = 55, units = "mm")




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


ggsave(paste0(goi, "_bwm_devexpl.pdf"),
       path = dir_fig_gam,
       width = 50, height = 45, units = "mm")










