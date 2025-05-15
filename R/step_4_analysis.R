

# Inits ----
library(tidyverse)
# library(Seurat)

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



out_dir <- "presentations/2503_figures/"


dir_assembled <- "intermediates/2502/250508_assembled/"


dir_step3 <- "intermediates/2502/250425_step3_genes_by_celltype/"
dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"
# dir_step2_noncent <- "intermediates/2502/250422_step2_nb/"




# Clusters and whether they are phasic according to bulk
bulk_cat <- qs::qread(file.path(dir_assembled, "250508_coherence_by_ct.qs")) |>
  mutate(bulk = case_when(
    p_adj < .05 & mean_coherence > 0.15 ~ "phasic",
    p_adj >= .05 | mean_coherence <= 0.15 ~ "nonphasic"
  )) |>
  select(cell_type, tissue, bulk)



osc_raw <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID),
         gene_name = i2s(gene_id, gids) )


osc_table <- osc_raw |>
  select(gene_name, bulk_class = Class, osc_amplitude = OscAmplitude)



# Load step 3 ----

cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name))



stopifnot(all.equal(
  cell_types_info$cell_type |> unique() |> sort(),
  all_genes$cell_type |> unique() |> sort()
))






# General plots ----




#~ Plot nb of puls ----



gg_puls_by_celltype <- all_genes |>
  summarize(nb_selected = sum(shape == "pulsatile"),
            nb_tested = n(),
            .by = cell_type) |>
  mutate(prop_peaky = nb_selected / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  mutate(cell_type_noneur = if_else(tissue == "neuron", "", cell_type)) |>
  ggplot() +
  theme_classic() +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0,.4)) +
  scale_x_continuous(labels = scales::label_comma(),
                     limits = c(0,NA)) +
  scale_shape_manual(values = c(19,8)) +
  scale_size_manual(values = c(1, 2.5)) +
  aes(x = nb_tested,
      y = prop_peaky,
      color = tissue) +
  xlab("Number of genes tested") +
  ylab("Proportion pulsatile") +
  geom_point(aes(shape = bulk, size = bulk)) +
  ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)

gg_puls_by_celltype

# ggsave("peaky_genes_by_celltype.pdf", gg_peaky_by_celltype,
#        path = out_dir,
#        width = 20, height = 15, units = "cm")



all_genes |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(prop_puls)

ordered_cells <- all_genes |>
  summarize(nb_puls = sum(shape == "pulsatile"),
            nb_tested = n(),
            .by = cell_type) |>
  mutate(prop_puls = nb_puls / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  arrange(tissue, desc(prop_puls)) |>
  pull(cell_type) |> fct_inorder() |> levels()

all_genes |>
  summarize(nb_puls = sum(shape == "pulsatile"),
            nb_nonpuls = sum(shape == "nonpulsatile"),
            .by = cell_type) |>
  left_join(bulk_cat, by = "cell_type") |>
  pivot_longer(cols = starts_with("nb_"),
               names_to = "type",
               values_to = "count",
               names_prefix = "nb_") |>
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



cell_types_info |>
  mutate(signif = case_when(
    p_diag_adj < .05 & p_coherence_adj < .05 ~ "both",
    p_diag_adj < .05 & p_coherence_adj >= .05 ~ "diag only",
    p_diag_adj >= .05 & p_coherence_adj < .05 ~ "bulk only",
    p_diag_adj >= .05 & p_coherence_adj >= .05 ~ "neither"
  )) |>
  ggplot() +
  theme_classic() +
  xlab("Mean local phase coherence (bulk)") +
  ylab("Diagonal similarity (sc)") +
  scale_shape_manual(values = c("both" = 8, "diag only" = 7, "bulk only" = 9, "neither" = 16)) +
  geom_point(aes(x = mean_coherence, y = similarity_diag, color = tissue,
                 shape = signif),
             size = 3) +
  ggrepel::geom_text_repel(aes(x = mean_coherence, y = similarity_diag, label = cell_type))


cell_types_osc_both <- cell_types_info |>
  filter(p_coherence_adj < .05,
         p_diag_adj < .05) |>
  pull(cell_type)










# Compare genes to known bulk ----

#~ overlap peaky/Osc ----

all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  select(gene_name, shape) |>
  summarize(is_puls = any(shape == "pulsatile"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(sc = if_else(is_puls, "pulsatile", "nonpulsatile")) |>
  (\(df) table(bulk = df$bulk_class, sc = df$sc))()



all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  select(gene_name, shape) |>
  summarize(is_puls = any(shape == "pulsatile"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(`single-cell` = if_else(is_puls, "pulsatile", "nonpulsatile")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(.8,.7)) +
  xlab("Osc amplitude (bulk)") +
  geom_density(aes(x = osc_amplitude, fill = `single-cell`),
               alpha = .5)
# 650 x 400





#~ abu genes in pha muscle ----

appg_genes <- c("abu-5", "abu-2", "pqn-76", "pqn-78", "pqn-79", "Y5H2A.4", "pqn-91",
                "pqn-90", "pqn-2", "abu-3", "abu-1", "abu-4", "F07H5.8", "pqn-71",
                "pqn-16", "M02G9.1", "abu-6", "abu-7", "abu-15", "abu-8",
                "M02G9.2", "pqn-54", "pqn-57", "abu-9", "M02G9.3", "abu-14",
                "abu-11", "pqn-13", "F35A5.4")

appg_genes_strict <- c("abu-11", "pqn-54", "pqn-2","abu-15","abu-1","abu-7","abu-8",
                       "abu-6","abu-14","abu-4","pqn-57", "pqn-71","pqn-13")


length(appg_genes)
length(appg_genes_strict)


all_genes |>
  filter(cell_type %in% c("pharyngeal_muscle", "pharynx_epithelial")) |>
  filter(gene_name %in% appg_genes_strict) |>
  pivot_wider(id_cols = gene_name,
              names_from = "cell_type",
              values_from = "shape") |>
  left_join(osc_table |> select(gene_name, bulk_class),
            by = "gene_name") |>
  relocate(bulk_class, .after = gene_name)




# check genes ----




#~ focus on nonOsc clusters ----

all_genes |>
  filter(! cell_type %in% cell_types_osc_both) |>
  summarize(nb_tested = n(),
            nb_puls = sum(shape == "pulsatile"),
            .by = cell_type) |>
  mutate(prop_peaky = round(100 * nb_puls / nb_tested)) |>
  arrange(desc(prop_peaky)) #|>  clipr::write_clip()





#~ focus on osc clusters ----

all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  summarize(nb_tested = n(),
            nb_puls = sum(shape == "pulsatile"),
            .by = cell_type) |>
  mutate(prop_peaky = round(100 * nb_puls / nb_tested)) |>
  arrange(cell_type) #|>  clipr::write_clip()




#~ intersections ----


#~ proportion of tested genes ----

in_osc_ct <- all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  mutate(is_puls = (shape == "pulsatile"))



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
  hclust(method = "ward.D") |>
  as.dendrogram() |>
  reorder(wts = diag(mat_intersections)) |>
  rev() |>
  as.hclust()


pheatmap::pheatmap(mat_intersections,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # width = 6, height = 5,
                   # filename = file.path(out_dir, "intersections.png"),
                   annotation_row = bulk_cat |> select(-bulk) |> column_to_rownames("cell_type"),
                   annotation_col = bulk_cat |> select(-bulk)  |> column_to_rownames("cell_type"),
                   annotation_colors = list(tissue = c(
                     glia = scales::hue_pal()(8)[[1]],
                     muscle = scales::hue_pal()(8)[[2]],
                     neuron = scales::hue_pal()(8)[[3]],
                     other = scales::hue_pal()(8)[[4]],
                     pharynx = scales::hue_pal()(8)[[5]],
                     reproductive = scales::hue_pal()(8)[[6]],
                     skin = scales::hue_pal()(8)[[7]]
                   )))


# dev.off()






#~ proportion of peaky genes in both (out of peaky genes in each) ----

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
                   # filename = file.path(out_dir, "intersections.png"),
                   annotation_row = bulk_cat |> select(-bulk) |> column_to_rownames("cell_type"),
                   annotation_col = bulk_cat |> select(-bulk)  |> column_to_rownames("cell_type"),
                   annotation_colors = list(tissue = c(
                     glia = scales::hue_pal()(8)[[1]],
                     muscle = scales::hue_pal()(8)[[2]],
                     neuron = scales::hue_pal()(8)[[3]],
                     other = scales::hue_pal()(8)[[4]],
                     pharynx = scales::hue_pal()(8)[[5]],
                     reproductive = scales::hue_pal()(8)[[6]],
                     skin = scales::hue_pal()(8)[[7]]
                   )))


# dev.off()


# # Check Hongyan's list of lysosomal genes in seams
# 
# genelist <- clipr::read_clip() |> unique()
# length(genelist)
# 
# head(i2s(genelist, gids))
# 
# xx <- all_res_DE |>
#   filter(cell_type == "seam") |>
#   left_join(osc_raw, by = c(gene_id = "WB_ID")) |>
#   filter(gene_id %in% genelist)
# 
# table(bulk = xx$Class == "Osc",
#             sc = xx$is_peaky)
# 
# xx |>
#   filter(Class == "Osc",
#          is_peaky) |> pull(gene_name) |> clipr::write_clip()






# Gene families ----


in_osc_ct <- all_genes |>
  filter(cell_type %in% cell_types_osc_both)


#~ GO ----

dict <- wormbaseEnrich::fetch_dictionary("go")

i <- 1
ct <- cell_types_osc_both[[i]]

all_go <- cell_types_osc_both |>
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


wormbaseEnrich::plot_enrichment_results(enr_go,labels_inside_columns = FALSE)

enr_go |>
  filter(observed > 5)





wormbaseEnrich::plot_enrichment_results(enr_go |>
                                          filter(observed > 5),labels_inside_columns = FALSE)







#~ PANTHER families ----

dict_panther <- qs::qread("data/gene_families/panther_dict.qs")

nb_genes_in_fam <- dict_panther[,-1] |> colSums()

dict_panther <- dict_panther[ , c(TRUE, nb_genes_in_fam > 5) ]


all_panther <- map_dfr(
  cell_types_osc_both,
  \(ct){
    
    ct_genes <- all_genes |>
      filter(cell_type == ct) |>
      mutate(gene_id = s2i(gene_name, gids))
    
    gene_list <- ct_genes$gene_id[ct_genes$shape == "pulsatile"] |> intersect(dict_panther$wbid)
    background_genes <- ct_genes$gene_id |> intersect(dict_panther$wbid)
    
    
    enr_res <- wormbaseEnrich::enrichment_analysis(
      gene_list = gene_list,
      dictionary = dict_panther,
      background_genes = background_genes
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
  })




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

panther_filt <- left_join(panther_filt,
                          manual_annot_panther_families,
                          by = c(family_id = "family")) |>
  select(cell_type, description, family_id, expected, observed, enrichment_fc, FDR)


# since both family and subfamily, only keep main family
subfams <- str_match(panther_filt$family_id, "^(PTHR[0-9]+)\\:SF[0-9]+$")
stopifnot(all(
  na.omit(subfams[,2]) %in% panther_filt$family_id
))

panther_filt <- panther_filt |> filter(! family_id %in% subfams[,1])


#~~ res ----

panther_filt |>
  count(description) |>
  arrange(desc(n))

panther_filt |>
  count(cell_type) |>
  arrange(desc(n))


panther_filt |>
  arrange(cell_type, FDR) |>
  mutate(description = fct_inorder(description)) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  facet_wrap(~ cell_type, scales = "free_y") +
  geom_col(aes(x = description,
               y = -log10(FDR))) +
  geom_text(aes(x = description, y = 0.1, label = description), 
            hjust = 0, color = "grey") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


map(unique(panther_filt$cell_type),
    ~ {
      panther_filt |>
        filter(cell_type == .x) |>
        arrange(desc(FDR)) |>
        mutate(description = fct_inorder(description)) |>
        ggplot() +
        theme_classic() +
        coord_flip() +
        geom_col(aes(x = description,
                     y = -log10(FDR))) +
        shadowtext::geom_shadowtext(aes(x = description, y = 0.1, label = description), 
                                    hjust = 0) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        ggtitle(NULL, subtitle = .x)
    }
) |> patchwork::wrap_plots()









#~ manually defined families ----


genelist_cutl_extended <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther), paste0("PTHR47327", " WBbt:0000000"))] == 1L |
                                              dict_panther[str_detect(colnames(dict_panther), paste0("PTHR22907", " WBbt:0000000"))] == 1L ]

genelist_cutl <- gids$gene_id[which(startsWith(gids$symbol, "cutl-"))]

genelist_nas_metalloproteases <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther),
                                                                           paste0("PTHR10127", " WBbt:0000000"))] == 1L  ]

genelist_hedgehog <- readxl::read_excel("data/gene_families/hedgehog_hao2006_table1.xlsx", skip = 1) |>
  pull(gene_name)

genelist_collagen_by_panther <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther),
                                                                          paste0("PTHR24637", " WBbt:0000000"))] == 1L  ]


genelist_patched <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther),
                                                              paste0("PTHR10796", " WBbt:0000000"))] == 1L  ]

genelist_abu <- str_match(read_lines("data/gene_families/panther_abu.txt"),
                          "WormBase=(WBGene[0-9]{8})")[,2] |>
  na.omit()

genelist_mam <- gids$name[which(startsWith(gids$symbol, "mam-"))] |>
  c("mlt-9", "CD4.11","R04B3.1", "R04B3.3", "Y43D4A.5") |>
  s2i(gids)



genelist_sundaram <- readxl::read_excel("data/gene_families/sundaram.xlsx",
                                        skip = 1L) |>
  pull(gene_id)







genelist <- genelist_abu
genelist <- appg_genes |> wb_clean_gene_names()
genelist <- s2i(genelist_hedgehog, gids)
genelist <- genelist_collagen_by_panther
genelist <- genelist_cutl_extended
genelist <- genelist_patched
genelist <-  genelist_nas_metalloproteases
genelist <- genelist_mam
genelist <- genelist_sundaram



if(any(is.na(genelist))) warning("NA")


length(genelist)

i2s(genelist, gids) |>
  paste0(collapse = ", ") |>
  message()

all_genes |>
  filter(gene_name %in% i2s(genelist, gids),
         cell_type %in% cell_types_osc_both) |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            `%` = round(100*n_puls / n_tot),
            .by = cell_type) |>
  arrange(desc(n_puls))




# Export ILso genes in families

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
  # writexl::write_xlsx("data/gene_families/250429_ILso_osc_genes_from_defined_families.xlsx")

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




# Illustrate GAM ----
dir_step2 <- "intermediates/2502/250512_step2/"




smooth_preds <- list.files(dir_step2,
                           pattern = "_preds\\.qs$") |>
  enframe(value = "filename",
          name = NULL) |>
  separate_wider_regex(filename,
                       patterns = c(
                         cell_type = "^.+",
                         "_preds\\.qs"
                       ),
                       cols_remove = FALSE) |>
  pmap(\(cell_type, filename){
    mat_preds <- qs::qread(file.path(dir_step2,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)



matplot(smooth_preds[,1:200], type = "l")


mods_uncentered <- qs::qread(file.path(dir_step2, paste0("ILso", "_mods_uncentered.qs")))


plot(mods_uncentered[[goi]])

goi <- "pugs-11"
mod <- mods_uncentered[[goi]]

ggplot() +
  theme_classic() +
  ylab("Expression (SCT)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod$model$pseudotime,
                           count = mod$y),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod$model$pseudotime,
                          prediction = mod$fitted.values),
            color = 'red3') +
  ggtitle(goi)

# ggsave(paste0("gam_ILso_", goi, ".pdf"),
#        path = "presentations/250514_gams",
#        width = 30, height = 30, units = "mm",
#        scale = 5)





dat <- data.frame(
  expr = mat_sct[goi,],
  pseudotime = pseudotime/max(pseudotime)
)


mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 4, bs = 'cc'),
          data = mod$model,
          family = gaussian())



ggplot() +
  theme_classic() +
  ylab("Expression (SCT)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = mod2$y),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod2$model$pseudotime,
                          prediction = mod2$fitted.values),
            color = 'red3') +
  ggtitle(goi)






mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 6, bs = 'cc'),
          data = data.frame(
            expr = mat_sct[.gene,],
            pseudotime = pseudotime/max(pseudotime)
          ),
          family = gaussian())


ggplot() +
  theme_classic() +
  ylab("Expression (SCT)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = mod2$y),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod2$model$pseudotime,
                          prediction =  mod2$fitted.values ),
            color = 'red3') +
  ggtitle(.gene)

# NB-GAM
mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 6, bs = 'cc'),
                  data = data.frame(
                    expr = GetAssayData(seu, assay = "SCT", layer = "data")[.gene,],
                    pseudotime = pseudotime/max(pseudotime)
                  ),
                  family = mgcv::nb(link = "log"))

mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 6, bs = 'cc'),
                  data = data.frame(
                    expr = GetAssayData(seu, assay = "RNA", layer = "counts")[.gene,],
                    pseudotime = pseudotime/max(pseudotime)
                  ),
                  family = mgcv::nb(link = "log"))



ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = log10(1 + mod2$model$expr  )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod2$model$pseudotime,
                          prediction = log10( 1 + mod2$fitted.values  )),
            color = 'red3') +
  ggtitle(.gene)


goi <- names(mods_uncentered)[[1]]
# with size factors
nf <- edgeR::calcNormFactors(mat_cnt)
size_factors <- colSums(mat_cnt) * nf

# goi<- sample(high_genes,1)
# goi

mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 6, bs = 'cc') + offset(log(size_factors)),
          data = data.frame(
            expr = mat_cnt[goi,],
            pseudotime = pseudotime/max(pseudotime)
          ),
          family = mgcv::nb(link = "log"))



ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = log10(1 + mod2$model$expr / exp(mod2$model$`offset(log(size_factors))`) )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod2$model$pseudotime,
                          prediction = log10( 1 + mod2$fitted.values / exp(mod2$model$`offset(log(size_factors))`) )),
            color = 'red3') +
  ggtitle(goi)




ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = log10(1 + mod2$model$expr / exp(mod2$model$`offset(log(size_factors))`) )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(
              pseudotime = (0:(len-1))/len,
              prediction = log10(1 + preds_uncentered[,goi] / mean_sf)
            ),
            color = 'red3')







# offset as argument
mod2 <- mgcv::gam(expr ~ s(pseudotime, k = 6, bs = 'cc'),
                  offset = offset(log(size_factors)),
                  data = data.frame(
                    expr = mat_cnt[goi,],
                    pseudotime = pseudotime/max(pseudotime)
                  ),
                  family = mgcv::nb(link = "log"))



ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = log10(1 + mod2$model$expr / exp(mod2$model$`(offset)`) )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod2$model$pseudotime,
                          prediction = log10( 1 + mod2$fitted.values / exp(mod2$model$`(offset)`) )),
            color = 'red3')



ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod2$model$pseudotime,
                           count = log10(1 + mod2$model$expr / exp(mod2$model$`(offset)`) )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(
              pseudotime = (0:(len-1))/len,
              prediction = preds_uncentered[,goi] / mean(size_factors)
            ),
            color = 'red3')




goi <- names(mods_uncentered)[[1]]













