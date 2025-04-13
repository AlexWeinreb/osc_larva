

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



out_dir <- "presentations/2503_figures/"


dir_assembled <- "intermediates/2502/250328_assembled/"
dir_step3 <- "intermediates/2502/250330_step3_genes_by_celltype//"

dir_step1 <- "intermediates/2502/250330_step1/"





# Clusters and whether they are phasic according to bulk
bulk_cat <- qs::qread(file.path(dir_assembled, "250330_coherence_by_ct.qs")) |>
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

all_genes <- qs::qread(file.path(dir_step3, "all_genes.qs"))

stopifnot(all( unique(cell_types_info$cell_type) |> sort() %in%
                 unique(all_genes$cell_type) |> sort() ))








# General plots ----


# Volcano
all_genes |>
  slice_sample(n = 500) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = peaky),
             alpha = .2)




#~ Plot nb of DEGs/peaky ----

# # not equal anymore
# stopifnot(all.equal(
#   sort(unique(bulk_cat$cell_type)),
#   sort(unique(all_genes$cell_type))
# ))



gg_peaky_by_celltype <- all_genes |>
  summarize(nb_selected = sum(peaky == "peak"),
            nb_tested = n(),
            .by = cell_type) |>
  mutate(prop_peaky = nb_selected / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  mutate(cell_type_noneur = if_else(tissue == "neuron", "", cell_type)) |>
  ggplot() +
  theme_classic() +
  xlab("Number of genes tested") +
  ylab("Proportion peaky") +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0,.8)) +
  scale_x_continuous(labels = scales::label_comma(),
                     limits = c(0,NA)) +
  scale_shape_manual(values = c(19,8)) +
  scale_size_manual(values = c(1, 2.5)) +
  aes(x = nb_tested,
      y = prop_peaky,
      color = tissue) +
  geom_point(aes(shape = bulk, size = bulk)) +
  ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)

gg_peaky_by_celltype

# ggsave("peaky_genes_by_celltype.pdf", gg_peaky_by_celltype,
#        path = out_dir,
#        width = 20, height = 15, units = "cm")



ordered_cells <- all_genes |>
  summarize(nb_peaky = sum(peaky == "peak"),
            nb_tested = n(),
            .by = cell_type) |>
  mutate(prop_peaky = nb_peaky / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  arrange(tissue, desc(prop_peaky)) |>
  pull(cell_type) |> fct_inorder() |> levels()

all_genes |>
  summarize(nb_peaky = sum(peaky == "peak"),
            nb_not_peaky = sum(peaky == "non-peak"),
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
  scale_fill_manual(values = c(peaky = "red4", not_peaky = "grey70")) +
  xlab(NULL) + ylab("Number of genes") +
  geom_col(aes(x = cell_type, y = count, fill = type),
           show.legend = FALSE) +
  geom_text(aes(x = cell_type,
                y = nb_peaky,
                label = prop_peaky),
            data = all_genes |>
              summarize(nb_peaky = sum(peaky == "peak"),
                        nb_tested = n(),
                        .by = cell_type) |>
              mutate(prop_peaky = round( 100 * nb_peaky / nb_tested) |> paste0("%") ),
            nudge_y = 100
            )



scale_y_continuous(labels = scales::label_percent(),
                   limits = c(0,.8)) +
scale_x_continuous(labels = scales::label_comma(),
                   limits = c(0,NA)) +
scale_shape_manual(values = c(19,8)) +
scale_size_manual(values = c(1, 2.5)) +
aes(x = nb_tested,
    y = prop_peaky,
    color = tissue) +
geom_point(aes(shape = bulk, size = bulk)) +
ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)





# Compare clusters to bulk ----



cell_types_cat <- left_join(cell_types_info, bulk_cat,
                            by = "cell_type") |>
  mutate(sc_heatmap = case_when(
    similarity_diag > .3 & p_adj < 0.05 ~ "phasic",
    similarity_diag <= .3 | p_adj >= .05 ~ "nonphasic"
  ))

table(bulk = cell_types_cat$bulk, sc_heatmap = cell_types_cat$sc_heatmap)

list(bulk = cell_types_cat$cell_type[cell_types_cat$bulk == "phasic"],
     sc_heatmap = cell_types_cat$cell_type[cell_types_cat$sc_heatmap == "phasic"]) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

cell_types_osc_both <- cell_types_cat |>
  filter(bulk == "phasic",
         sc_heatmap == "phasic") |>
  pull(cell_type)



# Compare genes to known bulk ----

#~ overlap peaky/Osc ----

all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  select(gene_name, peaky) |>
  summarize(is_peaky = any(peaky == "peak"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(sc = if_else(is_peaky, "peaky", "not peaky")) |>
  (\(df) table(bulk = df$bulk_class, sc = df$sc))()



all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  select(gene_name, peaky) |>
  summarize(is_peaky = any(peaky == "peak"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(`single-cell` = if_else(is_peaky, "peaky", "not peaky")) |>
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
              values_from = "peaky") |>
  left_join(osc_table |> select(gene_name, bulk_class),
            by = "gene_name") |>
  relocate(bulk_class, .after = gene_name)




# check genes ----




#~ focus on nonOsc clusters ----

all_genes |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "nonphasic") |>
  summarize(nb_tested = n(),
            nb_peaky = sum(peaky == "peak"),
            .by = c(cell_type, tissue)) |>
  mutate(prop_peaky = round(100 * nb_peaky / nb_tested)) |>
  relocate(tissue, .before = 1) |>
  arrange(desc(prop_peaky)) #|>  clipr::write_clip()





#~ focus on osc clusters ----

all_genes |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "phasic") |>
  summarize(nb_tested = n(),
            nb_peaky = sum(peaky == "peak"),
            .by = c(cell_type, tissue)) |>
  mutate(prop_peaky = round(100 * nb_peaky / nb_tested)) |>
  relocate(tissue, .before = 1) |>
  arrange(tissue, cell_type) #|>  clipr::write_clip()




#~ intersections ----


#~ proportion of tested genes ----

degs_of_interest <- all_genes |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "phasic") |>
  mutate(is_peaky = peaky == "peak")



# proportion of testable genes (high genes) that are peaky in both
intersections <- expand_grid(
  cell_1 = unique(degs_of_interest$cell_type),
  cell_2 = unique(degs_of_interest$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(degs_of_interest$gene_name[degs_of_interest$cell_type == .x &
                                                     degs_of_interest$is_peaky],
                        degs_of_interest$gene_name[degs_of_interest$cell_type == .y &
                                                     degs_of_interest$is_peaky]))),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(degs_of_interest$gene_name[degs_of_interest$cell_type == .x],
                    degs_of_interest$gene_name[degs_of_interest$cell_type == .y])))
  )

mat_intersections <- intersections |>
  mutate(prop_inter = 100 * intersection / union) |>
  select(cell_1, cell_2, prop_inter) |>
  pivot_wider(id_cols = cell_1,
              names_from = "cell_2",
              values_from = "prop_inter") |>
  column_to_rownames("cell_1") |>
  as.matrix()

# opar <- par(no.readonly = TRUE)
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
  cell_1 = unique(degs_of_interest$cell_type),
  cell_2 = unique(degs_of_interest$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(degs_of_interest$gene_name[degs_of_interest$cell_type == .x &
                                                     degs_of_interest$is_peaky],
                        degs_of_interest$gene_name[degs_of_interest$cell_type == .y &
                                                     degs_of_interest$is_peaky]))),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(degs_of_interest$gene_name[degs_of_interest$cell_type == .x &
                                                 degs_of_interest$is_peaky],
                    degs_of_interest$gene_name[degs_of_interest$cell_type == .y &
                                                 degs_of_interest$is_peaky])))
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
mat_intersections <- mat_intersections[!is.na(rowSums(mat_intersections)),
                                       !is.na(colSums(mat_intersections))]
diag(mat_intersections) <- NA_real_


# opar <- par(no.readonly = TRUE)
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


degs_of_interest <- all_genes |>
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
    
    gene_list <- ct_genes$gene_id[ct_genes$peaky == "peak"] |> intersect(dict$wbid)
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

dict_panther <- qs::qread(file.path(dir_step3, "panther_dict.qs"))

nb_genes_in_fam <- dict_panther[,-1] |> colSums()

dict_panther <- dict_panther[ , c(TRUE, nb_genes_in_fam > 5) ]


all_panther <- map_dfr(
  cell_types_osc_both,
  \(ct){
    
    ct_genes <- all_genes |>
      filter(cell_type == ct) |>
      mutate(gene_id = s2i(gene_name, gids))
    
    gene_list <- ct_genes$gene_id[ct_genes$peaky == "peak"] |> intersect(dict_panther$wbid)
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

panther_filt |>
  count(family_id) |>
  arrange(desc(n))

# check the families
i <- 0

i <- i+1
panther_filt[i,]
col_nb <- str_detect(colnames(dict_panther), paste0(panther_filt$family_id[i], " WBbt:0000000"))
colnames(dict_panther)[col_nb]


genelist <- dict_panther$wbid[dict_panther[col_nb] == 1L]
genelist |> i2s(gids)


dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther), paste0("PTHR10127", " WBbt:0000000"))] == 1L  ] |>
  i2s(gids)




#~ manually defined families ----


genelist_cutl_extended <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther), paste0("PTHR47327", " WBbt:0000000"))] == 1L |
                                              dict_panther[str_detect(colnames(dict_panther), paste0("PTHR22907", " WBbt:0000000"))] == 1L ]

genelist_cutl <- gids$gene_id[which(startsWith(gids$symbol, "cutl-"))]

genelist_nas_metalloproteases <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther),
                                                                           paste0("PTHR10127", " WBbt:0000000"))] == 1L  ]

genelist_hedgehog <- readxl::read_excel("data/gene_families/hedgehog_hao2006_table1.xlsx", skip = 1) |>
  pull(gene_name)

genelist_collagen_by_name <- gids$gene_id[which(startsWith(gids$symbol, "col-"))]
genelist_collagen_by_panther <- dict_panther$wbid[dict_panther[str_detect(colnames(dict_panther),
                                                                          paste0("PTHR24637", " WBbt:0000000"))] == 1L  ]

list(name = genelist_collagen_by_name,
     panther = genelist_collagen_by_panther) |>
  eulerr::euler() |>
  plot(quantities = TRUE)
genelist_collagen_by_name |> setdiff(genelist_collagen_by_panther) |> i2s(gids) |> paste(collapse = ", ")
genelist_collagen_by_panther |> setdiff(genelist_collagen_by_name) |> i2s(gids) |> paste(collapse = ", ")


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
  summarize(n_peaky = sum(peaky == "peak"),
            n_tot = n(),
            `%` = round(100*n_peaky / n_tot),
            .by = cell_type) |>
  arrange(desc(n_peaky))




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

# all_genes |>
#   filter(cell_type == "ILso") |>
#   inner_join(lookup_famid,
#             by = "gene_name") |>
#   relocate(family_description, .after = gene_name) |>
#   left_join(osc_table_nodup) |>
#   writexl::write_xlsx("data/gene_families/250413_ILso_osc_genes_from_defined_families.xlsx")

all_genes |>
  filter(cell_type == "ILso") |>
  inner_join(lookup_famid,
            by = "gene_name") |>
  relocate(family_description, .after = gene_name) |>
  left_join(osc_table_nodup) |>
  filter(peaky == "peak", bulk_class == "nonOsc")



#~ plot btw cell types ----
genelist |> i2s(gids) |> head()


all_genes |>
  filter(gene_name %in% i2s(genelist, gids),
         cell_type %in% cell_types_osc_both) |>
  summarize(n_peaky = sum(peaky == "peak"),
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



# Genes osc in one, expr in another ----

gene_expressions <- map_dfr(cell_types_osc_both,
                            ~ {
                              qs::qread(file.path(dir_step1, paste0(.x, "_gene_expressions.qs"))) |>
                                add_column(cell_type = .x,
                                           .before = 1) |>
                                as_tibble()
                            })

all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  left_join(gene_expressions,
            by = c("gene_name", "cell_type")) |>
  filter(prop_cells > .3, nb_cells > 50) |>
  summarize(prop = mean(peaky == "peak"),
            n = n(),
            .by = gene_name) |>
  filter(n > 3,
         prop > .3,
         prop < .7)


all_genes |>
  filter(cell_type %in% cell_types_osc_both) |>
  left_join(gene_expressions,
            by = c("gene_name", "cell_type")) |>
  filter(gene_name == "col-34") |>
  select(gene_name, cell_type, peaky, nb_cells, prop_cells)

# e.g. col-155 in ILso and seam



ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu.qs")) )


seam_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("seam", "_seu.qs")) )

goi <- "col-34"

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










