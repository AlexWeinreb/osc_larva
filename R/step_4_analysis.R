

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



dir_step3 <- "intermediates/2502/250325_step3_genes_by_celltype/"



osc_raw <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA")


# Clusters and whether they are phasic according to bulk
bulk_cat <- qs::qread("intermediates/240620_reanalysis_panglial_2021/241105_mean_dotprod_by_celltype_appr3_100k.qs") |>
  as_tibble() |>
  mutate(category = case_when(
    mean_mean_dotprod > 0.2 ~ "phasic",
    mean_mean_dotprod < 0.1 | p_adj > .05 ~ "nonphasic",
    .default = "intermediate"
  )) |>
  mutate(cell_type = str_remove(cell_type, "\\?$")) |>
  mutate(bulk = if_else(p_adj < 0.05, "oscillating", "non oscillating")) |>
  select(cell_type, tissue, bulk)






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

# not equal anymore
stopifnot(all.equal(
  sort(unique(bulk_cat$cell_type)),
  sort(unique(all_res_DE$cell_type))
))



gg_DEGs_by_celltype <- all_res_DE |>
  summarize(nb_signif = sum(p_adj < 0.05, na.rm = TRUE),
            nb_tested = n(),
            .by = cell_type) |>
  mutate(prop_de = nb_signif / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  mutate(cell_type_noneur = if_else(tissue == "neuron", "", cell_type)) |>
  ggplot((aes(x = nb_tested, y = prop_de, color = tissue))) +
  theme_classic() +
  xlab("Number of genes tested") + ylab("Proportion DEGs") +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::label_comma(),
                     limits = c(0,NA)) +
  scale_shape_manual(values = c(19,8)) +
  scale_size_manual(values = c(1, 2.5)) +
  geom_point(aes(shape = bulk, size = bulk)) +
  ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)

gg_DEGs_by_celltype

# ggsave("DEGs_by_celltype.pdf", gg_DEGs_by_celltype,
#        path = out_dir,
#        width = 20, height = 15, units = "cm")



gg_peaky_by_celltype <- all_res_DE |>
  summarize(nb_signif = sum(p_adj < 0.05, na.rm = TRUE),
            nb_tested = n(),
            nb_peaky = sum(is_peaky),
            .by = cell_type) |>
  mutate(prop_peaky = nb_peaky / nb_tested) |>
  left_join(bulk_cat, by = "cell_type") |>
  # mutate(cell_type_noneur = if_else(cell_type %in% c("intestine", "ACh", "AMsh"), cell_type, "")) |>
  mutate(cell_type_noneur = if_else(tissue == "neuron", "", cell_type)) |>
  ggplot((aes(x = nb_tested, y = prop_peaky, color = tissue))) +
  theme_classic() +
  xlab("Number of genes tested") + ylab("Proportion peaky") +
  scale_y_continuous(labels = scales::label_percent(),
                     limits = c(0,1)) +
  scale_x_continuous(labels = scales::label_comma(),
                     limits = c(1,NA),
                     transform = "sqrt") +
  scale_shape_manual(values = c(19,8)) +
  scale_size_manual(values = c(1, 2.5)) +
  geom_point(aes(shape = bulk, size = bulk)) +
  ggrepel::geom_text_repel(aes(label = cell_type_noneur), show.legend = FALSE)

gg_peaky_by_celltype

# ggsave("peaky_genes_by_celltype.pdf", gg_peaky_by_celltype,
#        path = out_dir,
#        width = 20, height = 15, units = "cm")





# Compare known bulk ----

#~ overlap peaky/Osc ----

all_res_DE |>
  select(gene_id, gene_name, is_peaky) |>
  summarize(is_peaky = any(is_peaky),
            .by = c(gene_id, gene_name)) |>
  left_join(osc_raw, by = c(gene_id = "WB_ID")) |>
  rename(bulk = Class) |>
  mutate(sc = if_else(is_peaky, "peaky", "not peaky")) |>
  (\(df) table(df$bulk, df$sc))()


all_res_DE |>
  select(gene_id, gene_name, is_peaky) |>
  summarize(is_peaky = any(is_peaky),
            .by = c(gene_id, gene_name)) |>
  left_join(osc_raw, by = c(gene_id = "WB_ID")) |>
  mutate(sc = if_else(is_peaky, "peaky", "not peaky")) |>
  ggplot() +
  theme_classic() +
  xlab("Osc amplitude (bulk)") +
  geom_density(aes(x = OscAmplitude, fill = sc),
               alpha = .5)






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


all_res_DE |>
  filter(cell_type == "muscle") |>
  select(gene_id, gene_name, is_peaky) |>
  full_join(osc_raw |> select(WB_ID,OscAmplitude,"Class"),
            by = c(gene_id = "WB_ID")) |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  filter(gene_name %in% appg_genes_strict) |>
  mutate(bulk = case_when(
    OscAmplitude >= 1.5 ~ "Osc",
    OscAmplitude < 1.5 ~ "~",
    is.na(OscAmplitude) ~ "non Osc"
  )) |>
  mutate(sc = if_else(is_peaky, "peaky", "not peaky")) |>
  select(gene_name, bulk, sc) #|> clipr::write_clip()


# check genes ----




#~ focus on nonOsc clusters ----

all_res_DE |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "non oscillating") |>
  summarize(nb_tested = n(),
            nb_signif = sum(p_adj < 0.05, na.rm = TRUE),
            nb_peaky = sum(is_peaky),
            .by = c(cell_type, tissue)) |>
  mutate(prop_de = round(100 * nb_signif / nb_tested),
         prop_peaky = round(100 * nb_peaky / nb_tested)) |>
  relocate(tissue, .before = 1) |>
  arrange(desc(prop_peaky)) #|>  clipr::write_clip()





#~ focus on osc clusters ----

all_res_DE |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "oscillating") |>
  summarize(nb_tested = n(),
            nb_signif = sum(p_adj < 0.05, na.rm = TRUE),
            nb_peaky = sum(is_peaky),
            .by = c(cell_type, tissue)) |>
  mutate(prop_de = round(100 * nb_signif / nb_tested),
         prop_peaky = round(100 * nb_peaky / nb_tested)) |>
  relocate(tissue, .before = 1) |>
  arrange(tissue, cell_type) #|>  clipr::write_clip()




#~ intersections ----


#~ proportion of tested genes ----
degs_of_interest <- all_res_DE |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "oscillating")

intersections <- expand_grid(
  cell_1 = unique(degs_of_interest$cell_type),
  cell_2 = unique(degs_of_interest$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(degs_of_interest$gene_id[degs_of_interest$cell_type == .x &
                                                   degs_of_interest$is_peaky],
                        degs_of_interest$gene_id[degs_of_interest$cell_type == .y &
                                                   degs_of_interest$is_peaky]))),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(degs_of_interest$gene_id[degs_of_interest$cell_type == .x],
                    degs_of_interest$gene_id[degs_of_interest$cell_type == .y])))
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






#~ proportion of peaky genes ----
degs_of_interest <- all_res_DE |>
  left_join(bulk_cat,
            by = "cell_type") |>
  filter(bulk == "oscillating")

intersections <- expand_grid(
  cell_1 = unique(degs_of_interest$cell_type),
  cell_2 = unique(degs_of_interest$cell_type)
) |>
  mutate(
    intersection = map2_int(
      cell_1, cell_2,
      ~length(intersect(degs_of_interest$gene_id[degs_of_interest$cell_type == .x &
                                                   degs_of_interest$is_peaky],
                        degs_of_interest$gene_id[degs_of_interest$cell_type == .y &
                                                   degs_of_interest$is_peaky]))),
    union = map2_int(
      cell_1, cell_2,
      ~length(union(degs_of_interest$gene_id[degs_of_interest$cell_type == .x &
                                               degs_of_interest$is_peaky],
                    degs_of_interest$gene_id[degs_of_interest$cell_type == .y &
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






