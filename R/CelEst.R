

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



out_dir <- "presentations/2503_figures/"
dir_step3 <- "intermediates/2502/250425_step3_genes_by_celltype/"
dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"
dir_step2_noncent <- "intermediates/2502/250422_step2_nb/"


# download.file("https://raw.githubusercontent.com/IBMB-MFP/CelEsT-MS/refs/heads/main/CelEsT_annotated_v1pt1.txt",
#               "data/CelEst_v1pt1.txt")




#~ load step 3 ----

cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name))


stopifnot(all.equal(
  unique(cell_types_info$cell_type) |> sort(),
  unique(all_genes$cell_type) |> sort()
))



# select cell types to test

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

cell_types_test <- cell_types_info |>
  filter(prop_puls >= 10,
         p_diag_adj < .05,
         p_coherence_adj < .05) |>
  pull(cell_type)








# Load CelEst ----

celest <- read_tsv("data/CelEst_v1pt1.txt")

#~ check weights ----
table(celest$weight)
all.equal(celest$weight,
          rowSums(celest[,c("with_motif" , "in_ChIP" , "in_eY1H")],
                  na.rm = TRUE) / 3)



#~ clean gene names ----

# # More comparison of gene correction
# comp <- tibble(
#   targets = unique(cellest$target),
#   updated = wb_clean_gene_names(targets, return_one = TRUE),
#   gene_id = wb_seq2id(targets, gids, warn_missing = TRUE)
# ) |>
#   mutate(compared = case_when(
#     is.na(updated) & is.na(gene_id) ~ "both NA",
#     is.na(updated) | is.na(gene_id) ~ "one NA",
#     updated == gene_id ~ "equal",
#     .default = "not equal"
#   ))
# 
# 
# comp$compared |> table()
# # not equal: wb_clean_name is right (dead genes)
# # one NA: wb_clean gets right when gene dead, seq2id gets right transposons
# # both NA: some dead unidentified genes


# # note, Wormbase temporarily unavailable
# wb_clean_gene_names <- function(...){
#   wbData::wb_clean_gene_names(...,dir_cache = ".", refresh = Inf)
# }

celest <- celest |>
  mutate(target_id = wb_clean_gene_names(target, warn_missing = FALSE, return_one = TRUE),
         source_id = wb_clean_gene_names(source)) |>
  filter(! is.na(target_id)) |>
  mutate(target_name = i2s(target_id, gids, warn_missing = TRUE),
         source_name = i2s(source_id, gids, warn_missing = TRUE))




# Enrichment ----

# Fisher exact test

all_tests <- map_dfr(
  cell_types_test,
  \(ct){
    
    all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
    puls_genes <- all_genes |> filter(cell_type == ct, shape == "pulsatile") |> pull(gene_name)
    nonpuls_genes <- all_genes_ct |> setdiff(puls_genes)
    
    stopifnot(all.equal(
      nonpuls_genes,
      all_genes |> filter(cell_type == ct, shape == "nonpulsatile") |> pull(gene_name)
    ))
    
    
    tfs_to_test_ct <- intersect(unique(celest$source_name),
                                all_genes_ct)
    
    
    res <- lapply(tfs_to_test_ct, function(tf) {
      
      targets_tf <- celest |> filter(source_name == tf) |> pull(target_name)
      
      
      contingency <- matrix(c(
        length( intersect(targets_tf, puls_genes) ),
        length( puls_genes |> setdiff(targets_tf) ),
        length( nonpuls_genes |> intersect(targets_tf) ),
        length( nonpuls_genes |> setdiff(targets_tf) )
      ), nrow = 2)
      
      
      data.frame(
        cell_type = ct,
        source_name = tf,
        p_val = fisher.test(contingency, alternative = "greater")$p.value
      )
    })
    
    res |> bind_rows()
    
  },
  .progress = TRUE) |>
  as_tibble() |>
  mutate(p_adj = p.adjust(p_val, method = "BH"))



hist(all_tests$p_val)
hist(all_tests$p_adj)

table(all_tests$p_adj < .05)

all_tests |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = cell_type)



#~ TFs enriched ----
# xx <- all_tests |>
#   filter(p_adj < 0.05)
# 
# split(xx, xx$cell_type) |>
#   map(~ .x |> pull(source_name)) |>
#   UpSetR::fromList() |>
#   UpSetR::upset(nsets = 50, nintersects = 200)


all_tests |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = cell_type) |>
  deframe()



# Display as table
df <- all_tests |>
  filter(p_adj < 0.05) |>
  select(cell_type, source_name) |>
  pivot_wider(id_cols = cell_type,
              names_from = source_name,
              values_from = source_name,
              values_fill = "")

df[order(rowSums(df == "")),order(colSums(df == ""))] |> as.data.frame() ##|> View()




# are these TFs pulsatile themselves?
all_tests |>
  filter(p_adj < 0.05) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  summarize(nrows = n(),
            nb_tested = sum(!is.na(shape)),
            nb_pulsatile = sum(shape == "pulsatile", na.rm = TRUE),
            .by = cell_type) |>
  arrange(nb_tested)

all_tests |>
  filter(p_adj < 0.1) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  select(source_name, cell_type, shape) |>
  distinct() |>
  View()



####

cell_types_info |>
  filter(prop_puls >= 10,
         p_diag_adj < .05,
         p_coherence_adj < .05) |>
  select(cell_type,
         n_puls) |>
  left_join(all_tests,
            by = "cell_type") |>
  filter(p_adj < 0.05) |>
  mutate(source_id = s2i(source_name, gids, warn_missing = TRUE)) |>
  filter(cell_type == "ILso") |> pull(source_name) |> paste(collapse = ", ")



all_tests |>
  filter(p_adj < 0.1) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  filter(cell_type == "ILso") %>%
  (\(.x) set_names(x = .x[["source_name"]], nm = .x[["shape"]]))()




# network ----

library(igraph)


gene_shape <- all_genes |>
  filter(cell_type == "ILso") |>
  select(gene_name, shape) |>
  distinct() |>
  column_to_rownames("gene_name")

sources_in_ilso <- all_tests |>
  filter(cell_type == "ILso") |>
  pull(source_name)

targets_of_source_in_ilso <- celest |>
  filter(source_name %in% sources_in_ilso) |>
  pull(target_name)

targets_in_ilso <- all_genes |>
  filter(
    cell_type == "ILso",
    shape == "pulsatile" | gene_name %in% targets_of_source_in_ilso
  ) |>
  pull(gene_name)

adj_list_ilso <- celest |>
  filter(source_name %in% sources_in_ilso,
         target_name %in% targets_in_ilso) |>
  left_join(all_genes |> filter(cell_type == "ILso") |> select(gene_name, shape),
            by = c(target_name = "gene_name")) |>
  select(source_name, target_name, weight, target_shape = shape)



gr <- graph_from_edgelist(as.matrix(adj_list_ilso[,1:2]))

V(gr)$is_puls <- gene_shape[names(V(gr)), "shape"]
V(gr)$type <- if_else(names(V(gr)) %in% adj_list_ilso$source_name,
                      "source",
                      "target")

graph_tbl <- as_tbl_graph(gr)

ggraph(graph_tbl, layout = "fr") +
  scale_size_manual(values = c(source = 2, target = .5)) +
  geom_node_point(aes(color = is_puls, shape = type, size = type)) +
  # geom_edge_link(alpha = .5) +
  theme_graph()

# subplot

adj_blmp1 <- adj_list_ilso |> filter(source_name == "blmp-1")
gr_blmp1 <- graph_from_edgelist(as.matrix(adj_blmp1[,1:2]))

V(gr_blmp1)$is_puls <- gene_shape[names(V(gr_blmp1)), "shape"]

# V(gr_blmp1)$is_puls <- case_match(
#   gene_shape[names(V(gr_blmp1)), "shape"],
#   "pulsatile" ~ "red3",
#   "nonpulsatile" ~ "blue3"
# )
V(gr_blmp1)$type <- if_else(names(V(gr_blmp1)) %in% adj_list_ilso$source_name,
                      "source",
                      "target")


graph_tbl <- as_tbl_graph(gr_blmp1)

ggraph(graph_tbl, layout = "fr") +
  scale_size_manual(values = c(source = 2, target = .8)) +
  geom_edge_link(alpha = .1) +
  geom_node_point(aes(color = is_puls, shape = type, size = type)) +
  theme_graph()


# qs::qsave(gr, "test_igraph.qs")



# time non-overlapping bins ----

list.files(dir_step2_noncent)

preds_ilso <- file.path(dir_step2_noncent, "ILso_preds.qs") |> qs::qread()

n_timepoints <- nrow(preds_ilso)

time_max <- apply(preds_ilso, 2, which.max) / n_timepoints

puls_genes <- all_genes |>
  filter(cell_type == "ILso",
         shape == "pulsatile") |>
  mutate(time_peak = time_max[gene_name],
         time_peak_binned = cut(time_peak,
                                 breaks = quantile(time_peak,
                                                   probs = seq(0, 1, 0.2)),
                                include.lowest = TRUE,
                                labels = FALSE))





ct <- "ILso"
all_tests_bin <- map_dfr(
  unique(puls_genes$time_peak_binned),
  \(bin){
    
    all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
    puls_genes_bin <- puls_genes |> filter(time_peak_binned == bin) |> pull(gene_name)
    nonpuls_genes_bin <- all_genes_ct |> setdiff(puls_genes_bin)
    
    stopifnot(all.equal(
      nonpuls_genes_bin |> sort(),
      union(
        all_genes |> filter(cell_type == ct, shape == "nonpulsatile") |> pull(gene_name),
        puls_genes |> filter(time_peak_binned != bin) |> pull(gene_name)
      ) |> sort()
    ))
    
    
    tfs_to_test_ct <- intersect(unique(celest$source_name),
                                all_genes_ct)
    
    
    res <- lapply(tfs_to_test_ct,
                  function(tf) {
                    
                    targets_tf <- celest |> filter(source_name == tf) |> pull(target_name)
                    
                    
                    contingency <- matrix(c(
                      length( targets_tf |> intersect(puls_genes_bin) ),
                      length( puls_genes_bin |> setdiff(targets_tf) ),
                      length( nonpuls_genes_bin |> intersect(targets_tf) ),
                      length( nonpuls_genes_bin |> setdiff(targets_tf) )
                    ), nrow = 2)
                    
                    
                    data.frame(
                      cell_type = ct,
                      time_bin = bin,
                      source_name = tf,
                      odds_ratio = fisher.test(contingency, alternative = "greater")$estimate,
                      p_val = fisher.test(contingency, alternative = "greater")$p.value
                    )
                  })
    
    res |> bind_rows()
    
  }) |>
  as_tibble() |>
  mutate(p_adj = p.adjust(p_val, method = "BH"))


#~ res ----

all_tests_bin



hist(all_tests_bin$p_val)
hist(all_tests_bin$p_adj)

table(all_tests_bin$p_adj < .05)

all_tests_bin |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = time_bin)



# TFs enriched
xx <- all_tests_bin |>
  filter(p_adj < 0.05)

split(xx, xx$time_bin) |>
  map(~ .x |> pull(source_name)) |>
  UpSetR::fromList() |>
  UpSetR::upset(nsets = 50, nintersects = 200)


all_tests_bin |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = time_bin) |>
  deframe()



# Display as table
df <- all_tests_bin |>
  filter(p_adj < 0.1) |>
  select(time_bin, source_name) |>
  pivot_wider(id_cols = time_bin,
              names_from = source_name,
              values_from = source_name,
              values_fill = "")

df[order(rowSums(df == "")),order(colSums(df == ""))] |> as.data.frame() ##|> View()

# heatmap of effect size ----
signif_sources <- all_tests |> filter(cell_type == "ILso", p_adj < .05) |> pull(source_name)
tf_by_time <- all_tests_bin |>
  filter(source_name %in% signif_sources) |>
  arrange(time_bin) |>
  pivot_wider(id_cols = source_name,
              names_from = time_bin,
              values_from = odds_ratio) |>
  column_to_rownames("source_name") |>
  as.matrix()

pheatmap::pheatmap(tf_by_time,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "none")




# time overlapping bins ----

ct <- "ILso"

preds_ct <- file.path(dir_step2_noncent, paste0(ct, "_preds.qs")) |> qs::qread()

n_timepoints <- nrow(preds_ct)

time_max <- apply(preds_ct, 2, which.max) / n_timepoints

puls_genes <- all_genes |>
  filter(cell_type == ct,
         shape == "pulsatile") |>
  mutate(time_peak = time_max[gene_name])

bin_width <- .1
bins_start <- seq(0, 1 - bin_width, .01)
bins_end <- bins_start + bin_width


all_tests_bin <- map_dfr(
  seq_along(bins_start),
  \(bin){
    
    all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
    
    puls_genes_bin <- puls_genes |>
      filter(time_peak >= bins_start[[bin]],
             time_peak < bins_end[[bin]] ) |>
      pull(gene_name)
    
    nonpuls_genes_bin <- all_genes_ct |> setdiff(puls_genes_bin)
    
    stopifnot(all.equal(
      nonpuls_genes_bin |> sort(),
      union(
        all_genes |> filter(cell_type == ct, shape == "nonpulsatile") |> pull(gene_name),
        puls_genes |> filter(time_peak >= bins_end[[bin]] | time_peak < bins_start[[bin]]) |> pull(gene_name)
      ) |> sort()
    ))
    
    
    tfs_to_test_ct <- intersect(unique(celest$source_name),
                                all_genes_ct)
    
    
    res <- lapply(tfs_to_test_ct,
                  function(tf) {
                    
                    targets_tf <- celest |> filter(source_name == tf) |> pull(target_name)
                    
                    
                    contingency <- matrix(c(
                      length( targets_tf |> intersect(puls_genes_bin) ),
                      length( puls_genes_bin |> setdiff(targets_tf) ),
                      length( nonpuls_genes_bin |> intersect(targets_tf) ),
                      length( nonpuls_genes_bin |> setdiff(targets_tf) )
                    ), nrow = 2)
                    
                    
                    data.frame(
                      cell_type = ct,
                      time_bin = bin,
                      source_name = tf,
                      odds_ratio = fisher.test(contingency, alternative = "greater")$estimate,
                      p_val = fisher.test(contingency, alternative = "greater")$p.value
                    )
                  })
    
    res |> bind_rows()
    
  },
  .progress = TRUE) |>
  as_tibble() |>
  mutate(p_adj = p.adjust(p_val, method = "BH"))


#~ res ----

all_tests_bin



hist(all_tests_bin$p_val)
hist(all_tests_bin$p_adj)

table(all_tests_bin$p_adj < .05)

all_tests_bin |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = time_bin)



# TFs enriched


all_tests_bin |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = time_bin) |>
  deframe()





#~ heatmap of effect size ----
signif_sources <- all_tests |> filter(cell_type == ct, p_adj < .05) |> pull(source_name)
tf_by_time <- all_tests_bin |>
  filter(source_name %in% signif_sources) |>
  arrange(time_bin) |>
  pivot_wider(id_cols = source_name,
              names_from = time_bin,
              values_from = odds_ratio) |>
  column_to_rownames("source_name") |>
  as.matrix()

pheatmap::pheatmap(tf_by_time,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "none")








#~ check the genes in regulated set ----


#~~ Block A ----
tf_by_time[c("blmp-1","nhr-85","nhr-23"), 1:10] |>
  pheatmap::pheatmap(cluster_rows = FALSE,
                     cluster_cols = FALSE)


puls_genes_in_bin <- puls_genes |>
  filter(time_peak >= bins_start[[1]],
         time_peak < bins_end[[5]] ) |>
  pull(gene_name)

length(puls_genes_in_bin)

celest |>
  filter(source_name %in% c("blmp-1","nhr-85","nhr-23"),
         target_name %in% puls_genes_in_bin) |>
  pull(target_name) |> unique() |> length() #paste(collapse = ", ")

celest |>
  filter(source_name %in% c("blmp-1","nhr-85","nhr-23"),
         target_name %in% puls_genes_in_bin) |>
  select(4:6,9:10) |>
  rowwise() |>
  mutate(evidence = list(c(
    "with_motif"[which(with_motif)],
                            "in_ChIP"[which(in_ChIP)],
                            "in_eY1H"[which(in_eY1H)]
                           ))) |>
  ungroup() |>
  summarize(TF = paste(source_name, collapse = ", "),
            evidence = paste(unique(unlist(evidence)), collapse = ", "),
            .by = target_name) |>
  flextable::flextable()





#~~ Block B ----
tf_by_time[c("eor-1","nhr-41","nhr-25"), 24:50] |> 
  pheatmap::pheatmap(cluster_rows = FALSE,
                     cluster_cols = FALSE)


puls_genes_in_bin <- puls_genes |>
  filter(time_peak >= bins_start[[28]],
         time_peak < bins_end[[44]] ) |>
  pull(gene_name)

length(puls_genes_in_bin)

celest |>
  filter(source_name %in% c("eor-1","nhr-41","nhr-25"),
         target_name %in% puls_genes_in_bin) |>
  pull(target_name) |> unique() |> paste(collapse = ", ")#length() #paste(collapse = ", ")

celest |>
  filter(source_name %in% c("eor-1","nhr-41","nhr-25"),
         target_name %in% puls_genes_in_bin) |>
  rowwise() |>
  mutate(evidence = list(c(
    "with_motif"[which(with_motif)],
    "in_ChIP"[which(in_ChIP)],
    "in_eY1H"[which(in_eY1H)]
  ))) |>
  ungroup() |>
  summarize(TF = paste(source_name, collapse = ", "),
            evidence = paste(unique(unlist(evidence)), collapse = ", "),
            .by = target_name) |>
  flextable::flextable()







# Regulators of particular genes ----
signif_sources <- all_tests |> filter(cell_type == "ILso", p_adj < .05) |> pull(source_name)

celest |>
  filter(target_name == "grl-18",
         source_name %in% signif_sources)



celest |>
  filter(target_name %in% c("wrt-3", "dyf-7","cutl-16"),
         source_name %in% signif_sources) |>
  rowwise() |>
  mutate(evidence = list(c(
    "with_motif"[which(with_motif)],
    "in_ChIP"[which(in_ChIP)],
    "in_eY1H"[which(in_eY1H)]
  ))) |>
  ungroup() |>
  summarize(evidence = paste(unique(unlist(evidence)), collapse = ", "),
            .by = c(target_name, source_name))








