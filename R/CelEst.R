
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


# note, Wormbase temporarily unavailable
wb_clean_gene_names <- function(...){
  wbData::wb_clean_gene_names(...,dir_cache = ".", refresh = Inf)
}

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
    
  }) |>
  as_tibble() |>
  mutate(p_adj = p.adjust(p_val, method = "BH"))



hist(all_tests$p_val)
hist(all_tests$p_adj)


all_tests |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = cell_type)

#####

table(all_tests$p_adj < .05)

all_tests |>
  split(all_tests$cell_type) |>
  map(~table(.x[["p_adj"]] < .05))

all_tests |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = cell_type)

# TFs enriched
xx <- all_tests |>
  filter(p_adj < 0.05)

split(xx, xx$cell_type) |>
  map(~ .x |> pull(source_name)) |>
  UpSetR::fromList() |>
  UpSetR::upset(nsets = 50, nintersects = 200)


all_tests |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = cell_type) |>
  pull(TFs)


df <- all_tests |>
  filter(p_adj < 0.1) |>
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
            nb_peaky = sum(shape == "pulsatile", na.rm = TRUE),
            .by = cell_type) |>
  arrange(nb_tested)

all_tests |>
  filter(p_adj < 0.1) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  select(source_name, shape) |>
  distinct() |>
  View()



####

cell_types_info |>
  filter(p_adj < 0.05) |>
  select(cell_type,
         nb_peaky_genes = nb_peaky) |>
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




# network

library(igraph)

## TODO
# adj_list_ilso <- 
celest |>
  filter(source_name %in% all_tests$source_name[all_tests$cell_type == "ILso"],
         target_name %in% all_genes$gene_name[all_genes$cell_type == "ILso"])
mutate(source_name = wb_clean_gene_names(source) |> i2s(gids, warn_missing = TRUE)) |>
  filter(source_name %in% xx$source_name[xx$cell_type == "ILso"]) |>
  select(source_name, target_name = gene_name) |>
  left_join(all_genes |> filter(cell_type == "ILso") |> select(gene_name, peaky),
            by = c(target_name = "gene_name"))

gr <- graph_from_edgelist(as.matrix(adj_list_ilso[,1:2]))

V(gr)$is_osc <- if_else(names(V(gr)) %in%
                          (all_genes |> filter(cell_type == "ILso") |> select(gene_name, peaky) |> pull(gene_name)),
                        "red",
                        "blue")

plot(gr, label = NA, label.cex = .1,labels=NA,
     vertex.color = V(gr)$is_osc)

# qs::qsave(gr, "test_igraph.qs")

# associated with time

all_genes |> filter(cell_type == "ILso")






