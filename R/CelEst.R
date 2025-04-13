
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
dir_step3 <- "intermediates/2502/250330_step3_genes_by_celltype//"


# download.file("https://raw.githubusercontent.com/IBMB-MFP/CelEsT-MS/refs/heads/main/CelEsT_annotated_v1pt1.txt",
#               "data/CelEst_v1pt1.txt")



# Load ----

cellest <- read_tsv("data/CelEst_v1pt1.txt")





cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- qs::qread(file.path(dir_step3, "all_genes.qs"))

stopifnot(all( unique(cell_types_info$cell_type) |> sort() %in%
                 unique(all_genes$cell_type) |> sort() ))


#~ check weights ----
table(cellest$weight)
all.equal(cellest$weight,
          rowSums(cellest[,c("with_motif" , "in_ChIP" , "in_eY1H")],
                  na.rm = TRUE) / 3)


#~ annotate gene names ----

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


cellest <- cellest |>
  mutate(target_id = wb_clean_gene_names(target, warn_missing = FALSE, return_one = TRUE),
         source_id = wb_clean_gene_names(source)) |>
  filter(! is.na(target_id)) |>
  mutate(target_name = i2s(target_id, gids, warn_missing = TRUE),
         source_name = i2s(source_id, gids, warn_missing = TRUE))


# Enrichment ----



cellest$source |> unique() |> length()
cellest$target |> unique() |> length()
expand_grid(unique(cellest$source), unique(cellest$target)) |> nrow()



all_tests <- map_dfr(all_genes$cell_type |> unique(),
                 \(ct){
                   message("---- ", ct)
                   
                   tested_genes <- all_genes$gene_name[all_genes$cell_type == ct]
                   
                   background_genes <- intersect(tested_genes,
                                                 cellest$target_name)
                   
                   message("background: ", length(background_genes)," genes")
                   
                   genes_osc <- all_genes |>
                     filter(cell_type == ct) |>
                     filter(peaky == "peak") |>
                     pull(gene_name) |>
                     intersect(background_genes)
                   
                   message("osc: ", length(genes_osc)," genes")
                   
                   
                   
                   
                   cellest_filt <- cellest |>
                     filter(target_name %in% background_genes,
                            source_name %in% tested_genes)
                   
                   
                   
                   candidates <- cellest_filt$source_name |> unique()
                   
                   p_vals <- map_dbl(candidates,
                                     \(candidate){
                                       
                                       nb_target_osc <- intersect(
                                         cellest_filt$target_name[cellest_filt$source_name == candidate],
                                         genes_osc
                                       ) |>
                                         length()
                                       
                                       nb_osc <- length(genes_osc)
                                       nb_nonosc <- background_genes |> setdiff(genes_osc) |> length()
                                       nb_target <- cellest_filt$target_name[cellest_filt$source_name == candidate] |> length()
                                       
                                       stats::phyper(q = nb_target_osc, m = nb_osc, n = nb_nonosc, k = nb_target,
                                                     lower.tail = FALSE)
                                     })
                   
                   tibble(cell_type = ct,
                          source_name = candidates,
                          p_val = p_vals)
                 }) |>
  mutate(p_adj = p.adjust(p_val, method = "BH"))

hist(all_tests$p_val)
hist(all_tests$p_adj)


table(all_tests$p_adj < .05)

all_tests |>
  split(all_tests$cell_type) |>
  map(~table(.x[["p_adj"]] < .05))

list(
  ILso = all_tests |>
    filter(cell_type == "ILso",
           p_adj < 0.05) |>
    pull(source_name),
  AM_PHso = all_tests |>
    filter(cell_type == "AM_PHso",
           p_adj < 0.05) |>
    pull(source_name),
  pharyngeal_muscle = all_tests |>
    filter(cell_type == "pharyngeal_muscle",
           p_adj < 0.05) |>
    pull(source_name),
  pharynx_epithelial = all_tests |>
    filter(cell_type == "pharynx_epithelial",
           p_adj < 0.05) |>
    pull(source_name)
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

xx <- all_tests |>
  filter(p_adj < 0.05) |>
  filter(! cell_type %in% c("mechanosensory",
                            "ACh_motoneuron","AIN",
                            "ASK","PVD","RMH","VD_DD"
                            ))

split(xx, xx$cell_type) |>
  map(~ .x |> pull(source_name)) |>
  UpSetR::fromList() |>
  UpSetR::upset(nsets = 50, nintersects = 200)

# are these TFs spiky themselves?
all_tests |>
  filter(p_adj < 0.05) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  summarize(nrows = n(),
            nb_tested = sum(!is.na(peaky)),
            nb_peaky = sum(peaky == "peak", na.rm = TRUE),
            .by = cell_type) |>
  arrange(nb_tested)


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
  filter(p_adj < 0.05) |>
  select(cell_type, source_name) |>
  left_join(all_genes,
            by = c("cell_type", source_name = "gene_name")) |>
  filter(cell_type == "ILso") %>%
  (\(.x) set_names(x = .x[["source_name"]], nm = .x[["peaky"]]))()
  



# network

library(igraph)

adj_list_ilso <- cellest_filt |>
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






