

# Inits ----
library(tidyverse)
library(wbData)
source("R/utils_heatmap_processing.R")

gids <- wb_load_gene_ids(295) |>
  add_row(X = "a",
          gene_id = "nsIs198",
          symbol = "GFP",
          sequence = "GFP",
          status = "Live",
          biotype = "protein_coding_gene",
          name = "GFP"
  )


dir_tf <- "intermediates/2502/250709_celest_tfs/"
dir_out <- "presentations/figures/250709_celest/"
dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"
dir_clust <- "intermediates/2502/250624_cluster"
dir_step2 <- "intermediates/2502/250624_step2/"
# dir_step2 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"


# download.file("https://raw.githubusercontent.com/IBMB-MFP/CelEsT-MS/refs/heads/main/CelEsT_annotated_v1pt1.txt",
#               "data/CelEst_v1pt1.txt")



# For gene name conversion, cf "step4" script
osc_table <- readxl::read_excel("data/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID, refresh = Inf),
         gene_name = i2s(gene_id, gids) ) |>
  filter(! is.na(gene_name)) |>
  mutate(osc_amplitude = if_else(is.na(OscAmplitude), 0, OscAmplitude)) |>
  select(gene_name, bulk_class = Class, osc_amplitude, bulk_peak = PeakPhase) |>
  group_by(gene_name) |>
  slice_max(osc_amplitude,
            with_ties = FALSE) |>
  ungroup()


stopifnot(!any(is.na(osc_table$gene_name)))
stopifnot(anyDuplicated(osc_table$gene_name) == 0L)




#~ load step 3 ----

cell_types_info <- qs::qread(file.path(dir_step3, "cell_types.qs"))

all_genes <- read_csv(file.path(dir_clust, "250624_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(cell_type %in% cell_types_info$cell_type)


stopifnot(all.equal(
  unique(cell_types_info$cell_type) |> sort(),
  unique(all_genes$cell_type) |> sort()
))



# select cell types to test

cell_types_osc <- cell_types_info |>
  filter(p_coherence_adj < .05,
         perplexity > 30) |>
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
  cell_types_osc,
  \(ct){
    
    all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
    puls_genes <- all_genes |> filter(cell_type == ct, shape == "pulsatile") |> pull(gene_name)
    nonpuls_genes <- all_genes_ct |> setdiff(puls_genes)
    
    stopifnot(all.equal(
      nonpuls_genes,
      all_genes |> filter(cell_type == ct, shape == "nonpulsatile" | shape == "low") |> pull(gene_name)
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
      
      observed <- contingency[1,1]
      expected <- contingency[1,2] * contingency[2,1]/sum(contingency)
      
      data.frame(
        cell_type = ct,
        source_name = tf,
        enrichment_fc = observed / expected,
        p_val = fisher.test(contingency, alternative = "greater")$p.value
      )
    })
    
    res |> bind_rows()
    
  },
  .progress = TRUE) |>
  as_tibble() |>
  mutate(p_adj = p.adjust(p_val, method = "BH"),
         signif = p_adj < .05)



hist(all_tests$p_val)
hist(all_tests$p_adj)

table(all_tests$p_adj < .05)

all_tests |>
  summarize(nb_tests = n(),
            nb_signif = sum(p_adj < .05),
            .by = cell_type)


# all_tests |>
#   # filter(cell_type == "hypodermis") |>
#   ggplot() +
#   theme_classic() +
#   scale_x_continuous(transform = "sqrt") +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .2)) +
#   scale_color_manual(values = c(`TRUE` = "red3", `FALSE` = "black")) +
#   facet_wrap(~cell_type) +
#   geom_hline(aes(yintercept = -log10(.05)),
#              linetype = "dashed", color = "grey") +
#   geom_point(aes(x = enrichment_fc, y = -log10(p_adj),
#                  alpha = signif, color = signif)) +
#   ggrepel::geom_text_repel(aes(x = enrichment_fc, y = -log10(p_adj),
#                                label = source_name),
#                            data = all_tests |>
#                              # filter(cell_type == "hypodermis") |>
#                              filter(signif)
#                            ,
#                            force_pull = .01,force = 10,
#                            max.overlaps = 10)
# 
# 
# all_tests |>
#   mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
#   ggplot() +
#   theme_classic() +
#   theme(legend.position = "none") +
#   scale_x_continuous(transform = "log2", labels = \(x) format(x, drop0trailing = TRUE)) +
#   scale_alpha_manual(values = c(`TRUE` = .8, `FALSE` = .2)) +
#   scale_color_manual(values = c(`TRUE` = "red3", `FALSE` = "black")) +
#   xlab("Fold Change (log)") +
#   facet_wrap(~cell_type) +
#   geom_hline(aes(yintercept = -log10(.05)),
#              linetype = "dashed", color = "grey") +
#   geom_point(aes(x = enrichment_fc, y = -log10(p_adj),
#                  alpha = signif, color = signif),
#              shape = 16, size = 1.5)





tfs_with_known_role <- tibble(
  source_name = c("nhr-23","grh-1","blmp-1","nhr-25","myrf-1","bed-3","nhr-85"),
  known_role = TRUE
)


tests_to_plot <- all_tests |>
  left_join(tfs_with_known_role,
            by = "source_name") |>
  left_join(all_genes |>
              select(cell_type, source_name = gene_name, shape),
            by = c("cell_type", "source_name")) |>
  left_join(cell_types_info |> select(cell_type, tissue),
            by = "cell_type") |>
  mutate(label = paste0('"', source_name, '"'),
         label = if_else(is.na(known_role),
                         label,
                         paste0("underline(", label,")")),
         label = if_else(shape == "pulsatile",
                         paste0("bolditalic(", label,")"),
                         paste0("italic(", label,")"))) |>
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
  mutate(tissue = factor(tissue,
                         levels = c("glia", "skin", "pharynx", "other"))) |>
  arrange(tissue, cell_type) |>
  mutate(cell_type = fct_inorder(cell_type) |> relevel(ref = "ILso"))

stopifnot(!any(is.na(tests_to_plot$tissue)))

# gg <- 
tests_to_plot |>
  # filter(cell_type == "glia sheath 2") |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, vjust = -5),
    strip.clip = "off"
  ) +
  scale_x_continuous(transform = "log2") +
  scale_alpha_manual(values = c(`TRUE` = .8, `FALSE` = .2)) +
  scale_fill_manual(values = c(`TRUE` = "orange", `FALSE` = "black")) +
  scale_color_manual(values = c(pulsatile = "red3", nonpulsatile = "black", low = "black")) +
  xlab("Fold Change (log)") +
  ylab(expression(-log[10](FDR))) +
  facet_wrap(~cell_type, axes = "all") +
  geom_hline(aes(yintercept = -log10(.05)),
             linetype = "dashed", color = "grey") +
  geom_point(aes(x = enrichment_fc, y = -log10(p_adj),
                 alpha = signif, fill = signif),
             shape = 21, stroke = NA, size = 2) +
  ggrepel::geom_text_repel(aes(x = enrichment_fc, y = -log10(p_adj),
                               label = label, color = shape),
                           data = tests_to_plot |>
                             # filter(cell_type == "glia sheath 2") |>
                             filter(signif)
                           ,
                           parse = TRUE,
                           size = 7/.pt, # convert mm to points
                           point.padding = unit(5, "mm"),
                           min.segment.length = unit(.5, "mm"),
                           force_pull = .005,force = 20,
                           direction = "x",
                           max.overlaps = 10)


# ggsave("volcano_TFs.pdf", path = dir_out,
#        width = 210, height = 150, units = "mm")


# For table
# tests_to_plot |>
#   mutate(
#     signif = if_else(signif, "*", ""),
#     known_role = if_else(!is.na(known_role) & known_role, "yes", ""),
#     TF_id = s2i(source_name, gids, warn_missing = TRUE)
#     ) |>
#   select(-label) |>
#   rename(TF_name = source_name) |>
#   relocate(tissue, TF_id, shape, .after = TF_name) |>
#   writexl::write_xlsx(file.path(dir_out, "table_s5_TFs_signif.xlsx"))




tests_to_plot |>
  filter(startsWith(cell_type, "pharyn")) |>
  filter(signif) |>
  summarize(FDR = min(p_adj),
            enr = min(enrichment_fc),
            n = n(),
            .by = source_name) |>
  filter(n > 1)

tests_to_plot |>
  filter(startsWith(cell_type, "pharyn")) |>
  filter(signif) |>
  pull(source_name) |>
  unique()

# tfs_signif <- all_tests |>
#   filter(any(signif),
#          .by = source_name) |>
#   pull(source_name) |>
#   unique()
# 
# stopifnot(identical(
#   tfs_signif,
#   all_tests |>
#     summarize(n_sig = sum(signif),
#               .by = source_name) |>
#     filter(n_sig > 0) |> pull(source_name)
# ))

mat_tfs <- all_tests |>
  filter(any(signif),
         .by = source_name) |>
  pivot_wider(id_cols = cell_type,
              names_from = source_name,
              values_from = p_val,
              values_fill = 1) |>
  column_to_rownames("cell_type") |>
  as.matrix()

pheatmap::pheatmap(mat_tfs)

hc_cts <- dist(mat_tfs) |> hclust()
hc_tfs <- dist(t(mat_tfs)) |> hclust()

all_tests |>
  filter(sum(p_adj < .05) > 0,
         .by = source_name) |>
  mutate(
    source_name = factor(source_name,
                         levels = rev(hc_tfs$labels[hc_tfs$order])),
    cell_type = factor(str_replace_all(cell_type, "_", " "),
                       levels = str_replace_all(hc_cts$labels, "_", " ")[hc_cts$order])
    ) |>
  ggplot() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    axis.text.y = element_text(face = "italic", size = 6),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(NULL) + ylab(NULL) +
  scale_color_gradient(low = "grey75", high = "red3", trans = c("log10", "reverse")) +
  scale_size_continuous(transform = "log2") +
  scale_alpha_continuous(transform = c("log10", "reverse")) +
  geom_point(aes(x = cell_type, y = source_name,
                 color = p_adj,
                 size = enrichment_fc,
                 alpha = p_adj),
             shape = 16)


# ggsave("dotplot_tfs_vs_celltypes.pdf",
#        path = dir_out,
#        width = 75, height = 175, units = "mm")



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
         p_coherence_adj < .05,
         perplexity > 30) |>
  select(cell_type,
         n_puls) |>
  left_join(all_tests,
            by = "cell_type") |>
  filter(p_adj < 0.05) |>
  # mutate(source_id = s2i(source_name, gids, warn_missing = TRUE)) |>
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





# Binned time enrichment ----




all(paste0(cell_types_osc, "_mods_uncentered.qs") %in% list.files(dir_step2))


bin_width <- 36
nb_bins <- 90

bins_start <- seq(0, 360 - bin_width, 360/nb_bins)
bins_end <- bins_start + bin_width

for(ct in cell_types_osc){
  
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
  
  png(file.path(dir_tf, paste0(ct, "_alignment.png")))
  plot(peak_times$bulk_peak,
       time_max_deg[peak_times$gene_name])
  dev.off()
  
  
  
  
  puls_genes <- all_genes |>
    filter(cell_type == ct,
           shape == "pulsatile") |>
    mutate(time_peak_deg = time_max_deg[gene_name])
  
  
  
  #~ overlapping bins ----
  
  all_tests_bin <- map_dfr(
    seq_along(bins_start),
    \(bin){
      
      all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
      
      puls_genes_bin <- puls_genes |>
        filter(time_peak_deg >= bins_start[[bin]],
               time_peak_deg < bins_end[[bin]] ) |>
        pull(gene_name)
      
      nonpuls_genes_bin <- all_genes_ct |> setdiff(puls_genes_bin)
      
      stopifnot(all.equal(
        nonpuls_genes_bin |> sort(),
        union(
          all_genes |> filter(cell_type == ct, shape == "nonpulsatile" | shape == "low") |> pull(gene_name),
          puls_genes |> filter(time_peak_deg >= bins_end[[bin]] | time_peak_deg < bins_start[[bin]]) |> pull(gene_name)
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
  
  qs::qsave(all_tests_bin,
            file.path(dir_tf, paste0(ct, "_all_tests_bin.qs")))
}




#~| Save heatmaps ----

for(ct in cell_types_osc){
  message(ct)
  
  all_tests_bin <- qs::qread(file.path(dir_tf, paste0(ct, "_all_tests_bin.qs")))
  
  signif_sources <- all_tests_bin |>
    filter(p_adj < 0.05) |>
    summarize(TFs = list(source_name),
              .by = time_bin) |>
    deframe() |> unlist() |> unique()
  
  if(length(signif_sources) < 2) next
  
  tf_by_time <- all_tests_bin |>
    filter(source_name %in% signif_sources) |>
    arrange(time_bin) |>
    mutate(signif = -log10(p_adj)) |>
    # mutate(signif = odds_ratio) |>
    pivot_wider(id_cols = source_name,
                names_from = time_bin,
                values_from = signif) |>
    column_to_rownames("source_name") |>
    as.matrix()
  
  
  pt_colnames <- round(bins_start + (bin_width/2), 1)
  pt_colnames[2 * (1:(length(pt_colnames)/2))] <- ""
  colnames(tf_by_time) <- pt_colnames
  
  
  # pheatmap::pheatmap(tf_by_time,
  #                    cluster_rows = TRUE,
  #                    clustering_distance_rows = "correlation",
  #                    cluster_cols = FALSE,
  #                    scale = "none")
  
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
  if(!is.null(dev.list()))  dev.off()
  
}



#### ____________________  ----

#~| res ----

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

all_tests_bin |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = time_bin) |>
  deframe() |> unlist() |> unique()




#~| heatmap of effect size ----
# signif_sources <- all_tests |> filter(cell_type == ct, p_adj < .05) |> pull(source_name)
signif_sources <- all_tests_bin |>
  filter(p_adj < 0.05) |>
  summarize(TFs = list(source_name),
            .by = time_bin) |>
  deframe() |> unlist() |> unique()

# tf_by_time <- all_tests_bin |>
#   filter(source_name %in% signif_sources) |>
#   arrange(time_bin) |>
#   pivot_wider(id_cols = source_name,
#               names_from = time_bin,
#               values_from = odds_ratio) |>
#   column_to_rownames("source_name") |>
#   as.matrix()
# 
# pheatmap::pheatmap(tf_by_time,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE,
#                    scale = "none")



tf_by_time <- all_tests_bin |>
  filter(source_name %in% signif_sources) |>
  arrange(time_bin) |>
  mutate(signif = -log10(p_adj)) |>
  pivot_wider(id_cols = source_name,
              names_from = time_bin,
              values_from = signif) |>
  column_to_rownames("source_name") |>
  as.matrix()


pt_colnames <- round(bins_start + (bin_width/2), 1)
pt_colnames[2 * (1:(length(pt_colnames)/2))] <- ""
colnames(tf_by_time) <- pt_colnames

pheatmap::pheatmap(tf_by_time,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "none")



# height 2.5 for ILso/main figure, height 3.5 for supp
pheatmap::pheatmap(
  tf_by_time,
  color = colorRampPalette(c("white", "#C994C7", "#DD1C77"))(100),
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  filename = paste0("presentations/figures/250616_celest/heatmap_timebins_", ct,".pdf"),
  width = 9, height = 3.5
)
dev.off()


# #~ non-overlapping bins ----
# 
# puls_genes_w_bin <- puls_genes |>
#   mutate(time_peak_binned = cut(time_peak_deg,
#                                  breaks = quantile(time_peak_deg,
#                                                    probs = seq(0, 1, 0.2)),
#                                 include.lowest = TRUE,
#                                 labels = FALSE))
# 
# 
# 
# 
# 
# all_tests_bin <- map_dfr(
#   unique(puls_genes_w_bin$time_peak_binned),
#   \(bin){
#     
#     all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
#     puls_genes_bin <- puls_genes_w_bin |> filter(time_peak_binned == bin) |> pull(gene_name)
#     nonpuls_genes_bin <- all_genes_ct |> setdiff(puls_genes_bin)
#     
#     stopifnot(all.equal(
#       nonpuls_genes_bin |> sort(),
#       union(
#         all_genes |> filter(cell_type == ct, shape == "nonpulsatile") |> pull(gene_name),
#         puls_genes_w_bin |> filter(time_peak_binned != bin) |> pull(gene_name)
#       ) |> sort()
#     ))
#     
#     
#     tfs_to_test_ct <- intersect(unique(celest$source_name),
#                                 all_genes_ct)
#     
#     
#     res <- lapply(tfs_to_test_ct,
#                   function(tf) {
#                     
#                     targets_tf <- celest |> filter(source_name == tf) |> pull(target_name)
#                     
#                     
#                     contingency <- matrix(c(
#                       length( targets_tf |> intersect(puls_genes_bin) ),
#                       length( puls_genes_bin |> setdiff(targets_tf) ),
#                       length( nonpuls_genes_bin |> intersect(targets_tf) ),
#                       length( nonpuls_genes_bin |> setdiff(targets_tf) )
#                     ), nrow = 2)
#                     
#                     
#                     data.frame(
#                       cell_type = ct,
#                       time_bin = bin,
#                       source_name = tf,
#                       odds_ratio = fisher.test(contingency, alternative = "greater")$estimate,
#                       p_val = fisher.test(contingency, alternative = "greater")$p.value
#                     )
#                   })
#     
#     res |> bind_rows()
#     
#   }) |>
#   as_tibble() |>
#   mutate(p_adj = p.adjust(p_val, method = "BH"))
# 
# 
# #~| res ----
# 
# all_tests_bin
# 
# 
# 
# hist(all_tests_bin$p_val)
# hist(all_tests_bin$p_adj)
# 
# table(all_tests_bin$p_adj < .05)
# 
# all_tests_bin |>
#   summarize(nb_tests = n(),
#             nb_signif = sum(p_adj < .05),
#             .by = time_bin)
# 
# 
# 
# # TFs enriched
# xx <- all_tests_bin |>
#   filter(p_adj < 0.05)
# 
# split(xx, xx$time_bin) |>
#   map(~ .x |> pull(source_name)) |>
#   UpSetR::fromList() |>
#   UpSetR::upset(nsets = 50, nintersects = 200)
# 
# 
# all_tests_bin |>
#   filter(p_adj < 0.05) |>
#   summarize(TFs = list(source_name),
#             .by = time_bin) |>
#   deframe()
# 
# 
# 
# # Display as table
# df <- all_tests_bin |>
#   filter(p_adj < 0.1) |>
#   select(time_bin, source_name) |>
#   pivot_wider(id_cols = time_bin,
#               names_from = source_name,
#               values_from = source_name,
#               values_fill = "")
# 
# df[order(rowSums(df == "")),order(colSums(df == ""))] |> as.data.frame() ##|> View()
# 
# 
# #~| heatmap of effect size ----
# signif_sources <- all_tests |> filter(cell_type == ct, p_adj < .05) |> pull(source_name)
# tf_by_time <- all_tests_bin |>
#   filter(source_name %in% signif_sources) |>
#   arrange(time_bin) |>
#   pivot_wider(id_cols = source_name,
#               names_from = time_bin,
#               values_from = odds_ratio) |>
#   column_to_rownames("source_name") |>
#   as.matrix()
# 
# pheatmap::pheatmap(tf_by_time,
#                    cluster_rows = TRUE,
#                    cluster_cols = FALSE,
#                    scale = "none")




# check the genes in regulated set ----


#~|~ Block A ----
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





#~|~ Block B ----
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

















