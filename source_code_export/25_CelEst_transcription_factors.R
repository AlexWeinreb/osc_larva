

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


dir_tf <- "intermediates/2502/260428_celest_tfs/"
dir_out <- "presentations/figures/260428_celest/"
dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"
dir_clust <- "intermediates/2502/250624_cluster"
dir_step2 <- "intermediates/2502/250624_step2/"
# dir_step2 <- "D:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250624_step2/"


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


# > +++ Fig. 6A +++ ----


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


# > +++ Table EV7 +++ ----

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
  mutate(cell_type = as.character(cell_type)) |>
  filter(startsWith(cell_type, "pharyn")) |>
  filter(signif) |>
  summarize(FDR = min(p_adj),
            enr = min(enrichment_fc),
            n = n(),
            .by = source_name) |>
  filter(n > 1)

tests_to_plot |>
  mutate(cell_type = as.character(cell_type)) |>
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


# > +++ Fig. EV6A +++ ----

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







# Binned time enrichment ----




all(paste0(cell_types_osc, "_mods_uncentered.qs") %in% list.files(dir_step2))


bin_width <- 10
nb_bins <- 100

bins_start <- seq(0, 100 - 100/nb_bins, 100/nb_bins)
bins_end <- bins_start + bin_width


origin_deg <- osc_table$bulk_peak[[ which( osc_table$gene_name == "dpy-6" ) ]]

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
  
  time_max_pct <- (100/360) * ( time_max_deg - origin_deg ) %% 360
  
  
  puls_genes <- all_genes |>
    filter(cell_type == ct,
           shape == "pulsatile") |>
    mutate(time_peak_pct = time_max_pct[gene_name])
  
  
  
  #~ overlapping bins ----
  
  all_tests_bin <- map_dfr(
    seq_along(bins_start),
    \(bin){
      
      all_genes_ct <- all_genes |> filter(cell_type == ct) |> pull(gene_name)
      
      puls_genes_bin <- puls_genes |>
        filter(time_peak_pct >= bins_start[[bin]],
               time_peak_pct < bins_end[[bin]] ) |>
        pull(gene_name)
      
      nonpuls_genes_bin <- all_genes_ct |> setdiff(puls_genes_bin)
      
      stopifnot(all.equal(
        nonpuls_genes_bin |> sort(),
        union(
          all_genes |> filter(cell_type == ct, shape == "nonpulsatile" | shape == "low") |> pull(gene_name),
          puls_genes |> filter(time_peak_pct >= bins_end[[bin]] | time_peak_pct < bins_start[[bin]]) |> pull(gene_name)
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
  
  
  pt_colnames <- bins_start
  pt_colnames[seq_along(pt_colnames) %% 5 != 1] <- ""
  colnames(tf_by_time) <- pt_colnames
  
  
  # pheatmap::pheatmap(tf_by_time,
  #                    cluster_rows = TRUE,
  #                    clustering_distance_rows = "correlation",
  #                    cluster_cols = FALSE,
  #                    scale = "none")
  
  n_tfs <- nrow(tf_by_time)
  
  # > +++ Fig. 6B +++ ----
  # > +++ Fig. EV6B +++ ----
  
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

# Check genes in regulated blocks ----

# This is included as a comment in Discussion, result not shown

#~ Method 1 ----

# only check genes pulsatile and target of block (do not look at the timing of peak)

#~| block A ----

tf_block_A <- c("blmp-1","nhr-85","nhr-23", "let-607", "klf-1", "nhr-41")
tf_block_A <- c("blmp-1","nhr-85","nhr-23")

list(
  puls = all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_A = celest |>
    filter(source_name %in% tf_block_A) |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

genes_puls_and_target_block_A <- intersect(
  all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    filter(source_name %in% tf_block_A) |>
    pull(target_name) |>
    unique()
)

# genes_puls_and_target_block_A |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_and_target_block_A.txt")



#~| block B ----
tf_block_B <- c("grh-1","lin-14","nhr-25")

list(
  puls = all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_B = celest |>
    filter(source_name %in% tf_block_B) |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


genes_puls_and_target_block_B <- intersect(
  all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    filter(source_name %in% tf_block_B) |>
    pull(target_name) |>
    unique()
)


# genes_puls_and_target_block_B |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_and_target_block_B.txt")


#~| any ----
# background list

list(
  puls = all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_any = celest |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


genes_puls_and_target_any <- intersect(
  all_genes |>
    filter(shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    pull(target_name) |>
    unique()
)


# genes_puls_and_target_any |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_and_target_any.txt")




#~~ compare ----

list(
  A = genes_puls_and_target_block_A,
  B = genes_puls_and_target_block_B,
  any = genes_puls_and_target_any
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)



# visually

setdiff(genes_puls_and_target_block_A,
        genes_puls_and_target_block_B) |>
  str_subset("-") |>
  sort()

setdiff(genes_puls_and_target_block_B,
        genes_puls_and_target_block_A) |>
  str_subset("-") |>
  sort()







#~ Method 1 ILso only ----
#~| block A ----

tf_block_A <- c("blmp-1","nhr-85","nhr-23", "let-607", "klf-1", "nhr-41")
tf_block_A <- c("blmp-1","nhr-85","nhr-23")

list(
  puls = all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_A = celest |>
    filter(source_name %in% tf_block_A) |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

genes_puls_in_ILso_and_target_block_A <- intersect(
  all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    filter(source_name %in% tf_block_A) |>
    pull(target_name) |>
    unique()
)

# genes_puls_in_ILso_and_target_block_A |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_in_ILso_and_target_block_A.txt")



#~| block B ----
tf_block_B <- c("grh-1","lin-14","nhr-25")

list(
  puls = all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_B = celest |>
    filter(source_name %in% tf_block_B) |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


genes_puls_in_ILso_and_target_block_B <- intersect(
  all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    filter(source_name %in% tf_block_B) |>
    pull(target_name) |>
    unique()
)


# genes_puls_in_ILso_and_target_block_B |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_in_ILso_and_target_block_B.txt")


#~| any ----
# background list

list(
  puls = all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  target_of_any = celest |>
    pull(target_name) |>
    unique()
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


genes_puls_in_ILso_and_target_any <- intersect(
  all_genes |>
    filter(cell_type == "ILso",
           shape == "pulsatile") |>
    pull(gene_name) |>
    unique(),
  celest |>
    pull(target_name) |>
    unique()
)


# genes_puls_in_ILso_and_target_any |> write_lines("intermediates/2502/250617_celest_tfs/genes_puls_in_ILso_and_target_any.txt")


#~~ compare ----

list(
  A = genes_puls_in_ILso_and_target_block_A,
  B = genes_puls_in_ILso_and_target_block_B,
  any = genes_puls_in_ILso_and_target_any
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


# visually

setdiff(genes_puls_in_ILso_and_target_block_A,
        genes_puls_in_ILso_and_target_block_B) |>
  str_subset("-") |>
  sort()

setdiff(genes_puls_in_ILso_and_target_block_B,
        genes_puls_in_ILso_and_target_block_A) |>
  str_subset("-") |>
  sort()





#~ Method 2 ----



# Reuse content of the previous loop to focus on a ct/timebin

ct <- "ILso"




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

plot(peak_times$bulk_peak,
     time_max_deg[peak_times$gene_name])

time_max_pct <- (100/360) * ( time_max_deg - origin_deg ) %% 360


puls_genes <- all_genes |>
  filter(cell_type == ct,
         shape == "pulsatile") |>
  mutate(time_peak_pct = time_max_pct[gene_name])





#~| Block A ----

time_bin_start <- 5
time_bin_end <- 30

tf_block_A <- c("blmp-1","nhr-85","nhr-23")




puls_genes_in_bin <- puls_genes |>
  filter(time_peak_pct >= time_bin_start,
         time_peak_pct <= time_bin_end ) |>
  pull(gene_name)

length(puls_genes_in_bin)


gene_puls_ILso_bin_and_block_A <- celest |>
  filter(source_name %in% tf_block_A,
         target_name %in% puls_genes_in_bin) |>
  pull(target_name) |> unique()


# celest |>
#   filter(source_name %in% c("blmp-1","nhr-85","nhr-23"),
#          target_name %in% puls_genes_in_bin) |>
#   select(4:6,9:10) |>
#   rowwise() |>
#   mutate(evidence = list(c(
#     "with_motif"[which(with_motif)],
#                             "in_ChIP"[which(in_ChIP)],
#                             "in_eY1H"[which(in_eY1H)]
#                            ))) |>
#   ungroup() |>
#   summarize(TF = paste(source_name, collapse = ", "),
#             evidence = paste(unique(unlist(evidence)), collapse = ", "),
#             .by = target_name) |>
#   flextable::flextable()





#~| Block B ----
time_bin_start <- 70
time_bin_end <- 90

tf_block_B <- c("grh-1","lin-14","nhr-25")




puls_genes_in_bin <- puls_genes |>
  filter(time_peak_pct >= time_bin_start,
         time_peak_pct <= time_bin_end ) |>
  pull(gene_name)

length(puls_genes_in_bin)


gene_puls_ILso_bin_and_block_B <- celest |>
  filter(source_name %in% tf_block_B,
         target_name %in% puls_genes_in_bin) |>
  pull(target_name) |> unique()



#~~ compare ----

list(
  A = gene_puls_ILso_bin_and_block_A,
  B = gene_puls_ILso_bin_and_block_B,
  any = genes_puls_in_ILso_and_target_any
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


# visually

setdiff(gene_puls_ILso_bin_and_block_A,
        gene_puls_ILso_bin_and_block_B) |>
  str_subset("-") |>
  sort()

setdiff(gene_puls_ILso_bin_and_block_B,
        gene_puls_ILso_bin_and_block_A) |>
  str_subset("-") |>
  sort()

# gene_puls_ILso_bin_and_block_A |> clipr::write_clip()
# gene_puls_ILso_bin_and_block_B |> clipr::write_clip()


#~~~ compare prev lists ----
list(
  bin = gene_puls_ILso_bin_and_block_A,
  all_time = genes_puls_in_ILso_and_target_block_A
) |>
  eulerr::euler() |>
  plot(quantities = TRUE)






# Regulators of particular genes ----
tfs_signif_in_ilso <- all_tests |> filter(cell_type == "ILso", p_adj < .05) |> pull(source_name)


celest |>
  filter(target_name %in% c("wrt-3", "dyf-7","cutl-16"),
         source_name %in% tfs_signif_in_ilso) |>
  rowwise() |>
  mutate(evidence = list(c(
    "with_motif"[which(with_motif)],
    "in_ChIP"[which(in_ChIP)],
    "in_eY1H"[which(in_eY1H)]
  ))) |>
  ungroup() |>
  summarize(evidence = paste(unique(unlist(evidence)), collapse = ", "),
            .by = c(target_name, source_name))







# ______________________ ----
# Relative timing TF vs targets ----


# For TFs that are themselves pulsatile, does their peak precede their targets?



#~ Collect peak times ----

all_peak_times <- map_dfr(cell_types_osc, \(ct) {
  
  message(ct)
  
  mods_uncentered <- qs::qread(file.path(dir_step2, paste0(ct, "_mods_uncentered.qs")))
  
  mean_sf <- lapply(mods_uncentered,
                    \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
    unlist() |> log() |> mean() |> exp()
  
  len <- 128
  
  preds_uncentered <- vapply(mods_uncentered,
                             \(.mod) predict(.mod,
                                             type = "response",
                                             newdata = data.frame(
                                               pseudotime = (0:(len - 1)) / len,
                                               size_factors = rep(mean_sf, len))
                             ),
                             FUN.VALUE = double(len))
  
  time_max <- apply(preds_uncentered, 2, which.max) / len
  
  # Align pseudotime to bulk phase
  peak_times <- left_join(
    enframe(time_max, name = "gene_name", value = "peak_pseudotime"),
    osc_table,
    by = "gene_name"
  ) |>
    filter(bulk_class == "Osc")
  
  alignment <- align_circular(peak_times$bulk_peak,
                              peak_times$peak_pseudotime * 360)
  
  time_max_deg <- if (alignment$invert) {
    ((360 - time_max * 360) - alignment$shift) %% 360
  } else {
    (time_max * 360 - alignment$shift) %% 360
  }
  
  
  # Convert to percent of stage, origin at dpy-6
  time_max_pct <- (100 / 360) * ((time_max_deg - origin_deg) %% 360)
  
  
  tibble(
    cell_type = ct,
    gene_name = names(time_max_pct),
    peak_pct = time_max_pct
  )
  
})

# qs::qsave(all_peak_times,
#           file.path(dir_tf, "260601_peak_times.qs"))

# all_peak_times <- qs::qread(file.path(dir_tf, "260601_peak_times.qs"))


signif_pulsatile_tfs <- all_tests |>
  filter(p_adj < 0.05) |>
  select(cell_type, source_name) |>
  inner_join(
    all_genes |> filter(shape == "pulsatile") |> select(cell_type, gene_name),
    by = c("cell_type", source_name = "gene_name")
  )


tf_target_timing <- signif_pulsatile_tfs |>
  pmap_dfr(\(cell_type, source_name) {
    
    # TF peak time
    tf_peak <- all_peak_times |>
      filter(cell_type == !!cell_type, gene_name == source_name) |>
      pull(peak_pct)
    
    if (length(tf_peak) == 0) return(tibble())
    
    # Pulsatile targets of this TF in this cell type
    targets <- celest |>
      filter(source_name == !!source_name) |>
      pull(target_name)
    
    pulsatile_targets <- all_genes |>
      filter(cell_type == !!cell_type, shape == "pulsatile",
             gene_name %in% targets) |>
      pull(gene_name)
    
    if (length(pulsatile_targets) == 0) return(tibble())
    
    # Target peak times
    target_peaks <- all_peak_times |>
      filter(cell_type == !!cell_type, gene_name %in% pulsatile_targets)
    
    target_peaks |>
      transmute(
        cell_type = !!cell_type,
        tf_name = source_name,
        target_name = gene_name,
        tf_peak_pct = tf_peak,
        target_peak_pct = peak_pct
      )
  })




tf_target_timing |>
  filter(cell_type == "ILso") |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = target_peak_pct),
                 binwidth = 5, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = tf_peak_pct),
             color = "red", linewidth = 1) +
  facet_wrap(~tf_name) +
  labs(x = "Peak time (% of larval stage)",
       y = "Number of pulsatile targets",
       subtitle = "Red line = TF peak") +
  xlim(0, 100)





tfs_to_show <- tf_target_timing |>
  summarize(n = n(), .by = c(tf_name, cell_type)) |>
  filter(n >= 100) |>
  count(tf_name) |>
  filter(n > 1) |>
  pull(tf_name) |>
  unique()

cts_to_show <- tf_target_timing |>
  summarize(n = n(), .by = c(tf_name, cell_type)) |>
  filter(n >= 100) |>
  count(cell_type) |>
  filter(n > 1) |>
  pull(cell_type) |>
  unique()

# plots_by_tf <- 


tf_target_timing |>
  rename(gene_name = tf_name) |>
  filter(gene_name %in% tfs_to_show,
         cell_type %in% cts_to_show) |>
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
  ggplot() +
  theme_minimal() +
  theme(panel.spacing.y = unit(0, "mm")) +
  theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 11, face = "italic")) +
  theme(panel.grid = element_blank()) +
  geom_density(aes(x = target_peak_pct),
               fill = "grey70") +
  geom_vline(data = tf_target_timing |>
               rename(gene_name = tf_name) |>
               filter(gene_name %in% tfs_to_show,
                      cell_type %in% cts_to_show) |>
               left_join(osc_table |>
                           filter(bulk_class == "Osc") |>
                           mutate(bulk_peak_pct = ((bulk_peak - origin_deg) %% 360) * 100 / 360) |>
                           select(gene_name, bulk_peak_pct),
                         by = "gene_name") |>
               select(gene_name, bulk_peak_pct) |>
               distinct(),
             aes(xintercept = bulk_peak_pct),
             color = "#457B9D", inewidth = 0.9, linetype = "25") +
  geom_vline(aes(xintercept = tf_peak_pct),
             color = "#E63946", linewidth = 1) +
  facet_grid(cell_type ~ gene_name, scale = "free_y")





#~ Test ----

# From Johnson 2023: https://journals.biologists.com/dev/article/150/10/dev201085/310520/NHR-23-activity-is-necessary-for-C-elegans
# "Most nhr-23-regulated genes involved in aECM structure/function, cholesterol metabolism,
# molting regulation, transcriptional regulation and signal transduction had peak amplitudes
# within 3 h of the nhr-23 expression peak (Fig. 3B)."


run_window_test <- function(before, after) {
  
  in_window <- function(target_pct, tf_pct) {
    diff <- (target_pct - tf_pct) %% 100
    diff <- if_else(diff > 50, diff - 100, diff)
    diff >= -before & diff <= after
  }
  
  signif_pulsatile_tfs |>
    pmap_dfr(\(cell_type, source_name) {
      
      tf_pct <- all_peak_times |>
        filter(cell_type == !!cell_type, gene_name == source_name) |>
        pull(peak_pct)
      
      if (length(tf_pct) == 0) return(tibble())
      
      # All pulsatile genes and their peaks in this cell type
      puls_ct <- all_genes |>
        filter(cell_type == !!cell_type, shape == "pulsatile") |>
        pull(gene_name)
      
      puls_peaks <- all_peak_times |>
        filter(cell_type == !!cell_type, gene_name %in% puls_ct)
      
      # TF targets among pulsatile genes
      targets <- celest |>
        filter(source_name == !!source_name) |>
        pull(target_name) |>
        intersect(puls_peaks$gene_name)
      
      non_targets <- setdiff(puls_peaks$gene_name, targets)
      
      if (length(targets) < 3) return(tibble())
      
      target_peaks <- puls_peaks |> filter(gene_name %in% targets)
      non_target_peaks <- puls_peaks |> filter(gene_name %in% non_targets)
      
      targets_in <- sum(in_window(target_peaks$peak_pct, tf_pct))
      targets_out <- nrow(target_peaks) - targets_in
      non_targets_in <- sum(in_window(non_target_peaks$peak_pct, tf_pct))
      non_targets_out <- nrow(non_target_peaks) - non_targets_in
      
      contingency <- matrix(c(targets_in, targets_out,
                              non_targets_in, non_targets_out), nrow = 2)
      
      ft <- fisher.test(contingency, alternative = "greater")
      
      tibble(
        cell_type = !!cell_type,
        tf_name = source_name,
        tf_peak_pct = tf_pct,
        n_targets = nrow(target_peaks),
        n_targets_in = targets_in,
        prop_targets_in = targets_in / nrow(target_peaks),
        prop_background_in = non_targets_in / nrow(non_target_peaks),
        odds_ratio = ft$estimate,
        p_val = ft$p.value
      )
    }) |>
    mutate(p_adj = p.adjust(p_val, method = "BH"))
}

windows <- list(
  narrow = c(before = 0,  after = 20),
  medium = c(before = 5, after = 30),
  wide   = c(before = 10, after = 40)
)

sensitivity <- map_dfr(names(windows), \(w) {
  message("Window: ", w)
  run_window_test(windows[[w]]["before"], windows[[w]]["after"]) |>
    mutate(window = w)
})

sensitivity |>
  summarize(
    n_tests = n(),
    n_signif = sum(p_adj < 0.05),
    median_odds = median(odds_ratio),
    median_prop_targets = median(prop_targets_in),
    median_prop_background = median(prop_background_in),
    .by = window
  )



tf_timings <- tf_target_timing |>
  rename(gene_name = tf_name) |>
  filter(gene_name %in% tfs_to_show,
         cell_type %in% cts_to_show) |>
  left_join(osc_table |>
              filter(bulk_class == "Osc") |>
              mutate(bulk_peak_pct = ((bulk_peak - origin_deg) %% 360) * 100 / 360) |>
              select(gene_name, bulk_peak_pct),
            by = "gene_name") |>
  mutate(cell_type = str_replace_all(cell_type, "_", " "))


tf_timings_segments <- tf_timings |>
  distinct(gene_name, cell_type, tf_peak_pct) |>
  crossing(
    do.call(rbind, windows) |>
      as.data.frame() |>
      rownames_to_column("window") |>
      mutate(
        # y_pos = c(0.0205, 0.019, 0.0175),
        y_pos = c(61, 0, 75)
        )
  ) |>
  left_join(sensitivity |>
              mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
              select(cell_type, gene_name = tf_name, window, p_adj),
            by = c("gene_name", "cell_type", "window")) |>
  mutate(x_start = (tf_peak_pct - before) %% 100,
         x_end   = (tf_peak_pct + after) %% 100,
         wraps = x_start > x_end)

tf_timings_segments_nowrap <- bind_rows(
  # Non-wrapping: keep as is
  tf_timings_segments |> filter(!wraps),
  # Wrapping: first part, from x_start to 100
  tf_timings_segments |> filter(wraps) |> mutate(x_end = 100),
  # Wrapping: second part, from 0 to x_end
  tf_timings_segments |> filter(wraps) |> mutate(x_start = 0)
) |>
  select(-wraps) |>
  filter(window == "medium")



# > +++ Fig. EV7 +++ ----

tf_timings |>
  ggplot() +
  theme_minimal() +
  theme(panel.spacing.y = unit(0, "mm")) +
  theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 11, face = "italic")) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 70)) +
  scale_color_manual(values = c(`TRUE` = "#2a9d8f", `FALSE` = "#f4a261"),
                     # guide = "none"
                     ) +
  labs(x = "Developmental progression (%)",
       y = "Number of target peaks") +
  geom_histogram(aes(x = target_peak_pct),
               fill = "grey70", color = "white", bins = 10) +
  geom_vline(data = tf_timings |>
               select(gene_name, bulk_peak_pct) |>
               distinct(),
             aes(xintercept = bulk_peak_pct),
             color = "#457B9D", linewidth = 0.3, linetype = "25") +
  geom_vline(aes(xintercept = tf_peak_pct),
             color = "#E63946", linewidth = 1) +
  facet_grid(cell_type ~ gene_name) +
  geom_segment(data = tf_timings_segments_nowrap |> rename(FDR = p_adj),
               aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos, color = FDR < 0.05),
               linewidth = 2)

# ggsave("260602_timing_tf_targets.pdf",
#        path = dir_out,
#        width = 8.5, height = 7, units = "in")
# 
# ggsave("260602_timing_tf_targets.png",
#        path = dir_out,
#        width = 8, height = 7, units = "in")




#~ No selection of ct/tf ----


tf_timings <- tf_target_timing |>
  rename(gene_name = tf_name) |>
  # filter(gene_name %in% tfs_to_show,
  #        cell_type %in% cts_to_show) |>
  left_join(osc_table |>
              filter(bulk_class == "Osc") |>
              mutate(bulk_peak_pct = ((bulk_peak - origin_deg) %% 360) * 100 / 360) |>
              select(gene_name, bulk_peak_pct),
            by = "gene_name") |>
  mutate(cell_type = str_replace_all(cell_type, "_", " "))


tf_timings_segments <- tf_timings |>
  distinct(gene_name, cell_type, tf_peak_pct) |>
  crossing(
    do.call(rbind, windows) |>
      as.data.frame() |>
      rownames_to_column("window") |>
      mutate(
        # y_pos = c(0.0205, 0.019, 0.0175),
        y_pos = c(61, 0, 75)
      )
  ) |>
  left_join(sensitivity |>
              mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
              select(cell_type, gene_name = tf_name, window, p_adj),
            by = c("gene_name", "cell_type", "window")) |>
  mutate(x_start = (tf_peak_pct - before) %% 100,
         x_end   = (tf_peak_pct + after) %% 100,
         wraps = x_start > x_end)

tf_timings_segments_nowrap <- bind_rows(
  # Non-wrapping: keep as is
  tf_timings_segments |> filter(!wraps),
  # Wrapping: first part, from x_start to 100
  tf_timings_segments |> filter(wraps) |> mutate(x_end = 100),
  # Wrapping: second part, from 0 to x_end
  tf_timings_segments |> filter(wraps) |> mutate(x_start = 0)
) |>
  select(-wraps) |>
  filter(window == "medium")


ct <- "pharynx epithelial"
tf_timings |>
  filter(cell_type == ct) |>
  ggplot() +
  theme_minimal() +
  theme(panel.spacing.y = unit(0, "mm")) +
  theme(strip.text.y = element_text(size = 10, angle = 0, hjust = 0),
        strip.text.x = element_text(size = 11, face = "italic")) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 70)) +
  scale_color_manual(values = c(`TRUE` = "#2a9d8f", `FALSE` = "#f4a261"),
                     # guide = "none"
  ) +
  labs(x = "Developmental progression (%)",
       y = "Number of target peaks") +
  geom_histogram(aes(x = target_peak_pct),
                 fill = "grey70", color = "white", bins = 10) +
  geom_vline(data = tf_timings |>
               filter(cell_type == ct) |>
               select(gene_name, bulk_peak_pct) |>
               distinct(),
             aes(xintercept = bulk_peak_pct),
             color = "#457B9D", linewidth = 0.3, linetype = "25") +
  geom_vline(aes(xintercept = tf_peak_pct),
             color = "#E63946", linewidth = 1) +
  facet_wrap( ~ gene_name) +
  geom_segment(data = tf_timings_segments_nowrap |>
                 filter(cell_type == ct) |> rename(FDR = p_adj),
               aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos, color = FDR < 0.05),
               linewidth = 2)


