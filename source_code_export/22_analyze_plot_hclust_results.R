# Inits ----

library(tidyverse)
library(ggrastr)
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

source("R/utils_fit.R")

dir_clust <- "intermediates/2502/250624_cluster"

dir_step2 <- "intermediates/2502/250624_step2"

# if working from external HDD
# dir_step2 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250624_step2"

dir_figures <- "presentations/figures/260907_clust_metrics/"
# dir.create(dir_figures)


# Descriptors ----

# used for clustering
mat_pred <- qs::qread(file.path(dir_clust, "mat_predictors.qs"))


# more descriptors
all_descriptors <- list.files(dir_step2,
                         pattern = "_descriptors\\.qs$") |>
  map_dfr(~qs::qread(file.path(dir_step2, .x))) |>
  as_tibble()




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







# clusters ----

hc <- qs::qread(file.path(dir_clust, "hclust_euclidean.qs"))

ncl <- 9
memberships <- cutree(hc, k = ncl)
table(memberships)


# reorder cluster so they start at 1 in plot
old_clusters_order <- memberships[hc$order] |> unique()
lookup_cluster_names <- setNames(as.character(seq_len(ncl)), as.character(old_clusters_order))
memberships_reordered <- lookup_cluster_names[as.character(memberships)] |>
  setNames(names(memberships))


table(memberships_reordered)


clusters <- enframe(memberships_reordered,
                     name = "cellgene",
                     value = "cluster") |>
  mutate(cluster = as.factor(cluster))



annot_clusts <- clusters |>
  mutate(cluster = fct_rev(cluster)) |>
  arrange(desc(cluster)) |>
  column_to_rownames("cellgene")





#~ plot metrics ----

# > +++ Fig. 4C +++ ----

# png("presentations/figures/250624_hclust/heatmap_genes_pred_wide.png",
#     width = 11.25, height = 4, units = "in", res = 500)
pheatmap::pheatmap(t(mat_pred),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # filename = "presentations/figures/250624_hclust/heatmap_genes_pred_wide.pdf",
                   # width = 11.25, height = 4,
                   show_colnames = FALSE,
                   annotation_col = annot_clusts,
                   annotation_colors = list(cluster = set_names(pals::alphabet(ncl),
                                                                seq_len(ncl))),
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 7,
                                                                     name = "PiYG"))(100))

# dev.off()



#~ plot tree ----
plot(hc, labels = FALSE)
# abline(h = mean(rev(hc$height)[(ncl-1):ncl]))
rect.hclust(hc, k = ncl)



#~ average metrics by cluster ----
stopifnot(identical(rownames(mat_pred),
                    clusters$cellgene))

cluster_means <- aggregate(mat_pred, by=list(cluster=clusters$cluster), FUN=mean) |>
  as_tibble()

cluster_means |>
  pivot_longer(-cluster,
               names_to = "metric",
               values_to = "value") |>
  mutate(metric = fct_inorder(metric) |> fct_rev(),
         cluster = cluster) |>
  ggplot() +
  theme_classic() +
  xlab(NULL) + ylab(NULL) +
  scale_fill_distiller(palette = "PiYG", direction = 1) +
  geom_tile(aes(x = cluster, y = metric, fill = value))

# ggsave("metrics_average.png",
#        path = "presentations/figures/250624_hclust/",
#        width = 120, height = 75, units = "mm",
#        scale = 1.5)
# ggsave("metrics_average.pdf",
#        path = "presentations/figures/250624_hclust/",
#        width = 120, height = 75, units = "mm",
#        scale = 1.5)


# sub <- t(mat_pred[sample(rownames(mat_pred), 10),])
# hc_sub <- dendextend::prune(hc, setdiff(hc$labels, colnames(sub)))
# 
# pheatmap::pheatmap(sub,
#                    cluster_rows = FALSE,
#                    cluster_cols = hc_sub,
#                    # scale = "row",
#                    # show_rownames = FALSE,
#                    show_colnames = FALSE,
#                    annotation_col = annot_clusts)





#~ heatmap pseudotime ----

# > +++ Fig. 4D +++ ----


# png("presentations/figures/250624_hclust/manh_heatmap_genes_time_wide.png",
#     width = 11.25, height = 4, units = "in", res = 500)
pheatmap::pheatmap(log1p(smooth_centered[,rownames(mat_pred)]),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # filename = "presentations/figures/250624_hclust/manh_heatmap_genes_time.pdf",
                   # width = 9, height = 4,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_colors = list(cluster = set_names(pals::alphabet(ncl), seq_len(ncl))),
                   annotation_col = annot_clusts)

# dev.off()




#~ plot average curves ----

len <- nrow(smooth_centered)

all_clustered_fits <- log1p(smooth_centered[,rownames(mat_pred)]) |>
  as.data.frame() |>
  rownames_to_column("pseudotime") |>
  mutate(pseudotime = as.numeric(pseudotime)/len) |>
  as_tibble() |>
  pivot_longer(-pseudotime,
               names_to = "cellgene",
               values_to = "log_cnt") |>
  left_join(clusters,
            by = join_by(cellgene)) |>
  separate_wider_delim(cellgene,
                       delim = "|",
                       names = c("cell_type", "gene_name"))


set.seed(123)
selected_fits <- all_clustered_fits |>
  select(cluster, cell_type, gene_name) |>
  distinct() |>
  group_by(cluster) |>
  slice_sample(n = 10) |>
  ungroup() |>
  inner_join(all_clustered_fits)

fits_averaged_by_clust <- all_clustered_fits |>
  summarize(
    average_signal = mean(log_cnt),
    sd_signal = sd(log_cnt),
    .by = c(cluster, pseudotime)
  ) |>
  mutate(cluster = paste("cluster ", cluster))


# > +++ Fig. EV3D +++ ----

selected_fits |>
  mutate(cluster = paste("cluster ", cluster)) |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 10),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  ylab("Expression") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_x_continuous(breaks = 0:1) +
  scale_y_continuous(n.breaks = 3) +
  facet_wrap(~ cluster, scales = "free_y") +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed', color = 'grey80') +
  geom_ribbon(
    aes(x = pseudotime, ymin = average_signal - sd_signal, ymax = average_signal + sd_signal),
    alpha = .2,
    fill = "orange2",
    data = fits_averaged_by_clust
  ) +
  geom_line(
    aes(x = pseudotime, y = log_cnt, group = interaction(cell_type, gene_name)),
    alpha = .4,
    linewidth = .2
  ) +
  geom_line(
    aes(x = pseudotime, y = average_signal),
    linewidth = 1.5,
    color = "orange2",
    data = fits_averaged_by_clust
  )


# ggsave("cluster_average.pdf",
#        path = "presentations/figures/250624_hclust/",
#        width = 105, height = 96, units = "mm")












# Save results table S3 ----



cluster_results <- clusters |>
  separate_wider_delim(cellgene,
                       delim = "|",
                       names = c("cell_type", "gene_name")) |>
  mutate(shape = case_match(
    as.numeric(cluster),
    7 ~ "pulsatile",
    c(8,9) ~ "low",
    .default = "nonpulsatile"
  ))

# cluster_results |>
#   write_csv(file.path(dir_clust, "250624_cluster_results.csv"))




#~~ table S4, all genes and cell types as matrix ----

# > +++ Table EV4 +++ ----

cellgenes_mat <- cluster_results |>
  mutate(shape = case_match(shape,
                            "pulsatile" ~ "p",
                            "low" ~ "l",
                            "nonpulsatile" ~ "e")) |>
  pivot_wider(id_cols = gene_name,
              names_from = "cell_type",
              values_from = "shape",
              values_fill = "n") |>
  mutate(gene_id = s2i(gene_name, gids, warn_missing = TRUE),
         .before = 2) |>
  arrange(gene_name)
# writexl::write_xlsx(cellgenes_mat,
#                     file.path(dir_figures, "cellgenes_mat.xlsx"))






#~ Compute time of peak ----

source("R/utils_heatmap_processing.R") # --> z03_utils_heatmap_processing


# For gene name conversion, cf "24" script
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



# align each cell type's pseudotime on a common origin (dpy-6 following Sonntag et al)
origin_deg <- osc_table$bulk_peak[[ which( osc_table$gene_name == "dpy-6" ) ]]

grid_len <- 128
cell_types <- unique(cluster_results$cell_type)

res_peak_times <- vector("list", length(cell_types)) |>
  set_names(cell_types)
for(ct in cell_types){
  
  message(ct)
  
  mods_uncentered <- qs::qread(file.path(dir_step2, paste0(ct, "_mods_uncentered.qs")))
  
  # computed same for all genes
  mean_sf <- lapply(mods_uncentered,
                    \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
    unlist() |>
    log() |>
    mean() |>
    exp()
  
  
  
  preds_uncentered <- vapply(mods_uncentered,
                             \(.mod) predict(.mod,
                                             type = "response",
                                             newdata = data.frame(
                                               pseudotime = (0:(grid_len-1))/grid_len ,
                                               size_factors = rep(mean_sf, grid_len))
                             ),
                             FUN.VALUE = double(grid_len))
  
  
  time_max_pt <- apply(preds_uncentered, 2, which.max) / grid_len
  
  
  # align to bulk phase
  
  peak_times <- left_join(
    enframe(time_max_pt,
            name = "gene_name",
            value = "peak_pseudotime"),
    osc_table,
    by = "gene_name"
  ) |>
    filter(bulk_class == "Osc")
  
  alignment <- align_circular(peak_times$bulk_peak,
                              peak_times$peak_pseudotime*360)
  
  time_max_deg <- if (alignment$invert) {
    ((360 - time_max_pt*360) - alignment$shift) %% 360
  } else {
    (time_max_pt*360 - alignment$shift) %% 360
  }
  
  
  time_max_pct <- (100/360) * ( time_max_deg - origin_deg ) %% 360
  
  stopifnot(identical(
    names(mods_uncentered),
    names(time_max_pct)
  ))
  
  res_peak_times[[ ct ]] <- tibble(
    cell_type = ct,
    gene_name = names(mods_uncentered),
    peak_time_percent = time_max_pct
  )
}

peaks <- bind_rows(res_peak_times)

all.equal(cluster_results |> select(cell_type, gene_name),
          peaks |> select(cell_type, gene_name))

cluster_results_timed <- left_join(cluster_results,
                                   peaks,
                                   by = c("cell_type", "gene_name"))



#~ add predictors ----
all.equal(cluster_results_timed |> select(cell_type, gene_name),
          all_descriptors |> select(cell_type, gene_name))

preds <- mat_pred |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  separate_wider_delim(cellgene,
                       delim = "|",
                       names = c("cell_type", "gene_name"))
  
all.equal(cluster_results_timed |> select(cell_type, gene_name),
          preds |> select(cell_type, gene_name))

cluster_results_full <- bind_cols(
  cluster_results_timed,
  preds |> select(-cell_type, -gene_name) |> rename_with(~paste0("pred_",.x)),
  all_descriptors |> select(-cell_type, -gene_name) |> rename_with(~paste0("desc_",.x))
)


# > +++ Table EV3 +++ ----

# Export. Without rounding, 27 MB
# cluster_results_full |>
#   mutate(
#     across(peak_time_percent, ~ round(.x, 1)),
#     across(starts_with("pred"), ~ round(.x, 2)),
#     across(starts_with("desc"), ~ round(.x, 4))
#   ) |>
#   writexl::write_xlsx(file.path(dir_figures, "table_S3_cellgene_clusters.xlsx"))






# Agreement with bulk ----


dir_step3 <- "intermediates/2502/250624_step3_genes_by_celltype/"

cell_types_osc <- qs::qread(file.path(dir_step3, "cell_types.qs")) |>
  filter(p_coherence_adj < .05) |>
  pull(cell_type)

all_genes <- read_csv(file.path(dir_clust, "250624_cluster_results.csv")) |>
  filter(cell_type %in% cell_types_osc)



# see step_4 script for deduplication comments
osc_table <- readxl::read_excel("data/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wbData::wb_clean_gene_names(WB_ID),
         gene_name = wbData::i2s(gene_id, wbData::wb_load_gene_ids(295)) ) |>
  filter(! is.na(gene_name)) |>
  mutate(osc_amplitude = if_else(is.na(OscAmplitude), 0, OscAmplitude)) |>
  select(gene_name, bulk_class = Class, osc_amplitude) |>
  group_by(gene_name) |>
  slice_max(osc_amplitude,
            with_ties = FALSE) |>
  ungroup()



osc_genes_compare <- all_genes |>
  filter(cell_type %in% cell_types_osc) |>
  select(gene_name, shape) |>
  summarize(is_puls = any(shape == "pulsatile"),
            .by = gene_name) |>
  left_join(osc_table, by = "gene_name") |>
  mutate(`single-cell` = if_else(is_puls, "pulsatile", "nonpulsatile"))

osc_genes_compare |>
  (\(df) table(bulk = df$bulk_class, `single-cell` = df$`single-cell`, useNA = 'ifany'))()


# > +++ Fig. 4E +++ ----

# pdf(file.path(dir_figures, "osc_vs_pulsatile_euler.pdf"),
#     width = 2, height = 2)
list(
  bulk = osc_genes_compare$gene_name[which(osc_genes_compare$bulk_class == "Osc")],
  sc = osc_genes_compare$gene_name[which(osc_genes_compare$`single-cell` == "pulsatile")]
) |>
  eulerr::euler() |>
  plot(quantities = TRUE,
       edges = FALSE)

# dev.off()




## Compare OscAmplitude

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
#        path = dir_figures,
#        width = 75, height = 45, units = "mm")


## stat test

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  wilcox.test(osc_amplitude ~ `single-cell`, data = _)

osc_genes_compare |>
  filter(bulk_class == "Osc") |>
  summarize(median_amplitude = median(osc_amplitude),
            mad = mad(osc_amplitude),
            .by = `single-cell`)










# illustrate metrics ----

# > +++ Fig. EV3C +++ ----

dir_step1 <- "intermediates/2502/250609_step1"
# dir_step1 <- "E:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250609_step1"
# dir_step1 <- "D:/2025-06-27/Projects/glia/osc_larva/intermediates/2502/250609_step1/"

ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu_unsmoothed.qs")) )



mods_centered <- qs::qread(file.path(dir_step2, paste0("ILso", "_mods_centered.qs")))

len <- nrow(smooth_centered)


# computed same for all genes
mean_sf <- lapply(mods_centered,
                  \(.mod) exp(.mod$model$`offset(log(size_factors))`)) |>
  unlist() |>
  log() |>
  mean() |>
  exp()


ref <- dnorm( (seq_len(len) - 1)/len , mean = .5, sd = .01)
ref <- (ref - min(ref))/max(ref - min(ref))
ref <- circ_perm_mat(matrix(ref, ncol = 1))



#~ gene ----

goi <- "col-33"
goi <- "his-35"
goi <- "grl-18"


mod <- mods_centered[[goi]]

expr_smooth <- smooth_centered[,paste0("ILso|",goi)]




#~| cells ----


dat <- Seurat::FetchData(ilso_subseu, vars = c("PC_1","PC_2",goi))

# for ILso, invert axes for easier interpretation
dat$PC_1 <- -dat$PC_1



dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  labs(x = "PC 1", y = "PC 2") +
  scale_color_gradient(low = alpha("grey", .2), high = alpha("blue2", .8)) +
  # ggtitle(goi) +
  ggrastr::geom_point_rast(aes(x = PC_1, y = PC_2,
                               color = .data[[goi]]),
                           shape = 16,
                           size = 1,
                           raster.dpi = 500)



# ggsave(paste0(goi, "_expr.pdf"),
#        path = dir_figures,
#        width = 50, height = 45, units = "mm")




#~| peak/baseline ----


dat <- tibble(
  pseudotime = (seq_len(len) - 1)/(len - 1),
  expr = log10( 1 + expr_smooth )
) |>
  arrange(expr) |>
  mutate(is_min = rep(c(T,F), times = c(floor(.2*len), len - floor(.2*len)))) |>
  arrange(pseudotime)

dat2 <- dat |>
  mutate(change_region = is_min != lag(is_min, default = FALSE),
         region_nb = cumsum(change_region)) |>
  slice( c(1, n()),
         .by = region_nb) |>
  filter(is_min)

dat |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  labs(x = "pseudotime",
       y = "Expression") +
  scale_x_continuous(breaks = 0:1) +
  scale_y_continuous() +
  geom_ribbon(aes(x = pseudotime, group = region_nb,
                  ymin = -Inf, ymax = Inf),
              dat = dat2,
              alpha = .2, fill = "brown3") +
  geom_line(aes(x = pseudotime, y = expr), 
            linewidth = .8) +
  geom_hline(
    aes(
      yintercept = log10( 1 + all_descriptors$baseline[all_descriptors$cell_type == "ILso" &
                                                         all_descriptors$gene_name == goi])
    ),
    linetype = "dashed", color = "grey"
  ) +
  geom_segment(
    aes(x = pseudotime,
        y = 0,
        yend = expr),
    color = "cyan3",
    arrow = arrow(ends = "both",
                  length = unit(3, "mm"),
                  angle = 20),
    data = dat |> filter(expr == max(expr))
  )



# ggsave(paste0(goi, "_baseline.pdf"),
#        path = dir_figures,
#        width = 50, height = 45, units = "mm")




#~| dtw ----

tibble(
  pseudotime = (seq_len(len) - 1)/(len - 1),
  sig = expr_smooth,
  sig_norm = sig / max(sig),
  ref = ref
) |>
  ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  labs(x = "pseudotime",
       y = "Scaled expression") +
  scale_x_continuous(breaks = 0:1) +
  geom_ribbon(aes(x = pseudotime, ymin = ref, ymax = sig_norm), 
              fill = alpha("cornsilk", 0.8), 
              color = NA) +
  geom_line(aes(x = pseudotime, y = sig_norm), 
            linewidth = .8) +
  geom_line(aes(x = pseudotime, y = ref), 
            linewidth = .5, 
            linetype = c("22"))



# ggsave(paste0(goi, "_dtw.pdf"),
#        path = dir_figures,
#        width = 50, height = 45, units = "mm")





#~| deviance explained ----




# clipped

clip <- log10( 1 + mean_sf * mod$model$expr / exp(mod$model$`offset(log(size_factors))`) ) |>
  quantile(probs = .95)

ggplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  ylab("Expression") +
  scale_x_continuous(breaks = 0:1) +
  scale_y_continuous(limits = c(0, clip),
                     oob = scales::squish) +
  geom_point_rast(
    aes(x = pseudotime,
        y = count),
    data = tibble(pseudotime = mod$model$pseudotime,
                  count = log10(1 + mean_sf * mod$model$expr / exp(mod$model$`offset(log(size_factors))`) )),
    alpha = .3,
    size = 1,
    shape = 16,
    raster.dpi = 500
  ) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(
              pseudotime = (0:(len-1))/len,
              prediction = log10( 1 + expr_smooth  )
            ),
            color = 'orange2',
            linewidth = .8)




  
# ggsave(paste0(goi, "_devexpl.pdf"),
#        path = dir_figures,
#        width = 50, height = 45, units = "mm")




#~| metrics ----

all_descriptors |>
  filter(cell_type == "ILso",
         gene_name %in% c("grl-18", "his-35")) |>
  select(gene_name,
         baseline, max_peak, dist_dtw, dev_expl) |>
  mutate(baseline = log10(1 + baseline),
         max_peak = log10(1 + max_peak))


mat_pred[c("ILso|grl-18", "ILso|his-35"),]






