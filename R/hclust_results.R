# Inits ----

library(tidyverse)


source("R/utils_fit.R")

dir_clust <- "intermediates/2502/250624_cluster"

dir_step2 <- "intermediates/2502/250624_step2"

# if working from external HDD
# dir_step2 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250624_step2"

dir_figures <- "presentations/figures/250624_clust_metrics/"



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
# png("presentations/figures/250624_hclust/heatmap_genes_pred.png",
#     width = 9, height = 4, units = "in", res = 500)
pheatmap::pheatmap(t(mat_pred),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # filename = "presentations/figures/250624_hclust/heatmap_genes_pred.pdf",
                   # width = 9, height = 4,
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
all_clustered_fits <- log1p(smooth_centered[,rownames(mat_pred)]) |>
  as.data.frame() |>
  rownames_to_column("time") |>
  as_tibble() |>
  pivot_longer(-time,
               names_to = "cellgene",
               values_to = "log_cnt") |>
  left_join(clusters,
            by = join_by(cellgene)) |>
  separate_wider_delim(cellgene,
                       delim = "|",
                       names = c("cell_type", "gene_name")) |>
  mutate(time = as.numeric(time))


set.seed(123)
selected <- all_clustered_fits |>
  select(cluster, cell_type, gene_name) |>
  distinct() |>
  group_by(cluster) |>
  slice_sample(n = 10) |>
  ungroup()

all_clustered_fits |>
  inner_join(selected) |>
  ggplot() +
  theme_classic() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  # facet_grid(rows = vars(cluster)) +
  facet_wrap(~ cluster) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed', color = 'grey80') +
  geom_ribbon(
    aes(x = time, ymin = average_signal - sd_signal, ymax = average_signal + sd_signal),
    alpha = .2,
    fill = "orange2",
    data = all_clustered_fits |>
      summarize(average_signal = mean(log_cnt),
                sd_signal = sd(log_cnt),
                .by = c(cluster, time))
  ) +
  geom_line(
    aes(x = time, y = log_cnt, group = interaction(cell_type, gene_name)),
    alpha = .4,
    linewidth = .2
  ) +
  geom_line(
    aes(x = time, y = average_signal),
    linewidth = 1.5,
    color = "orange2",
    data = all_clustered_fits |>
      summarize(average_signal = mean(log_cnt),
                sd_signal = sd(log_cnt),
                .by = c(cluster, time))
  )

# ggsave("cluster_average.png",
#        path = "presentations/figures/250624_hclust/",
#        width = 120, height = 75, units = "mm",
#        scale = 1.5)
# ggsave("cluster_average.pdf",
#        path = "presentations/figures/250624_hclust/",
#        width = 120, height = 75, units = "mm",
#        scale = 2)

# ggsave("cluster_average_vert.pdf",
#        path = "presentations/figures/250624_hclust/",
#        width = 70, height = 9*24, units = "mm",
#        scale = 2)





#~ save results ----

cluster_results <- all_clustered_fits |>
  select(cell_type, gene_name, cluster) |>
  distinct() |>
  mutate(shape = case_match(
    as.numeric(cluster),
    7 ~ "pulsatile",
    c(8,9) ~ "low",
    .default = "nonpulsatile"
  ))

# cluster_results |>
#   write_csv(file.path(dir_clust, "250624_cluster_results.csv"))





# # Compare euclidean and manhattan ----
# res_eucl <- read_csv(file.path(dir_clust, "250610_cluster_results.csv")) |>
#   mutate(cell_gene = paste0(cell_type, "_", gene_name))
# res_manh <- read_csv(file.path(dir_clust, "250610_manh_cluster_results.csv")) |>
#   mutate(cell_gene = paste0(cell_type, "_", gene_name))
# 
# list(eucl = res_eucl$cell_gene[res_eucl$shape == "pulsatile"],
#      manh = res_manh$cell_gene[res_manh$shape == "pulsatile"]) |>
#   eulerr::euler() |> plot(quantities = TRUE)
# 
# all.equal(res_eucl |> select(cell_type, gene_name), res_manh |> select(cell_type, gene_name))
# res <- cbind(
#   res_eucl |> select(cell_type, gene_name, cluster_eucl = cluster, shape_eucl = shape),
#   res_manh |> select(cluster_manh = cluster, shape_manh = shape)
# ) |>
#   as_tibble()
# 
# res |>
#   count(cluster_eucl, shape_eucl, cluster_manh, shape_manh) |>
#   mutate(cluster_eucl = as.factor(cluster_eucl),
#          cluster_manh = as.factor(cluster_manh)) |>
#   ggplot() +
#   theme_classic() +
#   scale_fill_viridis_c() +
#   geom_tile(aes(x = cluster_eucl, y = cluster_manh, fill = (n) ))
# 
# table(res$cluster_eucl, res$cluster_manh)
# 
# # examine genesets where clusterings disagree
# geneset <- res |>
#   filter(cluster_eucl == 4,
#          cluster_manh == 9) |>
#   mutate(cell_gene = paste0(cell_type, "|", gene_name)) |>
#   pull(cell_gene)
# 
# 
# length(geneset)
# 
# 
# all_clustered_fits <- log1p(smooth_centered[,geneset]) |>
#   as.data.frame() |>
#   rownames_to_column("time") |>
#   as_tibble() |>
#   pivot_longer(-time,
#                names_to = "cellgene",
#                values_to = "log_cnt") |>
#   separate_wider_delim(cellgene,
#                        delim = "|",
#                        names = c("cell_type", "gene_name")) |>
#   mutate(time = as.numeric(time))
# 
# 
# selected <- all_clustered_fits |>
#   select(cell_type, gene_name) |>
#   distinct() |>
#   slice_sample(n = 15)
# 
# all_clustered_fits |>
#   inner_join(selected) |>
#   ggplot() +
#   theme_classic() +
#   scale_color_brewer(type = "qual", palette = "Set2") +
#   scale_fill_brewer(type = "qual", palette = "Set2") +
#   # facet_grid(rows = vars(cluster)) +
#   geom_hline(aes(yintercept = 0),
#              linetype = 'dashed', color = 'grey80') +
#   geom_ribbon(
#     aes(x = time, ymin = average_signal - sd_signal, ymax = average_signal + sd_signal),
#     alpha = .2,
#     fill = "orange2",
#     data = all_clustered_fits |>
#       summarize(average_signal = mean(log_cnt),
#                 sd_signal = sd(log_cnt),
#                 .by = c(time))
#   ) +
#   geom_line(
#     aes(x = time, y = log_cnt, group = interaction(cell_type, gene_name)),
#     alpha = .4,
#     linewidth = .2
#   ) +
#   geom_line(
#     aes(x = time, y = average_signal),
#     linewidth = 1.5,
#     color = "orange2",
#     data = all_clustered_fits |>
#       summarize(average_signal = mean(log_cnt),
#                 sd_signal = sd(log_cnt),
#                 .by = c(time))
#   )
# 
# #> Looking at some of the disagreements:
# #> eucl  manh
# #>    5     9  some good, some bad; keep
# #>    7     3  more bad, discard
# #>    2     6  discard
# #>    2     7  discard
# #>    7     7  discard
# #>    4     9  keep
# 
# #>> keep the Euclidean





cluster_results |> filter(cell_type == "ILso") |> View()





# plot raw genes sc/PCA ----
dir_step1 <- "E:/backups/Projects_june2025/glia/osc_larva/intermediates/2502/250609_step1/"

ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu.qs")) )


goi <- "Y43F4A.1"


Seurat::FeaturePlot(ilso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "ILso")





geneset <- cluster_results |>
  filter(shape == "pulsatile",
         cell_type == "ILso") |>
  mutate(cell_gene = paste0(cell_type, "|", gene_name)) |>
  pull(cell_gene)


log1p(smooth_centered[,geneset]) |>
  matplot(type = "l")




# illustrate metrics ----


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


mod <- mods_centered[[goi]]

expr_smooth <- smooth_centered[,paste0("ILso|",goi)] / mean_sf




#~| cells ----

Seurat::FeaturePlot(ilso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "ILso")

# ggsave(paste0(goi, "_expr.png"),
#        path = dir_figures,
#        width = 50, height = 50, units = "mm",
#        scale = 2)
# ggsave(paste0(goi, "_expr.pdf"),
#        path = dir_figures,
#        width = 50, height = 50, units = "mm",
#        scale = 2)


#~| peak/baseline ----


dat <- tibble(
  pseudotime = seq_len(len) - 1,
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
  labs(title = goi,
       x = "pseudotime",
       y = "Smoothed expression") +
  geom_line(aes(x = pseudotime, y = expr), 
            linewidth = 1.2) +
  # scale_alpha_manual(values = c(`FALSE` = 0, `TRUE` = .5)) +
  # scale_fill_manual(values = c(`FALSE` = "green3", `TRUE` = "red4")) +
  geom_ribbon(aes(x = pseudotime, group = region_nb,
                  ymin = -Inf, ymax = Inf),
              dat = dat2,
              alpha = .2, fill = "brown3") +
  geom_hline(aes(
    yintercept = log10( 1 + all_descriptors$baseline[all_descriptors$cell_type == "ILso" &
                                                       all_descriptors$gene_name == goi]/mean_sf)
                 ),
    linetype = "dashed", color = "grey")


# ggsave(paste0(goi, "_baseline.png"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)
# ggsave(paste0(goi, "_baseline.pdf"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)




#~| dtw ----

tibble(
  x = seq_len(len) - 1,
  sig = expr_smooth,
  sig_norm = sig / max(sig),
  ref = ref
) |>
  ggplot() +
  theme_classic() +
  labs(title = goi,
       x = "pseudotime",
       y = "Smoothed expression") +
  geom_ribbon(aes(x = x, ymin = ref, ymax = sig_norm), 
              fill = alpha("cornsilk", 0.8), 
              color = NA) +
  geom_line(aes(x = x, y = sig_norm), 
            linewidth = 1.2) +
  geom_line(aes(x = x, y = ref), 
            linewidth = 1, 
            linetype = c("22"))

# ggsave(paste0(goi, "_dtw.png"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)
# ggsave(paste0(goi, "_dtw.pdf"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)





#~| deviance explained ----




# clipped

clip <- log10( 1 + mod$model$expr / exp(mod$model$`offset(log(size_factors))`) ) |>
  quantile(probs = .95)

ggplot() +
  theme_classic() +
  ylab("Expression (log, normalized)") +
  scale_y_continuous(limits = c(0, clip),
                     oob = scales::squish) +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod$model$pseudotime,
                           count = log10(1 + mod$model$expr / exp(mod$model$`offset(log(size_factors))`) )),
             alpha = .4,
             size = 2,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(
              pseudotime = (0:(len-1))/len,
              prediction = log10( 1 + expr_smooth )
            ),
            color = 'orange2',
            linewidth = 1.5) +
  ggtitle(goi)


# ggsave(paste0(goi, "_devexpl.png"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)
# ggsave(paste0(goi, "_devexpl.pdf"),
#        path = dir_figures,
#        width = 60, height = 50, units = "mm",
#        scale = 2)












