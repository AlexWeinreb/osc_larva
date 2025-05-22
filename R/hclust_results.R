# Inits ----

library(tidyverse)


source("R/utils_fit.R")

dir_step2 <- "intermediates/2502/250519_step2"


# mat_pred <- qs::qread(file.path("intermediates/2502/250516_step2", "250516_mat_predictors.qs"))
# hc <- qs::qread(file.path("intermediates/2502/250516_step2", "250516_hclust.qs"))



# Descriptors ----

# used for clustering
mat_pred <- qs::qread(file.path(dir_step2, "250519_mat_predictors.qs"))


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

hc <- qs::qread(file.path(dir_step2, "250519_hclust_euclidean.qs"))
hc <- qs::qread(file.path(dir_step2, "250519_hclust_manhattan.qs"))

ncl <- 10
memberships <- cutree(hc, k = ncl)




table(memberships)


annot_clusts <- enframe(memberships,
                        name = "cellgene",
                        value = "cluster") |>
  mutate(cluster = as.factor(cluster) |> fct_rev()) |>
  # mutate(is_cl=ifelse(cluster == 2, "y", "n")) |>
  arrange(desc(cluster)) |>
  column_to_rownames("cellgene")





# based on hc
pheatmap::pheatmap(t(mat_pred),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # scale = "row",
                   # show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_clusts,
                   annotation_colors = list(cluster = set_names(pals::alphabet(ncl), seq_len(ncl))),
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name =
                                                             "PiYG"))(100))


plot(hc, labels = FALSE)
abline(h = mean(rev(hc$height)[(ncl-1):ncl]))
rect.hclust(hc, k = ncl)



# average by cluster

cluster_means <- aggregate(mat_pred, by=list(cluster=memberships), FUN=mean) |>
  as_tibble()

cluster_means |>
  pivot_longer(-cluster,
               names_to = "metric",
               values_to = "value") |>
  mutate(metric = fct_inorder(metric) |> fct_rev(),
         cluster = factor(cluster, levels = unique( memberships[hc$labels[hc$order]] ))) |>
  ggplot() +
  theme_classic() +
  xlab(NULL) + ylab(NULL) +
  scale_fill_distiller(palette = "PiYG", direction = 1) +
  geom_tile(aes(x = cluster, y = metric, fill = value))



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


pheatmap::pheatmap(log1p(smooth_centered[,rownames(mat_pred)]),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # scale = "row",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_colors = list(cluster = set_names(pals::alphabet(ncl), seq_len(ncl))),
                   annotation_col = annot_clusts)




#~ plot average ----
all_clustered_fits <- log1p(smooth_centered[,rownames(mat_pred)]) |>
  as.data.frame() |>
  rownames_to_column("time") |>
  as_tibble() |>
  pivot_longer(-time,
               names_to = "cell_gene",
               values_to = "log_cnt") |>
  left_join(enframe(memberships,
                    name = "cell_gene",
                    value = "cluster"),
            by = join_by(cell_gene)) |>
  separate_wider_delim(cell_gene,
                       delim = "|",
                       names = c("cell_type", "gene_name")) |>
  mutate(cluster = as.factor(cluster),
         time = as.numeric(time))


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
  facet_wrap(~cluster) +
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




# save results ----

cluster_results <- all_clustered_fits |>
  select(cell_type, gene_name, cluster) |>
  distinct() |>
  mutate(shape = case_match(
    as.numeric(cluster),
    c(3,4,8,9)~ "pulsatile",
    .default = "nonpulsatile"
  ))

# cluster_results |>
#   write_csv(file.path(dir_step2, "250519_cluster_results.csv"))





cluster_results |> filter(cell_type == "ILso") |> View()





# plot raw ----
dir_step1 <- "intermediates/2502/250330_step1/"

ilso_subseu <- qs::qread( file.path(dir_step1,
                                    paste0("ILso", "_seu.qs")) )


goi <- "dpy-19"


Seurat::FeaturePlot(ilso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "ILso")




table(memberships == 10)

names(memberships)[memberships == 10] |> str_split_i(fixed("|"), 1) |> table() |> sort()


names(memberships)[memberships == 14] |> str_subset("^ILso")






log1p(smooth_centered[,names(memberships)[memberships == 14] |> str_subset("^ILso") ]) |>
  matplot(type = "l")
mat_sc[names(memberships)[memberships == 14] |> str_subset("^ILso"), ] |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |> View()


# illustrate metrics ----

colnames(mat_pred)

all_descriptors |>
  filter(cell_type == "ILso") |>
  arrange(desc(dist_dtw)) |> View()

goi <- "col-33"
plot(
  log10(1+smooth_centered[,paste0("ILso|",goi)]),
  type = "l",
  main = goi, xlab = "pseudotime", ylab = "log10( expression + 1 )",
  ylim = c(0,1.1*max(log10(1 + smooth_centered[,paste0("ILso|",goi)])))
)


#~ baseline ----
goi <- "act-4"
plot(
  log10(1+smooth_centered[,paste0("ILso|",goi)]),
  type = "l",
  main = goi, xlab = "pseudotime", ylab = "log10( expression + 1 )",
  ylim = c(0,1.1*max(log10(1 + smooth_centered[,paste0("ILso|",goi)])))
)

abline(h = log10( 1 + all_descriptors$baseline[all_descriptors$cell_type == "ILso" &
                                      all_descriptors$gene_name == goi]),
       lty = "dashed",
       col = "grey"
       )



#~ dtw ----


ref <- dnorm( (seq_len(len) - 1)/len , mean = .5, sd = .01)
ref <- (ref - min(ref))/max(ref - min(ref))
ref <- circ_perm_mat(matrix(ref, ncol = 1))

goi <- "col-33"
sig <- smooth_centered[,paste0("ILso|",goi)]
plot(
  sig/max(sig),
  type = "l",
  main = goi, xlab = "pseudotime", ylab = "log10( expression + 1 )",
  ylim = c(0,1)
)

polygon(x = c(seq_len(len) - 1, rev(seq_len(len))),
        y = c(ref, rev(sig/max(sig))),
        border = NA,
        col = alpha("cornsilk", .5))

lines((seq_len(len) - 1), ref, lwd = 2.5, lty = 'dotted')


#~ deviance explained ----

all_descriptors |>
  filter(cell_type == "ILso") |>
  arrange(dev_expl) |> View()

mods_centered <- qs::qread(file.path(dir_step2, paste0("ILso", "_mods_centered.qs")))



goi <- "his-35"
goi <- "col-33"

Seurat::FeaturePlot(ilso_subseu,
                    features = goi,
                    reduction = "pca",
                    pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                    alpha = .5) +
  ggtitle(goi, "ILso")


mod <- mods_centered[[goi]]


plot(
  log10(1+smooth_centered[,paste0("ILso|",goi)]),
  type = "l",
  main = goi, xlab = "pseudotime", ylab = "log10( expression + 1 )",
  ylim = c(0,1.1*max(log10(1 + smooth_centered[,paste0("ILso|",goi)])))
)



ggplot() +
  theme_classic() +
  ylab("Expression (log-count)") +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod$model$pseudotime,
                           count = log10(1 + mod$model$expr / exp(mod$model$`offset(log(size_factors))`) )),
             alpha = .2,
             size = 3,
             shape = 16) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod$model$pseudotime,
                          prediction = log10( 1 + mod$fitted.values / exp(mod$model$`offset(log(size_factors))`) )),
            color = 'red3') +
  ggtitle(goi)













