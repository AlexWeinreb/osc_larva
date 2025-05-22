# test clustering of genes to find pulsatile ones


# Inits ----

library(tidyverse)

source("R/utils_fit.R")


opar <- par(no.readonly = TRUE)


#~ load ----
dir_step2 <- "intermediates/2502/250514_step2/"

dir_out_scvelo <- "intermediates/2502/250501_scvelo"




pred_manual <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

annot_pred_manual_amphso <- pred_manual |>
  dplyr::select(manual, gene_name) |>
  mutate(gene_name = paste0("AM_PHso|", gene_name)) |>
  column_to_rownames("gene_name")

manual_amphso <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

manual_several <- readxl::read_excel("intermediates/2502/250330_step2/manual_annotation.xlsx") |>
  dplyr::select(-amplitude, -dev_expl)




smooth_preds <- list.files(dir_step2,
                           pattern = "_preds\\.qs$") |>
  enframe(value = "filename",
          name = NULL) |>
  separate_wider_regex(filename,
                       patterns = c(
                         cell_type = "^.+",
                         "_preds\\.qs"
                       ),
                       cols_remove = FALSE) |>
  pmap(\(cell_type, filename){
    mat_preds <- qs::qread(file.path(dir_step2,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)


smooth_centered <- circ_perm_mat(smooth_preds)



# matplot(log1p(smooth_preds)[,1:200], type = "l")
# matplot(log1p(smooth_centered)[,1:200], type = "l")
# printMat::matimage(log1p(smooth_preds)[,1:100])
# printMat::matimage(log1p(smooth_centered)[,1:100])

stopifnot(all(paste0("AM_PHso|", pred_manual$gene_name) %in% colnames(smooth_centered)))

# pheatmap::pheatmap(log1p(smooth_centered)[,paste0("AM_PHso|", pred_manual$gene_name)],
#                    cluster_rows = FALSE,
#                    cluster_cols = TRUE,
#                    # scale = "column",
#                    # clustering_distance_cols = as.dist(dist_mat),
#                    annotation_col = pred_manual |>
#                      mutate(cellgene = paste0(cell_type,"|",gene_name)) |>
#                      dplyr::select(manual, cellgene) |>
#                      column_to_rownames("cellgene"))









# Prepare features ----

# note this part is slow, we can skip to cached

predictors <- list.files(dir_step2,
                      pattern = "_descriptors\\.qs$") |>
  map_dfr(~qs::qread(file.path(dir_step2, .x))) |>
  as_tibble()




# Summarize s1, s2, s3 (quite correlated)
predictors <- predictors |>
  mutate(cellgene = paste0(cell_type, "|", gene_name),
         .before = 3) |>
  rowwise() |>
  mutate(var_s = var(c(s1, s2, s3, s4)),
         max_s = max(c(s1, s2, s3, s4))) |>
  ungroup()

predictors


stopifnot(all.equal(predictors$cellgene,
                    colnames(smooth_centered)))

len <- nrow(smooth_centered)
predictors$baseline <- apply(smooth_centered, 2,
                   \(expr){
                     
                     starts <- seq_len(.1*len)
                     ends <- rev(len - starts + 1)
                     
                     x <- c(starts, ends)
                     y <- expr[x]
                     
                     mean(y)
                   })


#~ add scVelo ----

res_scvelo <- list.files(dir_out_scvelo, pattern = "_scvelo_fit\\.qs$") |>
  str_remove("_scvelo_fit\\.qs$") |>
  map_dfr(\(ct){
    qs::qread(file.path(dir_out_scvelo, paste0(ct, "_scvelo_fit.qs")))
  }) |>
  select(-gene_id)



predictors <- predictors |>
  left_join(res_scvelo,
            by = join_by(cell_type, gene_name)) |>
  mutate(fit_likelihood = replace_na(fit_likelihood, 0))


#~ add autoencoder ----

# res_autoencoder <- qs::qread(file.path(dir_step2, "250515_autoencoder.qs"))
# colnames(res_autoencoder) <- paste0("ae_", 1:ncol(res_autoencoder))
# 
# stopifnot(all.equal(
#   predictors$cellgene,
#   rownames(res_autoencoder)
# ))
# 
# predictors <- cbind(predictors,
#                     res_autoencoder) |>
#   remove_rownames()


#as matrix

mat_predictors <- predictors |>
  select(-gene_name, -cell_type) |>
  column_to_rownames("cellgene") |>
  as.matrix()

head(mat_predictors)
mat_predictors[1:2,]


# check scales
par(mfrow = c(4,4), mar = c(3, 2, 2, 1) + 0.1)

for(i in seq_len(ncol(mat_predictors))){
  hist(mat_predictors[,i], main = colnames(mat_predictors)[[i]])
}
par(opar)


mat_sc <- apply(mat_predictors, 2, DescTools::Winsorize, na.rm = TRUE) |> scale()


# check scaling
par(mfrow = c(4,4), mar = c(3, 2, 2, 1) + 0.1)

for(i in seq_len(ncol(mat_sc))){
  hist(mat_sc[,i], main = colnames(mat_sc)[[i]])
}
par(opar)



stopifnot( !any(is.na(mat_sc)) )

# for comparisons
# qs::qsave(mat_sc,
#           file.path(dir_step2, "250516_mat_predictors_scaled.qs"))



# Correlated features? ----
cor(mat_sc) |> abs() |> pheatmap::pheatmap()
cor(mat_sc) |> pheatmap::pheatmap()

cor(mat_predictors) |> abs() |> pheatmap::pheatmap()
cor(mat_predictors) |> pheatmap::pheatmap()




# ## Using MSE or dev_expl as fit quality metric?
# c(mean( residuals(.mod)^2 )) in `mat`
# plot(mat_sc[,c("dev_expl", "mse")], cex= .1,
#      col = case_match(annot_pred_manual_amphso[rownames(mat_sc),],
#                       NA ~ "grey",
#                       "yes" ~ "green3",
#                       "no" ~ "red3",
#                       "maybe" ~ "blue3"))
# lm(dev_expl ~ mse, data = as.data.frame(mat_sc[,c("dev_expl", "mse")]))
# cor.test(mat_sc[,"dev_expl"], mat_sc[, "mse"])
# 
# as.data.frame(mat_sc[,c("dev_expl", "mse")]) |>
#   rownames_to_column("cellgene") |>
#   as_tibble() |>
#   filter(startsWith(cellgene, "AM_PHso")) |>
#   left_join(annot_pred_manual_amphso |> rownames_to_column("cellgene")) |>
#   ggplot() +
#   scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .1)) +
#   scale_color_manual(values = c(`NA` = 'grey', "yes" = "green3", "no" = "red3", "maybe" = "blue3")) +
#   geom_point(aes(x = dev_expl, y = mse,
#                  size = !is.na(manual),
#                  alpha = !is.na(manual),
#                  color = manual))

# ## Using AUC or intercept as fit quality metric?
# 
# plot(mat_sc[,c("auc", "intercept")], cex= .1,
#      col = case_match(annot_pred_manual_amphso[rownames(mat_sc),],
#                       NA ~ "grey",
#                       "yes" ~ "green3",
#                       "no" ~ "red3",
#                       "maybe" ~ "blue3"))
# lm(auc ~ intercept, data = as.data.frame(mat_sc[,c("auc", "intercept")]))
# cor.test(mat_sc[,"auc"], mat_sc[, "intercept"])
# 
# as.data.frame(mat_sc[,c("auc", "intercept")]) |>
#   rownames_to_column("cellgene") |>
#   as_tibble() |>
#   filter(startsWith(cellgene, "AM_PHso")) |>
#   left_join(annot_pred_manual_amphso |> rownames_to_column("cellgene")) |>
#   ggplot() +
#   scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .1)) +
#   scale_color_manual(values = c(`NA` = 'grey', "yes" = "green3", "no" = "red3", "maybe" = "blue3")) +
#   geom_point(aes(x = auc, y = intercept,
#                  size = !is.na(manual),
#                  alpha = !is.na(manual),
#                  color = manual))




#~ kmeans ----

set.seed(123)

# mat_sc1 <- mat_sc[ sample(nrow(mat_sc), .05*nrow(mat_sc)), subset_list[[cols]] ]
pca_res <- prcomp(mat_sc)
pcs_preds <- pca_res$x[,seq_len(4)]

km <- kmeans(pcs_preds, centers = 7, nstart = 100)

# table(km$cluster)

memberships <- km$cluster

acc(memberships)
acc3(memberships)
acc4(memberships)



# Check


cluster_labels <- enframe(memberships,
        name = "cell_gene",
        value = "cluster") |>
  separate_wider_delim(cell_gene,
                       delim = "|",
                       names = c("cell_type", "gene_name")) |>
  right_join(manual,
             by = join_by(cell_type, gene_name)) |>
  summarise(predicted_label = names(which.max(table(manual))),
            .by = cluster)

enframe(memberships,
        name = "cell_gene",
        value = "cluster") |>
  separate_wider_delim(cell_gene,
                       delim = "|",
                       names = c("cell_type", "gene_name")) |>
  filter(cluster %in% cluster_labels$cluster[cluster_labels$predicted_label == "yes"]) |>
  count(gene_name) |>
  pull(n) |>
  hist()


enframe(memberships,
        name = "gene_name",
        value = "cluster") |>
  filter(startsWith(gene_name, "AM_PHso") |
         startsWith(gene_name, "ILso")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "nas\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0) |>
  arrange(desc(guess))


enframe(memberships,
        name = "gene_name",
        value = "cluster") |>
  filter(startsWith(gene_name, "pharyngeal_muscle") |
           startsWith(gene_name, "pharynx_epithelial")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "abu\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0) |>
  arrange(desc(guess))




#~ PCA/UMAP ----

pca_res <- prcomp(mat_sc)
plot(pca_res$x[,c(1,2)])



pca_res$x |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  # filter(startsWith(cellgene, "AM_PHso")) |>
  left_join(annot_pred_manual_amphso |> rownames_to_column("cellgene")) |>
  ggplot() +
  theme_classic() +
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .1)) +
  scale_color_manual(values = c(`NA` = 'grey', "yes" = "green3", "no" = "red3", "maybe" = "blue3")) +
  geom_point(aes(x = PC1, y = PC2,
                 size = !is.na(manual),
                 alpha = !is.na(manual),
                 color = manual))


pca_res$x |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  mutate(cell_type = str_split_i(cellgene, "\\|", 1)) |>
  ggplot() +
  theme_classic() +
  ggsci::scale_color_d3(palette = "category20") +
  geom_point(aes(x = PC1, y = PC2,
                 color = cell_type),
             alpha = .2)



pca_res$x |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  left_join(enframe(memberships,
                    name = "cellgene",
                    value = "cluster"),
            by = "cellgene") |>
  mutate(cluster = as.factor(cluster)) |>
  ggplot() +
  theme_classic() +
  ggsci::scale_color_d3(palette = "category20") +
  geom_point(aes(x = PC1, y = PC2,
                 color = cluster),
             alpha = .6)





# sloooow on full matrix

# mat_sc1 <- mat_sc[startsWith(rownames(mat_sc), "AM_PHso"),]

um <- umap::umap(mat_sc1)
# qs::qsave()

um$layout |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  # filter(startsWith(cellgene, "AM_PHso")) |>
  left_join(annot_pred_manual_amphso |> rownames_to_column("cellgene")) |>
  ggplot() +
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .1)) +
  scale_color_manual(values = c(`NA` = 'grey', "yes" = "green3", "no" = "red3", "maybe" = "blue3")) +
  geom_point(aes(x = V1, y = V2,
                 size = !is.na(manual),
                 alpha = !is.na(manual),
                 color = manual)) +
  theme(legend.position = "none")




um$layout |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  mutate(guess = case_when(
    str_detect(cellgene, "col\\-[0-9]+") ~ "puls",
    str_detect(cellgene, "cutl\\-[0-9]+") ~ "puls",
    str_detect(cellgene, "grl\\-[0-9]+") ~ "puls",
    str_detect(cellgene, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(cellgene, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  # filter(startsWith(cellgene, "AM_PHso")) |>
  ggplot() +
  theme_classic() +
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .1)) +
  scale_color_manual(values = c(`NA` = 'grey', "puls" = "green3", "nonpuls" = "red3")) +
  geom_point(aes(x = V1, y = V2,
                 size = !is.na(guess),
                 alpha = !is.na(guess),
                 color = guess)) +
  theme(legend.position = "none")



um$layout |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  mutate(cell_type = str_split_i(cellgene, "\\|", 1)) |>
  # filter(startsWith(cellgene, "AM_PHso")) |>
  ggplot() +
  theme_classic() +
  ggsci::scale_color_d3(palette = "category20") +
  geom_point(aes(x = V1, y = V2,
                 color = cell_type),
             alpha = .8,
             size = .1)



um$layout |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  left_join(enframe(memberships,
                    name = "cellgene",
                    value = "cluster"),
            by = "cellgene") |>
  mutate(cluster = as.factor(cluster)) |>
  ggplot() +
  theme_classic() +
  ggsci::scale_color_d3(palette = "category20") +
  geom_point(aes(x = V1, y = V2,
                 color = cluster),
             alpha = .6)







#~ heatmaps ----

#~~ AM/PHso manually annotated only ----
table(rownames(annot_pred_manual_amphso)  %in% rownames(mat_sc))
mat_am_phso <- mat_sc[ rownames(mat_sc) %in% rownames(annot_pred_manual_amphso),  ]

# self-cluster (hclust)
pheatmap::pheatmap(t(mat_am_phso),
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   # scale = "row",
                   # show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_pred_manual_amphso)

membs1 <- memberships[names(memberships) %in% rownames(mat_am_phso)]
mat_am_phso_memberships <- mat_am_phso[names(membs1)[order(membs1)], ]
annot_pred_manual_amphso_membs <- annot_pred_manual_amphso
annot_pred_manual_amphso_membs$memb <- as.factor(membs1[rownames(annot_pred_manual_amphso_membs)])
pheatmap::pheatmap(t(mat_am_phso_memberships),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   # scale = "row",
                   # show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_pred_manual_amphso_membs,
                   annotation_colors = list(manual = c(yes = "green3",
                                                       no = "red3",
                                                       maybe = "blue3",
                                                       `NA` = "grey")))

smooth_amph <- smooth_centered[, rownames(annot_pred_manual_amphso_membs |> arrange(memb))]
pheatmap::pheatmap(log1p(smooth_amph),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   # scale = "row",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_pred_manual_amphso_membs,
                   annotation_colors = list(manual = c(yes = "green3",
                                                       no = "red3",
                                                       maybe = "blue3",
                                                       `NA` = "grey")))




#~~ AM/PHso, all genes ----

# only AM/PHso, all genes

mat_sc1 <- mat_sc[startsWith(rownames(mat_sc), "AM_PHso"),]

membs1 <- memberships[names(memberships) %in% rownames(mat_sc1)]
mat_sc1_ordered <- mat_sc1[names(membs1)[order(membs1)], ]


annot_amphso_all <- enframe(membs1,
                            name = "cellgene",
                            value = "cluster") |>
  mutate(cluster = as.factor(cluster)) |>
  left_join(pred_manual |>
              dplyr::select(manual, gene_name) |>
              mutate(gene_name = paste0("AM_PHso|", gene_name)),
            by = join_by(cellgene == gene_name)) |>
  column_to_rownames("cellgene")

pheatmap::pheatmap(t(mat_sc1_ordered),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   # scale = "row",
                   # show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_amphso_all,
                   annotation_colors = list(manual = c(yes = "green3",
                                                       no = "red3",
                                                       maybe = "blue3",
                                                       `NA` = "grey")))

smooth_amph <- smooth_centered[, rownames(annot_amphso_all |> arrange(cluster))]
pheatmap::pheatmap(log1p(smooth_amph),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   # scale = "row",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = annot_amphso_all,
                   annotation_colors = list(manual = c(yes = "green3",
                                                       no = "red3",
                                                       maybe = "blue3",
                                                       `NA` = "grey")))






#~ plot genes ----

# # plot examples from some clusters
# par(mfrow = c(2,2), mar = c(3, 2, 2, 1) + 0.1)
# for(i in c(1,2,8,11)){
#   exple_cl <- sample(names(memberships)[memberships == i], 10)
#   matplot(log1p(all_preds_centered)[,exple_cl], type = "l")
# }
# par(opar)
# 
# 
# # plot average per cluster
# for(cl in unique(memberships)){
#   # cl = cl+1
#   
#   genes_in_clust <- names(memberships)[memberships == cl]
#   
#   average_signal <- log1p(all_preds_centered)[,genes_in_clust] |> rowMeans()
#   signal_sd <- log1p(all_preds_centered)[,genes_in_clust] |> matrixStats::rowSds()
#   
#   plot(average_signal, type = "l", lwd = 1.5, ylim = c(0, max(3, 1.1*max(average_signal))))
#   lines(average_signal - signal_sd , lwd = .5, lty = 'dashed')
#   lines(average_signal + signal_sd , lwd = .5, lty = 'dashed')
# }

#~ plot average ----
all_clustered_fits <- log1p(smooth_centered) |>
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

# all_clustered_fits |>
#   summarize(average_signal = mean(log_cnt),
#             sd_signal = sd(log_cnt),
#             .by = c(cluster, time)) |>
#   ggplot() +
#   theme_classic() +
#   facet_wrap(~cluster) +
#   geom_hline(aes(yintercept = 0),
#              linetype = 'dashed', color = 'grey80') +
#   geom_ribbon(aes(
#     x = time,
#     ymin = average_signal - sd_signal,
#     ymax = average_signal + sd_signal
#   ),
#   alpha = .2) +
#   geom_line(aes(x = time, y = average_signal))


selected <- all_clustered_fits |>
  select(cluster, cell_type, gene_name) |>
  distinct() |>
  group_by(cluster) |>
  slice_sample(n = 30) |>
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
    c(1, 2, 4, 5) ~ "nonpulsatile",
    c(3,6)~ "pulsatile"
    ))

# cluster_results |>
#   write_csv(file.path(dir_step2, "250424_cluster_results.csv"))




cluster_results |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(prop_puls)


# par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)
# 
# 
# exple_cl1 <- sample(names(memberships)[memberships == 1], 10)
# matplot(log1p(all_preds)[,exple_cl1], type = "l")
# 
# exple_cl2 <- sample(names(memberships)[memberships == 3], 10)
# matplot(log1p(all_preds)[,exple_cl2], type = "l")
# par(opar)



# Compare clusterings ----

cluster_results <- read_csv(file.path(dir_step2, "250424_cluster_results.csv")) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(cellgene %in% names(memberships))







enframe(memberships,
        name = "gene_name",
        value = "cluster") |>
  filter(startsWith(gene_name, "AM_PHso")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0)


cluster_results |>
  filter(startsWith(cellgene, "AM_PHso")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0)



# Other clustering algorithms ----


#~ hclust ----


hc <- mat_sc |>
  dist() |>
  fastcluster::hclust()
hccl <- cutree(hc, k = 4)


table(hccl)


enframe(hccl,
        name = "gene_name",
        value = "cluster") |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0)


#~ Leiden ----


get_knn_graph <- function(mat, k){
  
  # build Annoy index
  ann <- new(RcppAnnoy::AnnoyEuclidean, ncol(mat))
  
  for(i in 1:nrow(mat)) ann$addItem(i - 1, mat[i,]) # RcppAnnoy 0-indexed
  
  ann$build(50)
  
  # get knns
  ann.idx <- lapply(seq_len(nrow(mat)),
                    \(i){
                      res <- ann$getNNsByVectorList(mat[i,], k, -1, TRUE)
                      # discard res$distance, only use res$item (+1 since 0-indexed)
                      res$item + 1L
                    })
  
  nns.idx <- do.call(rbind, ann.idx)
  
  # we have a matrix where each row = source (n rows), k columns = targets (k cols)
  # we transform to edglist with one row for each pair (i.e. n*k rows)
  row_idx <- rep(seq_len(nrow(nns.idx)), times = ncol(nns.idx))
  col_idx <- as.vector(nns.idx)
  
  # keep only one direction and remove any self-loop
  valid_edges <- row_idx < col_idx
  
  edge_matrix <- cbind(row_idx[valid_edges], col_idx[valid_edges])
  
  igraph::graph_from_edgelist(edge_matrix, directed = FALSE)
  
}

# singleton detection adapted from Seurat
remove_singletons <- function(ids, SNN){
  
  singletons <- names(which(table(ids) == 1))
  cluster_names <- as.character(unique(ids)) |> setdiff(singletons)
  for (i in singletons) {
    i.cells <- which(ids == i)
    
    connectivity <- vapply(cluster_names,
                           \(cl){
                             
                             subSNN <- SNN[i.cells, which(ids == cl)]
                             
                             mean(subSNN)
                           },
                           double(1L))
    
    
    m <- max(connectivity, na.rm = T)
    mi <- which(connectivity == m, arr.ind = TRUE)
    closest_cluster <- names(connectivity[mi]) |> sample(1)
    ids[i.cells] <- closest_cluster
  }
  ids
}

gr <- get_knn_graph(mat_sc, 30)




ld_cl <- leidenbase::leiden_find_partition(gr,
                                           partition_type = "RBConfigurationVertexPartition",
                                           resolution_parameter = 5e-1) |>
  chuck("membership") |>
  remove_singletons(SNN = igraph::as_adjacency_matrix(gr))





length(unique(ld_cl))


table(ld_cl)



tibble(gene_name = rownames(mat_sc),
       cluster = ld_cl) |>
  filter(startsWith(gene_name, "ILso")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0)


memberships <- ld_cl |> set_names(rownames(mat_sc))



pca <- prcomp(mat_sc)
pca$x |>
  as.data.frame() |>
  rownames_to_column("cellgene") |>
  as_tibble() |>
  add_column(cluster = as.factor(memberships)) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PC1, y = PC2, color = cluster),
             alpha = .2)



#~ louvain ----
adj_mat <- igraph::as_adjacency_matrix(gr)
rownames(adj_mat) <- colnames(adj_mat) <- rownames(mat_sc)

cl_louv <- Seurat:::RunModularityClustering(SNN = Seurat::as.Graph(adj_mat), 
                                            resolution = .8) |>
  remove_singletons(SNN = igraph::as_adjacency_matrix(gr))

length(unique(cl_louv))

table(cl_louv)



tibble(gene_name = rownames(mat_sc),
       cluster = cl_louv) |>
  # filter(startsWith(gene_name, "AM_PHso")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0)


memberships <- cl_louv |> set_names(rownames(mat_sc))




#~ kohonen ----

som <- kohonen::som(mat_sc1, grid = kohonen::somgrid(6,6,topo = "hexagonal"))

som$unit.classif |> table()
# 
# cl <- cutree(hclust(dist(som$codes[[1]])), k = 5)
# 
# plot(som, type = "mapping", bgcol = rainbow(5)[cl])
# kohonen::add.cluster.boundaries(som, cl)

cl <- kmeans(som$codes[[1]], centers = 3)
plot(som, type = "mapping", bgcol = rainbow(5)[cl$cluster])
kohonen::add.cluster.boundaries(som, cl$cluster)


memberships <- som$unit.classif |> set_names(rownames(mat_sc1))
# memberships <- cl$cluster[som$unit.classif] |> set_names(rownames(mat_sc))

table(memberships)





all_clustered_fits |>
  inner_join(selected) |>
  ggplot() +
  theme_classic() +
  scale_color_brewer(type = "qual", palette = "Set2") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  facet_wrap(~cluster,nrow = 6, ncol = 6) +
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
    linewidth = 1,
    color = "orange2",
    data = all_clustered_fits |>
      summarize(average_signal = mean(log_cnt),
                sd_signal = sd(log_cnt),
                .by = c(cluster, time))
  )






guesses <- enframe(memberships,
        name = "gene_name",
        value = "cluster") |>
  # filter() |>
  mutate(guess = case_when(
    (startsWith(gene_name, "AM_PHso") |
      startsWith(gene_name, "ILso")) & str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    (startsWith(gene_name, "AM_PHso") |
       startsWith(gene_name, "ILso")) & str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    (startsWith(gene_name, "AM_PHso") |
       startsWith(gene_name, "ILso")) & str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    (startsWith(gene_name, "AM_PHso") |
       startsWith(gene_name, "ILso")) & str_detect(gene_name, "nas\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls",
    (startsWith(gene_name, "pharyngeal_muscle") |
       startsWith(gene_name, "pharynx_epithelial")) &str_detect(gene_name, "abu\\-[0-9]+") ~ "puls"
  )) |>
  mutate(color = case_match(guess,
                            "puls" ~ "green3",
                            "nonpuls" ~ "red3",
                            .default = alpha("grey", .01)))

plot(som, type = "mapping", col = guesses$color)


#manually annotated
guesses <- enframe(memberships,
                   name = "gene_name",
                   value = "cluster") |>
  separate_wider_delim(gene_name, delim = "|",
                       names = c("cell_type","gene_name")) |>
  left_join(manual_several) |>
  mutate(color = case_match(manual,
                            "yes" ~ "green3",
                            "no" ~ "red3",
                            .default = alpha("grey", .01)))

plot(som, type = "mapping", col = guesses$color)
text(som$grid$pts[,1], som$grid$pts[,2],
     labels = 1:nrow(som$grid$pts),
     adj = c(0,0))


guesses |>
  filter(!is.na(manual),
         cluster == 24)

plot(som, type = "mapping", col = 'grey90')
text(som$grid$pts[,1], som$grid$pts[,2], labels = 1:nrow(som$grid$pts))


guesses |>
  filter(!is.na(guess),
         cluster == 17)

enframe(memberships,
        name = "gene_name",
        value = "cluster") |>
  filter(startsWith(gene_name, "pharyngeal_muscle") |
           startsWith(gene_name, "pharynx_epithelial")) |>
  mutate(guess = case_when(
    str_detect(gene_name, "abu\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n",
              values_fill = 0) |>
  arrange(desc(guess))





# _____ ----


#~ Embedding AE ----


library(keras)
len <- nrow(smooth_centered)

clusters <- map(1:10,
                ~{
                  
                  genes <- colnames(smooth_centered) |> sample(10000)
                  
                  x_train <- smooth_centered[, sample(genes, .8*length(genes))] |> t()
                  x_test <- smooth_centered[, setdiff(genes, rownames(x_train))] |> t()
                  
                  x_train <- x_train
                  x_test <- x_test
                  
                  stopifnot(all.equal(
                    genes |> sort(),
                    union(rownames(x_train), rownames(x_test)) |> sort()
                  ))
                  
                  # dim(x_train);dim(x_test)
                  
                  # matplot(x_train[sample(nrow(x_train), 5),] |> t(), type = "l")
                  
                  
                  
                  
                  n_bottleneck <- 8
                  
                  # Define the encoder
                  encoder <- keras_model_sequential() |>
                    layer_dense(units = 256, input_shape = len, activation = 'relu') |>
                    layer_dense(units = n_bottleneck, activation = NULL, name = "bottleneck") |>
                    layer_activation_leaky_relu(alpha = .6)
                  
                  # Define the decoder
                  decoder <- keras_model_sequential() |>
                    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
                    layer_dense(units = len, activation = "linear")
                  
                  # Connect them to create the autoencoder
                  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
                  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
                  
                  
                  fitted <- autoencoder |>
                    fit(x_train,
                        x_train,
                        epochs = 10,
                        batch_size = 16,
                        validation_data = list(x_test, x_test),
                        verbose = FALSE)
                  
                  
                  # extract the bottleneck layer
                  intermediate_layer_model <- keras_model(inputs = autoencoder$input,
                                                          outputs = get_layer(autoencoder, "bottleneck")$output)
                  intermediate_output <- predict(intermediate_layer_model, t(smooth_centered))
                  
                  kmeans(apply(intermediate_output, 1, \(x) x/max(x)) |> t(), centers = 5)
                },
                .progress = TRUE)



cons <- clue::cl_consensus(clusters,
                           method = "HE")[[".Data"]]



clusters <- tibble(
  names = colnames(smooth_centered),
  cluster = apply(cons, 1, function(row_vals) which(row_vals == 1))
) |>
  separate_wider_delim(names,
                       delim = "|",
                       names = c("cell_type", "gene_name"))

clusters |>
  mutate(guess = case_when(
    str_detect(gene_name, "^col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "^cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "^grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "^nas\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "^rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "^rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess) |>
  pivot_wider(id_cols = guess, names_from = "cluster", values_from = "n")





# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = autoencoder$input,
                                        outputs = get_layer(autoencoder, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, t(smooth_centered))

reconstructed <- predict(autoencoder, t(smooth_centered))


rownames(intermediate_output) <- rownames(reconstructed) <- colnames(smooth_centered)

# intermediate_output <- intermediate_output[rownames(intermediate_output) %in% paste0("AM_PHso|", pred_manual$gene_name),]
# reconstructed <- reconstructed[rownames(reconstructed) %in% paste0("AM_PHso|", pred_manual$gene_name),]
# preds_cent <- smooth_centered[,colnames(smooth_centered) %in% paste0("AM_PHso|", pred_manual$gene_name)]

# opar <- par(no.readonly = TRUE)

par(mfrow = c(1,3), mar = c(3, 2, 2, 1) + 0.1)
layout(matrix(c(1,3,2,3), nrow = 2))
exple <- sample(rownames(x_train), 5)


matplot(t(x_train)[,exple], type = "l")
matplot(t(reconstructed)[,exple], type = "l")
matplot(t(intermediate_output)[,exple], type = "p")
par(opar)


hc <- hclust(dist(apply(intermediate_output, 1, \(x) x/max(x)) |> t()))

pheatmap::pheatmap(smooth_centered,
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # scale = "column",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = pred_manual |>
                     dplyr::select(manual, gene_name) |>
                     mutate(gene_name = paste0("AM_PHso_", gene_name)) |>
                     column_to_rownames("gene_name"))

# pheatmap::pheatmap(apply(intermediate_output, 1, \(x) x/max(x)),
#                    cluster_rows = FALSE,
#                    cluster_cols = hc,
#                    # scale = "column",
#                    annotation_col = pred_manual |>
#                      dplyr::select(manual, gene_name) |>
#                      column_to_rownames("gene_name"))

# cl <- cluster::pam(apply(intermediate_output, 1, \(x) x/max(x)) |> t(), k = 2)
# cl$clustering |> table(useNA = 'ifany')

km <- kmeans(apply(intermediate_output, 1, \(x) x/max(x)) |> t(), centers = 2)

km$cluster |> table(useNA = 'ifany')

enframe(km$cluster,
        name = "gene_name",
        value = "cluster") |>
  mutate(guess = case_when(
    str_detect(gene_name, "_col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "_cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "_grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "_nas\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "_rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "_rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  filter(!is.na(guess)) |>
  count(cluster, guess)


res <- pred_manual |>
  mutate(gene_name = paste0("AM_PHso_", gene_name)) |>
  left_join(
    tibble(gene_name = names(km$cluster),
           cluster = km$cluster),
    by = "gene_name"
  ) |>
  arrange(desc(manual), cluster)

table(res$manual, res$cluster)


# examples from each cluster

par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)

exple_cl1 <- sample(names(km$cluster)[km$cluster == 1], 10)
matplot(log1p(all_preds_centered)[,exple_cl1], type = "l")

exple_cl2 <- sample(names(km$cluster)[km$cluster == 2], 10)
matplot(log1p(all_preds_centered)[,exple_cl2], type = "l")
par(opar)


# loss based on manual

res |>
  filter(manual %in% c("no", "yes")) |>
  count(manual, cluster)




#~ previous CNN model ----



# autoencoder
n_bottleneck <- 8

# Define the encoder
encoder <- keras_model_sequential() |>
  layer_conv_1d(filters = 64, kernel_size = 7, activation = 'relu', input_shape = c(len, 1),
                name = "inputConv") |>
  layer_max_pooling_1d(pool_size = 2,
                       name = "Pool") |>
  # layer_conv_1d(filters = 64, kernel_size = 3, activation = 'relu') |>
  # layer_max_pooling_1d(pool_size = 2) |>
  layer_flatten() |>
  layer_dense(units = n_bottleneck, activation = 'relu', name = "bottleneck")

# Define the decoder
decoder <- keras_model_sequential() |>
  layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck,
              name = "denseFromBottleneck") |>
  layer_reshape(target_shape = c(64, 1),
                name = "reshape") |>
  layer_conv_1d_transpose(filters = 64, kernel_size = 7,
                          strides = 2, activation = 'relu',
                          padding = 'same',
                          name = "deconv") |>
  # layer_conv_1d_transpose(filters = 32, kernel_size = 3,
  #                         strides = 2, activation = 'relu',
  #                         padding = 'same') |>
  layer_conv_1d(filters = 1, kernel_size = 3,
                activation = 'linear', padding = 'same',
                name = "output")

# Connect them to create the autoencoder
autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')







#~ clust ts dist ----

library(dtwclust)
dist_cort <- manual_dist_CORT(t(all_preds_centered), k = 0.1)

dist_dtw <- dtw::dtwDist(t(all_preds_centered))


proxy::dist()



dist_jeffreys <- philentropy::dist_many_many(t(log1p(all_preds_centered)),
                                             t(log1p(all_preds_centered)),
                                             method = "jeffreys")
colnames(dist_jeffreys) <- rownames(dist_jeffreys) <- colnames(all_preds_centered)

# kl_log <- sign(kl_dist_int) * log(abs(kl_dist_int))
# diag(kl_log) <- 0
# dist_kl <- ( kl_log - min(kl_log) )/max( kl_log - min(kl_log) )

hist(dist_jeffreys)

tsc <- dtwclust::tsclust(t(log1p(all_preds_centered)),
                         type = "partition",
                         control = dtwclust::partitional_control(
                           distmat = dist_cort
                         ),
                         k = 10L)



(tsc |> plot(type = "centroids", labels = list())) + ggplot2::facet_wrap(~cl, scales = "free_x")


res <- pred_manual |>
  left_join(
    tibble(gene_name = colnames(all_preds_centered),
           cluster = tsc@cluster),
    by = "gene_name"
  ) |>
  arrange(desc(manual), cluster)



table(res$manual, res$cluster)



View(res)

dtwclust::cvi(tsc, type = c("CH","SF"))


#~ compare clusterings ----
cfgs <- compare_clusterings_configs(
  types = c("p"),
  k = 4:10,
  controls = list(
    partitional = partitional_control(
      iter.max = 30L,
      nrep = 1L
    )
  ),
  centroids = pdc_configs(
    type = "centroid",
    partitional = list(
      pam = list()
    )
  )
)





cmp <- dtwclust::compare_clusterings(t(log1p(all_preds_centered)),
                                     types = "p",
                                     configs = cfgs,
                                     score.clus = cvi_evaluators(type = c("SF","CH"))$score
)

cmp

cvi(tsc, type = c("CH","SF"))


# wavelets ----



wavelet_trans <- apply(all_preds_centered, 2, \(ts) {
  
  
  wt_res <- waved::WaveD(ts, MC = TRUE, SOFT = TRUE)
  wt_res
  
  # dwt_result <- wavelets::dwt(ts, filter = "c24")
  # unlist(dwt_result@W)
  
})

xx <- scale(t(wavelet_trans))

kmeans_result <- kmeans(xx, centers = 5L)
table(kmeans_result$cluster)


res <- pred_manual |>
  left_join(
    tibble(gene_name = names(kmeans_result$cluster),
           cluster = kmeans_result$cluster),
    by = "gene_name"
  ) |>
  arrange(desc(manual), cluster)



table(res$manual, res$cluster)
###


















#### Diff centering? ----

# does centering before fitting make a difference?

res_noncent <- qs::qread("intermediates/2502/250422_step2_nb/AM_PHso_res_gam.qs")
res_cent <- qs::qread("intermediates/2502/250422_step2_nb_centered/AM_PHso_res_gam.qs")



mat_noncent <- vapply(res_noncent$gam_fit,
                      \(.mod){
                        coef(.mod)
                      },
                      double(5L)) |> t()
colnames(mat_noncent) <- c("intercept", "s1", "s2", "s3", "s4")
mat_noncent <- cbind(
  mat_noncent,
  var_s = matrixStats::rowVars(mat_noncent[,c("s1","s2","s3", "s4")]),
  max_s = matrixStats::rowMaxs(mat_noncent[,c("s1","s2","s3", "s4")])
)
mat_noncent <- mat_noncent[,c("intercept", "var_s", "max_s")]


mat_cent <- vapply(res_cent$gam_fit,
                   \(.mod){
                     coef(.mod)
                   },
                   double(5L)) |> t()
colnames(mat_cent) <- c("intercept", "s1", "s2", "s3", "s4")
mat_cent <- cbind(
  mat_cent,
  var_s = matrixStats::rowVars(mat_cent[,c("s1","s2","s3", "s4")]),
  max_s = matrixStats::rowMaxs(mat_cent[,c("s1","s2","s3", "s4")])
)
mat_cent <- mat_cent[,c("intercept", "var_s", "max_s")]


printMat::matimage(mat_noncent)
printMat::matimage(mat_cent)

all.equal(rownames(mat_cent),
          rownames(mat_noncent))

ex_rows <- rownames(mat_cent)[startsWith(rownames(mat_cent), "cutl") |
                                startsWith(rownames(mat_cent), "col") |
                                startsWith(rownames(mat_cent), "grl") |
                                startsWith(rownames(mat_cent), "rps") |
                                startsWith(rownames(mat_cent), "rpl")]
annot <- tibble(gene_name = ex_rows) |>
  mutate(guess = case_when(
    str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
    str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
    str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  )) |>
  column_to_rownames("gene_name")

pheatmap::pheatmap(apply(mat_noncent[ex_rows,], 2, DescTools::Winsorize) |> scale(),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   fontsize_row = 5,
                   main = "noncent",
                   annotation_row = annot)
pheatmap::pheatmap(apply(mat_cent[ex_rows,], 2, DescTools::Winsorize) |> scale(),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   fontsize_row = 5,
                   annotation_row = annot,
                   main = "centered")





