# Inits ----

library(tidyverse)




dir_step2 <- "intermediates/2502/250514_step2"

dir_out <- "intermediates/2502/250512_clustering"

mat_sc <- qs::qread(file.path(dir_step2, "250515_mat_predictors_scaled.qs"))



manual <- bind_rows(
  # readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx") |>
  #   select(cell_type, gene_name, manual),
  # readxl::read_excel("intermediates/2502/250330_step2/manual_annotation.xlsx") |>
  #   dplyr::select(-amplitude, -dev_expl),
  readxl::read_excel("data/manual_annotations/250507_manual.xlsx") |>
    select(2:4)
) |>
  filter(manual %in% c("yes", "no"))


# Functions ----

# #~ accuracy ----
# 
# # based on manual
# acc <- function(membs){
#   comp_manual <- enframe(membs,
#                          name = "cell_gene",
#                          value = "cluster") |>
#     separate_wider_delim(cell_gene,
#                          delim = "|",
#                          names = c("cell_type", "gene_name")) |>
#     right_join(manual,
#                by = join_by(cell_type, gene_name))
#   
#   
#   
#   cluster_labels <- comp_manual |>
#     summarise(predicted_label = names(which.max(table(manual))),
#               .by = cluster)
#   
#   comp_manual <- comp_manual |>
#     left_join(cluster_labels,
#               by = join_by(cluster))
#   
#   # accuracy
#   round( 100 * mean(comp_manual$manual == comp_manual$predicted_label) )
# }
# 
# 
# 
# # based on gene name (guess)
# 
# acc_guess <- function(membs){
#   comp_manual <- enframe(membs,
#                          name = "cell_gene",
#                          value = "cluster") |>
#     separate_wider_delim(cell_gene,
#                          delim = "|",
#                          names = c("cell_type", "gene_name")) |>
#     mutate(guess = case_when(
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "col-") ~ "puls",
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "cutl-") ~ "puls",
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "grl-") ~ "puls",
#       startsWith(gene_name, "rps-") ~ "nonpuls",
#       startsWith(gene_name, "rpl-") ~ "nonpuls",
#       (startsWith(cell_type, "pharyngeal_muscle") |
#          startsWith(cell_type, "pharynx_epithelial")) &startsWith(gene_name, "abu-") ~ "puls"
#     )) |>
#     filter(!is.na(guess))
#   
#   # lapply(1:4,
#   #        \(.cl) {comp_manual |> filter(cluster == .cl) |> pull(guess) |> table()})
#   
#   cluster_labels <- comp_manual |>
#     summarise(predicted_label = names(which.max(table(guess))),
#               .by = cluster)
#   
#   comp_manual <- comp_manual |>
#     left_join(cluster_labels,
#               by = join_by(cluster))
#   
#   # comp_manual |>
#   #   count(guess, predicted_label)
#   
#   # accuracy
#   round( 100 * mean(comp_manual$guess == comp_manual$predicted_label) )
# }
# 
# 
# # use gene families to select good/bad clusters, check the accuracy of manual "yes"
# acc3 <- function(membs){
#   cluster_labels <- enframe(membs,
#                          name = "cell_gene",
#                          value = "cluster") |>
#     separate_wider_delim(cell_gene,
#                          delim = "|",
#                          names = c("cell_type", "gene_name")) |>
#     mutate(guess = case_when(
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "col-") ~ "puls",
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "cutl-") ~ "puls",
#       (startsWith(cell_type, "AM_PHso") |
#          startsWith(cell_type, "ILso")) & startsWith(gene_name, "grl-") ~ "puls",
#       startsWith(gene_name, "rps-") ~ "nonpuls",
#       startsWith(gene_name, "rpl-") ~ "nonpuls",
#       (startsWith(cell_type, "pharyngeal_muscle") |
#          startsWith(cell_type, "pharynx_epithelial")) &startsWith(gene_name, "abu-") ~ "puls"
#     )) |>
#     filter(!is.na(guess)) |>
#     summarise(predicted_label = names(which.max(table(guess))),
#               .by = cluster)
#   
#   
#   comp_manual <- enframe(membs,
#                          name = "cell_gene",
#                          value = "cluster") |>
#     separate_wider_delim(cell_gene,
#                          delim = "|",
#                          names = c("cell_type", "gene_name")) |>
#     inner_join(manual,
#                by = join_by(cell_type, gene_name)) |>
#     left_join(cluster_labels,
#               by = join_by(cluster))
#   
#   
#   
#   # comp_manual |>
#   #   count(manual, predicted_label)
#   
#   # accuracy
#   round( 100 * sum(comp_manual$manual == "yes" & comp_manual$predicted_label == "puls") / sum(comp_manual$manual == "yes") )
# }

# split manual in train/test
acc4 <- function(membs){
  
  comp_manual <- enframe(membs,
                         name = "cell_gene",
                         value = "cluster") |>
    separate_wider_delim(cell_gene,
                         delim = "|",
                         names = c("cell_type", "gene_name")) |>
    inner_join(manual,
               by = join_by(cell_type, gene_name))
  
  replicate(20,
            {
              train <- sample(nrow(comp_manual), .5*nrow(comp_manual))
              test <- seq_len(nrow(comp_manual)) |> setdiff(train)
              
              # table(comp_manual[train,]$manual, comp_manual[train,]$cluster)
              
              cluster_labels <- comp_manual[train,] |>
                summarise(predicted_label = names(which.max(table(manual))),
                          .by = cluster)
              
              
              test_manual <- comp_manual[test,] |>
                left_join(cluster_labels,
                          by = join_by(cluster))
              
              # # accuracy (total)
              # round( 100 * mean(test_manual$manual == test_manual$predicted_label) )
              
              # test_manual |>
              #   count(manual, predicted_label)
              
              # accuracy (yes only)
              round( 100 * sum(test_manual$manual == "yes" & test_manual$predicted_label == "yes") / sum(test_manual$manual == "yes") )
            }) |>
    mean(na.rm = TRUE)
  
}

run_once_algos <- function(cols, method, nb_clust, nb_pcs, replicate){
  
  mat_sc1 <- mat_sc[ sample(nrow(mat_sc), .1*nrow(mat_sc)),
                     subset_list[[cols]] ]
  
  nb_pcs <- min(nb_pcs, ncol(mat_sc1))
  
  if(method == "pca_km"){
    
    pca_res <- prcomp(mat_sc1)
    pcs_preds <- pca_res$x[,seq_len(nb_pcs)]
    
    km <- kmeans(pcs_preds, centers = nb_clust, nstart = 2)
    
    memberships <- km$cluster
    
    memberships <- km$cluster
    wss <- km$tot.withinss
    centroids <- km$centers
    n <- nrow(pcs_preds)
    p <- ncol(pcs_preds)
    
    # sil <- cluster::silhouette(
    #   memberships,
    #   stats::dist(pcs_preds)
    # )
    
  } else if(method == "cor_pam"){
    
    d <- as.dist(1 - cor(t(mat_sc1)))
    
    pam <- cluster::pam(d, k = nb_clust)
    
    memberships <- pam$clustering
    
    
    # wss
    n <- nrow(mat_sc1)
    p <- ncol(mat_sc1)
    centroids <- matrix(0, nrow = nb_clust, ncol = p)
    for(i in 1:nb_clust) {
      if(sum(memberships == i) > 0) {
        centroids[i,] <- colMeans(mat_sc1[memberships == i, , drop = FALSE])
      }
    }
    
    # Calculate WSS in original feature space
    wss <- 0
    for(i in 1:n) {
      cluster_i <- memberships[i]
      wss <- wss + sum((mat_sc1[i,] - centroids[cluster_i,])^2)
    }
    
    
    # sil <- cluster::silhouette(
    #   memberships,
    #   d
    # )
    
  } else if(method == "cor_hc"){
    
    d <- as.dist(1 - cor(t(mat_sc1)))
    hc <- fastcluster::hclust(d)
    
    memberships <- stats::cutree(hc, k = nb_clust)
    
    
    
    # For hierarchical clustering, compute centroids and WSS manually
    n <- nrow(mat_sc1)
    p <- ncol(mat_sc1)
    
    # Calculate cluster centroids
    centroids <- matrix(0, nrow = nb_clust, ncol = p)
    for(i in 1:nb_clust) {
      if(sum(memberships == i) > 0) {
        centroids[i,] <- colMeans(mat_sc1[memberships == i, , drop = FALSE])
      }
    }
    
    # Calculate WSS
    wss <- 0
    for(i in 1:n) {
      cluster_i <- memberships[i]
      wss <- wss + sum((mat_sc1[i,] - centroids[cluster_i,])^2)
    }
    
    # sil <- cluster::silhouette(
    #   memberships,
    #   d
    # )
    
  } else{
    stop("Invalid method")
  }
  
  
  dunn <- clValid::dunn(clusters = memberships, Data = mat_sc1)
  
  
  m <- ifelse(method == "pca_km", nb_pcs, ncol(mat_sc1))
  k <- m * nb_clust + nb_clust
  
  bic <- k * log(n) + n * log(wss/n)
  aic <- 2 * k + n * log(wss/n)
  
  data.frame(cols = cols,
             method = method,
             nb_clust = nb_clust,
             nb_pcs = nb_pcs,
             replicate = replicate,
             accur = acc4(memberships),
             # silhouette = mean(sil[, 3]),
             dunn = dunn,
             aic = aic,
             bic = bic)
}


# several kmeans
run_once_km <- function(cols, method, nb_clust, nb_comp, replicate){
  
  mat_sc1 <- mat_sc[ sample(nrow(mat_sc), .1*nrow(mat_sc)),
                     subset_list[[cols]] ]
  
  
  if(method == "pca_km"){
    
    
    nb_pcs <- min(nb_comp, ncol(mat_sc1))
    
    pca_res <- prcomp(mat_sc1)
    pcs_preds <- pca_res$x[,seq_len(nb_pcs)]
    
    km <- kmeans(pcs_preds, centers = nb_clust, nstart = 2)
    
    
    memberships <- km$cluster
    
    # for AIC
    reduced_dim <- nb_pcs
    
    
    
  } else if(method == "som_km"){
    
    
    som <- kohonen::som(mat_sc1,
                        grid = kohonen::somgrid(nb_comp,
                                                nb_comp,
                                                topo = "hexagonal"))
    
    
    km <- kmeans(som$codes[[1]], centers = nb_clust, nstart = 2)
    # plot(som, type = "mapping", bgcol = rainbow(5)[km$cluster])
    # kohonen::add.cluster.boundaries(som, km$cluster)
    
    
    # memberships <- som$unit.classif |> set_names(rownames(mat_sc1))
    memberships <- km$cluster[som$unit.classif] |> set_names(rownames(mat_sc1))
    
    # table(memberships)
    
    
    # For AIC
    reduced_dim <- nrow(som$grid$pts)
    
    
  } else{
    stop("Invalid method")
  }
  
  
  dunn <- clValid::dunn(clusters = memberships, Data = mat_sc1)
  
  
  # AIC
  n <- nrow(mat_sc1)
  p <- ncol(mat_sc1)
  
  centroids_orig <- seq_len(nb_clust) |>
    lapply(
      \(.i) which(memberships == .i)
    ) |>
    vapply(
      \(.idx) colMeans(mat_sc1[.idx, , drop = FALSE]),
      FUN.VALUE = double(p)
    ) |>
    t()
  
  wss <- seq_len(n) |>
    vapply(
      \(.i) ( mat_sc1[.i,] - centroids_orig[ memberships[[.i]], ] )^2,
      FUN.VALUE = double(p)
    ) |>
    sum()
  k <- reduced_dim * nb_clust + nb_clust
  
  bic <- k * log(n) + n * log(wss/n)
  aic <- 2 * k + n * log(wss/n)
  
  data.frame(cols = cols,
             method = method,
             nb_clust = nb_clust,
             nb_comp = nb_comp,
             replicate = replicate,
             accur = acc4(memberships),
             # silhouette = mean(sil[, 3]),
             dunn = dunn,
             aic = aic,
             bic = bic)
}

# run_once <- function(cols, method, nb_clust, nb_pcs, replicate){
#   
#   mat_sc1 <- mat_sc[ sample(nrow(mat_sc), .1*nrow(mat_sc)),
#                      subset_list[[cols]] ]
#   
#   nb_pcs <- min(nb_pcs, ncol(mat_sc1))
#   
#   if(method == "pca_km"){
#     
#     pca_res <- prcomp(mat_sc1)
#     pcs_preds <- pca_res$x[,seq_len(nb_pcs)]
#     
#     km <- kmeans(pcs_preds, centers = nb_clust, nstart = 2)
#     
#     memberships <- km$cluster
#     
#     
#     
#     # sil <- cluster::silhouette(
#     #   memberships,
#     #   stats::dist(pcs_preds)
#     # )
#     
#   } else if(method == "cor_pam"){
#     
#     d <- as.dist(1 - cor(t(mat_sc1)))
#     
#     pam <- cluster::pam(d, k = nb_clust)
#     
#     memberships <- pam$clustering
#     
#     # sil <- cluster::silhouette(
#     #   memberships,
#     #   d
#     # )
#     
#   } else if(method == "cor_hc"){
#     
#     d <- as.dist(1 - cor(t(mat_sc1)))
#     hc <- fastcluster::hclust(d)
#     
#     memberships <- stats::cutree(hc, k = nb_clust)
#     
#     # sil <- cluster::silhouette(
#     #   memberships,
#     #   d
#     # )
#     
#   } else{
#     stop("Invalid method")
#   }
#   
#   
#   dunn <- clValid::dunn(clusters = memberships, Data = mat_sc1)
#   
#   data.frame(cols = cols,
#              method = method,
#              nb_clust = nb_clust,
#              nb_pcs = nb_pcs,
#              replicate = replicate,
#              accur = acc4(memberships),
#              # silhouette = mean(sil[, 3]),
#              dunn = dunn)
# }




# Configurations ----


all_cols <- colnames(mat_sc)

# all_cols |> paste0(collapse = "', '")

subset_list <- list(
  all =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'area_under_curve', 'dist_dtw', 'fit_likelihood','ae_1', 'ae_2', 'ae_3', 'ae_4'),
  no_fit =   c(                                                  'amplitude', 'area_under_curve', 'dist_dtw', 'fit_likelihood','ae_1', 'ae_2', 'ae_3', 'ae_4'),
  no_curve = c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl',                                  'dist_dtw', 'fit_likelihood','ae_1', 'ae_2', 'ae_3', 'ae_4'),
  no_dtw =   c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'area_under_curve',             'fit_likelihood','ae_1', 'ae_2', 'ae_3', 'ae_4'),
  no_velo =  c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'area_under_curve', 'dist_dtw',                  'ae_1', 'ae_2', 'ae_3', 'ae_4'),
  no_ae =    c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'area_under_curve', 'dist_dtw', 'fit_likelihood'                               )
)



#~ methods ----
param_grid <- expand_grid(
  cols = names(subset_list)[[1]],
  method = c("pca_km", "cor_pam", "cor_hc"),
  nb_clust = 4:5,
  nb_pcs = 3:4,
  replicate = 1:10
)

#~ run
results <- pmap_df(param_grid,
                   \(cols, method, nb_clust, nb_pcs, replicate) {
                     
                     run_once_with_bic(cols,
                                       method,
                                       nb_clust,
                                       nb_pcs,
                                       replicate)
                   },
                   .progress = TRUE) |>
  as_tibble()



# qs::qsave(results, file.path(dir_out, "250512_methods_aic.qs"))
# qs::qsave(results, file.path(dir_out, "250512_methods_bic.qs"))
# results <- qs::qread(file.path(dir_out, "250512_methods_bic.qs"))




#~~ plot across methods ----


results

res_long <- results |>
  mutate(nb_clust = as.factor(nb_clust),
         nb_pcs = as.factor(nb_pcs)) |>
  pivot_longer(6:8,
               names_to = "metric")

res_long$cols |> table()

res_long |>
  ggplot() +
  theme_classic() + theme(legend.position = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = method, y = value, fill = interaction(nb_pcs, nb_clust)),
               position = position_dodge(.8))








#~ preprocessing ----
param_grid <- expand_grid(
  cols = names(subset_list)[[1]],
  method = c("pca_km", "som_km"),
  nb_clust = c(4,8),
  nb_comp = c(4,9),
  replicate = 1:5
)

#~ run
results <- pmap_df(param_grid,
                   \(cols, method, nb_clust, nb_comp, replicate) {
                     
                     run_once2(cols,
                               method,
                               nb_clust,
                               nb_comp,
                               replicate)
                   },
                   .progress = TRUE) |>
  as_tibble()



# qs::qsave(results, file.path(dir_out, "250512_methods2.qs"))
# results <- qs::qread(file.path(dir_out, "250512_methods2.qs"))




#~~ plot methods2 ----

results

results |>
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = aic, y = bic))



res_long <- results |>
  mutate(nb_clust = as.factor(nb_clust),
         nb_comp = as.factor(nb_comp)) |>
  pivot_longer(6:9,
               names_to = "metric")


res_long |>
  ggplot() +
  theme_classic() + theme(legend.position = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = method, y = value, fill = interaction(nb_comp, nb_clust)),
               position = position_dodge(.8))








#~ nb clusters ----
param_grid <- expand_grid(
  cols = names(subset_list)[[1]],
  method = c("pca_km"),
  nb_clust = 3:9,
  nb_pcs = 2:5,
  replicate = 1:10
)

#~ run
results <- pmap_df(param_grid,
                   \(cols, method, nb_clust, nb_pcs, replicate) {
                     
                     run_once_with_bic(cols,
                                       method,
                                       nb_clust,
                                       nb_pcs,
                                       replicate)
                   },
                   .progress = TRUE) |>
  as_tibble()

qs::qsave(results, file.path(dir_out, "250512_clust_nb_aic.qs"))
# results <- qs::qread(file.path(dir_out, "250512_clust_nb_aic.qs"))

table(results$cols)
table(results$method)
table(results$nb_clust, results$nb_pcs)




#~~ plot nb clusters ----

results <- qs::qread(file.path(dir_out, "250513_clust_nb.qs"))
results

results |>
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = aic, y = bic))



res_long <- results |>
  mutate(nb_clust = as.factor(nb_clust),
         nb_pcs = as.factor(nb_pcs)) |>
  pivot_longer(6:8,
               names_to = "metric")


res_long |>
  ggplot() +
  theme_classic() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = nb_clust, y = value, fill = nb_pcs))


res_long |>
  # filter(nb_pcs == 3) |>
  # filter(nb_clust == 7) |>
  ggplot() +
  theme_classic() +
  facet_grid(cols = vars(cols),
             rows = vars(metric),
             scales = "free_y") +
  geom_boxplot(aes(x = nb_clust, y = value, fill = nb_pcs))


res_long |>
  filter(cols == "all") |>
  summarize(mean_acc = mean(value),
            sd_acc = sd(value),
            .by = c(nb_clust, nb_pcs, metric)) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_point(aes(x = nb_clust, y = mean_acc, color = nb_pcs),
             position = position_dodge(.8)) +
  geom_errorbar(aes(x = nb_clust,
                    ymin = mean_acc-sd_acc,
                    ymax = mean_acc+sd_acc,
                    color = nb_pcs),
                position = position_dodge(.8)) +
  geom_line(aes(x = nb_clust, y = mean_acc, color = nb_pcs,
                group = nb_pcs),
            position = position_dodge(.8))


res_long |>
  filter(as.numeric(as.character(nb_clust)) > 2) |>
  summarize(mean_acc = mean(value),
            sd_acc = sd(value),
            .by = c(cols, nb_clust, nb_pcs, metric)) |>
  ggplot() +
  theme_classic() +
  facet_grid(cols = vars(cols),
             rows = vars(metric),
             scales = "free_y") +
  geom_point(aes(x = nb_clust, y = mean_acc, color = nb_pcs),
             position = position_dodge(.8)) +
  geom_errorbar(aes(x = nb_clust,
                    ymin = mean_acc-sd_acc,
                    ymax = mean_acc+sd_acc,
                    color = nb_pcs),
                position = position_dodge(.8)) +
  geom_line(aes(x = nb_clust, y = mean_acc, color = nb_pcs,
                group = nb_pcs),
            position = position_dodge(.8))



#~~~ vary one at a time ----

res_long |>
  filter(nb_pcs == 3) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = nb_clust, y = value, fill = nb_pcs))





#~ columns ----
param_grid <- expand_grid(
  cols = names(subset_list),
  method = c("pca_km"),
  nb_clust = 6,
  nb_pcs = 3,
  replicate = 1:10
)

#~ run
results <- pmap_df(param_grid,
                   \(cols, method, nb_clust, nb_pcs, replicate) {
                     
                     run_once_with_bic(cols,
                                       method,
                                       nb_clust,
                                       nb_pcs,
                                       replicate)
                   },
                   .progress = TRUE) |>
  as_tibble()

# qs::qsave(results, file.path(dir_out, "250512_columns.qs"))



#~~ plot columns ----
results


results |>
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(aes(x = aic, y = bic))

res_long <- results |>
  mutate(nb_clust = as.factor(nb_clust),
         nb_pcs = as.factor(nb_pcs)) |>
  pivot_longer(6:8,
               names_to = "metric")


res_long |>
  ggplot() +
  theme_classic() +
  facet_grid(cols = vars(cols),
             rows = vars(metric),
             scales = "free_y") +
  geom_boxplot(aes(x = nb_clust, y = value, fill = nb_pcs))

res_long |>
  ggplot() +
  theme_classic() + theme(legend.position = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = cols, y = value, fill = interaction(nb_pcs, nb_clust)),
               position = position_dodge(.8))



res_long |>
  summarize(mean_acc = mean(value),
            sd_acc = sd(value),
            .by = c(cols, nb_clust, nb_pcs, metric)) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_point(aes(x = cols, y = mean_acc, color = metric),
             position = position_dodge(.8)) +
  geom_errorbar(aes(x = cols,
                    ymin = mean_acc-sd_acc,
                    ymax = mean_acc+sd_acc,
                    color = metric),
                position = position_dodge(.8)) +
  geom_line(aes(x = cols, y = mean_acc, color = metric,
                group = metric),
            position = position_dodge(.8))

res_long |>
  ggplot() +
  theme_classic() + theme(legend.position = "none") +
  facet_wrap(~ metric, scales = "free_y") +
  geom_boxplot(aes(x = cols, y = value, fill = interaction(nb_pcs, nb_clust)),
               position = position_dodge(.8))







