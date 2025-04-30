
# Inits ----
library(tidyverse)

dir_step2 <- "intermediates/2502/250422_step2_nb_centered/"

mat <- qs::qread(file.path(dir_step2, "mat_predictors.qs"))
mat <- mat[!is.na(mat[,"dtw"]),]

mat_sc <- apply(mat, 2, DescTools::Winsorize) |> scale()





# Testing function ----
evaluate_clustering <- function(data, cols_nm, method = c("kmeans", "hclust", "leiden"),
                                pca = FALSE,
                                centers = 6, sil_sample_size = 1000,
                                knn_k = 30, leiden_res = 0.5, verbose = FALSE) {
  
  
  cols <- subset_list[[cols_nm]]
  
  subset_data <- data[, cols, drop = FALSE]
  method <- match.arg(method)
  
  
  if(verbose){
    message("---------------------- ")
    message("   pca: ",pca)
    message("   method: ",method)
    message("   centers: ",centers)
    message("   leiden_res: ",leiden_res)
  }
  
  
  if(pca){
    pca_res <- prcomp(subset_data)
    subset_data <- pca_res$x[,1:3]
  }
  
  # Cluster assignment
  if (method == "kmeans") {
    km <- stats::kmeans(subset_data, centers = centers, nstart = 25)
    clusters <- km$cluster
    
  } else if (method == "hclust") {
    hc <- fastcluster::hclust(stats::dist(subset_data))
    clusters <- stats::cutree(hc, k = centers)
    
  } else if (method == "leiden") {
    gr <- get_knn_graph(subset_data, k = knn_k)
    part <- leidenbase::leiden_find_partition(
      gr,
      partition_type = "RBConfigurationVertexPartition",
      resolution_parameter = leiden_res
    )
    clusters <- purrr::chuck(part, "membership")
    clusters <- remove_singletons(clusters, SNN = igraph::as_adjacency_matrix(gr))
  }
  
  # Internal metrics
  int_metrics <- clusterCrit::intCriteria(
    as.matrix(subset_data),
    as.integer(clusters),
    c("Calinski_Harabasz", "Davies_Bouldin")
  )
  
  # Silhouette (on sampled subset only)
  sil_score <- NA
  if (nrow(subset_data) > sil_sample_size) {
    idx <- sample(seq_len(nrow(subset_data)), sil_sample_size)
    sil_data <- subset_data[idx, , drop = FALSE]
    sil_clusters <- clusters[idx]
    sil <- tryCatch(cluster::silhouette(sil_clusters, stats::dist(sil_data)),
                    error = \(e) matrix(NA, ncol = 3))
    sil_score <- mean(sil[, 3])
  }
  
  data.frame(
    columns = cols_nm,
    method = method,
    pca = pca,
    centers = centers,
    leiden_res = leiden_res,
    silhouette = sil_score,
    calinski_harabasz = int_metrics$calinski_harabasz,
    davies_bouldin = int_metrics$davies_bouldin
  )
}


# Define conditions ----

all_cols <- colnames(mat_sc)

subset_list <- list(
  all =         c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci', 's1', 's2', 's3', 's4'),
  coefs =       c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw',           's1', 's2', 's3', 's4'),
  no_coef =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
  dtw_t0 =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
  dtw_ci =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc',        'dtw_ci'),
  no_dev =      c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw', 'dtw_ci'),
  no_mse =      c('intercept',        'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
  no_var_s =    c('intercept', 'mse',          'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
  no_max_s =    c('intercept', 'mse', 'var_s',          'dev_expl', 'amplitude', 'auc', 'dtw'),
  nodev_dtwt0 = c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
  no_ampl =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl',              'auc', 'dtw'),
  no_auc =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude',        'dtw'),
  no_maxvar_s = c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw'),
  no_dev_mse =  c('intercept',        'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
  no_ampl_mse = c('intercept',        'var_s', 'max_s', 'dev_expl',              'auc', 'dtw')
)

method <- c("kmeans", "hclust", "leiden")
pca = c(TRUE, FALSE)

param_grid <- expand_grid(
  cols_nm = names(subset_list),
  method,
  pca,
  replicate = 1:5,
  centers = 2:10,
  leiden_res = seq(.1, 1.5, .2)
)

results <- pmap_df(param_grid,
                   \(cols_nm, method, pca, replicate, centers, leiden_res) {
                     
                     mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
                     
                     evaluate_clustering(mat_sc1,
                                         cols_nm,
                                         method = method,
                                         pca = pca,
                                         centers = centers,
                                         leiden_res = leiden_res,
                                         verbose = TRUE)
                   })

qs::qsave(results, file.path(dir_step2, "compare_clusterings.qs"))

sessionInfo()




# #___Individual___ ----
# 
# # separate passes to vary fewer parameters at a time
# 
# 
# 
# 
# #~ 1: algos ----
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   all = c('intercept', 's1', 's2', 's3', 's4', 'mse', 'var_s', 'max_s',
#           'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
#   no_coef = c('intercept', 'mse', 'var_s', 'max_s',
#           'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci')
#   # dtw_t0 = c('intercept', 'mse', 'var_s', 'max_s',
#   #            'dev_expl', 'amplitude', 'auc', 'dtw'),
#   # dtw_ci = c('intercept', 'mse', 'var_s', 'max_s',
#   #         'dev_expl', 'amplitude', 'auc', 'dtw_ci'),
#   # nodev = c('intercept', 'mse', 'var_s', 'max_s',
#   #         'amplitude', 'auc', 'dtw', 'dtw_ci')
# )
# 
# method <- c("kmeans", "hclust", "leiden")
# pca = c(TRUE, FALSE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5
# )
# 
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca)
#                    },
#                    .progress = TRUE)
# 
# results |>
#   ggplot() +
#   theme_classic() +
#   facet_wrap(~ method) +
#   geom_boxplot(aes(x = columns, fill = pca, y = silhouette))
# 
# patchwork::wrap_plots(
#   results |>
#     ggplot() +
#     theme_classic() +
#     theme(legend.position = "none") +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, fill = pca, y = calinski_harabasz))
#   ,
#   results |>
#     ggplot() +
#     theme_classic() +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, fill = pca, y = davies_bouldin))
# )
# 
# 
# 
# results |>
#   mutate(do_pca = if_else(pca, "PC", "")) |>
#   filter(
#     # method == "kmeans",
#     columns != "all"
#     # startsWith(columns, "dtw")
#          # columns %in% c("all", "no_coef")
#          ) |>
#   ggplot() +
#   theme_classic() +
#   geom_point(aes(x = silhouette, y = davies_bouldin)) +
#   ggrepel::geom_text_repel(aes(x = silhouette,
#                                y = davies_bouldin,
#                                label = paste(columns, method, do_pca)))
# 
# 
# 
# 
# 
# 
# 
# #~ 2: vary columns for kmeans ----
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   all = c('intercept', 's1', 's2', 's3', 's4', 'mse', 'var_s', 'max_s',
#           'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
#   no_coef = c('intercept', 'mse', 'var_s', 'max_s',
#               'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
#   dtw_t0 = c('intercept', 'mse', 'var_s', 'max_s',
#              'dev_expl', 'amplitude', 'auc', 'dtw'),
#   dtw_ci = c('intercept', 'mse', 'var_s', 'max_s',
#           'dev_expl', 'amplitude', 'auc', 'dtw_ci'),
#   nodev = c('intercept', 'mse', 'var_s', 'max_s',
#           'amplitude', 'auc', 'dtw', 'dtw_ci')
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5
# )
# 
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca)
#                    },
#                    .progress = TRUE)
# 
# 
# 
# results |>
#   ggplot() +
#   theme_classic() +
#   facet_wrap(~ method) +
#   geom_boxplot(aes(x = columns, fill = pca, y = silhouette))
# 
# patchwork::wrap_plots(
#   results |>
#     ggplot() +
#     theme_classic() +
#     theme(legend.position = "none") +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, fill = pca, y = calinski_harabasz))
#   ,
#   results |>
#     ggplot() +
#     theme_classic() +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, fill = pca, y = davies_bouldin))
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #~ 3: vary more columns for kmeans ----
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_coef = c('intercept', 'mse', 'var_s', 'max_s',
#               'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
#   dtw_t0 = c('intercept', 'mse', 'var_s', 'max_s',
#              'dev_expl', 'amplitude', 'auc', 'dtw'),
#   nodev = c('intercept', 'mse', 'var_s', 'max_s',
#             'amplitude', 'auc', 'dtw', 'dtw_ci'),
#   nodev_dtwt0 = c('intercept', 'mse', 'var_s', 'max_s',
#                   'amplitude', 'auc', 'dtw'),
#   no_max_s = c('intercept', 'mse', 'var_s',
#                'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_var_s = c('intercept', 'mse', 'max_s',
#                'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_maxvar_s = c('intercept', 'mse',
#                'dev_expl', 'amplitude', 'auc', 'dtw'),
#   coefs = c('intercept', 's1', 's2', 's3', 's4', 'mse',
#           'dev_expl', 'amplitude', 'auc', 'dtw')
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca)
#                    },
#                    .progress = TRUE)
# 
# 
# 
# results |>
#   pivot_longer(cols = 4:6,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = columns, y = value))
# 
# patchwork::wrap_plots(
#   results |>
#     ggplot() +
#     theme_classic() +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, y = silhouette))
#   ,
#   results |>
#     ggplot() +
#     theme_classic() +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, y = calinski_harabasz))
#   ,
#   results |>
#     ggplot() +
#     theme_classic() +
#     facet_wrap(~ method) +
#     geom_boxplot(aes(x = columns, y = davies_bouldin))
# )
# 
# 
# 
# 
# #~ 4: vary ampl columns for kmeans ----
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_coef = c('intercept', 'mse', 'var_s', 'max_s',
#               'dev_expl', 'amplitude', 'auc', 'dtw'),
#   nodev = c('intercept', 'mse', 'var_s', 'max_s',
#             'amplitude', 'auc', 'dtw'),
#   no_ampl = c('intercept', 'mse', 'var_s', 'max_s',
#               'dev_expl', 'auc', 'dtw'),
#   no_auc = c('intercept', 'mse', 'var_s', 'max_s',
#              'dev_expl', 'amplitude', 'dtw'),
#   no_mse = c('intercept', 'var_s', 'max_s',
#              'dev_expl', 'amplitude', 'auc', 'dtw')
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca)
#                    },
#                    .progress = TRUE)
# 
# 
# 
# results |>
#   pivot_longer(cols = 4:6,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = columns, y = value))
# 
# 
# 
# #~ 5: vary no mse ----
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_mse =  c('intercept', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_dev =  c('intercept', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
#   no_ampl = c('intercept', 'var_s', 'max_s', 'dev_expl',              'auc', 'dtw')
#   
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca)
#                    },
#                    .progress = TRUE)
# 
# 
# 
# results |>
#   pivot_longer(cols = 4:6,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = columns, y = value))
# 
# 
# 
# 
# 
# 
# #~ 6: nb clusters kmeans ----
# 
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_mse =  c('intercept', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw')
#   
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5,
#   centers = 3:10
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate, centers) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca,
#                                          centers = centers)
#                    },
#                    .progress = TRUE)
# 
# results |> head()
# 
# 
# results |>
#   mutate(centers = as.factor(centers)) |>
#   pivot_longer(cols = 5:7,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = centers, y = value))
# 
# 
# 
# 
# 
# #~ 7&8: higher subsample ----
# 
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_mse =  c('intercept', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw')
#   
# )
# 
# method <- c("kmeans")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5,
#   centers = 3:9
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate, centers) {
#                      
#                      mat_sc1 <- mat_sc[sample(.6*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca,
#                                          centers = centers)
#                    },
#                    .progress = TRUE)
# 
# results |> head()
# 
# 
# results |>
#   mutate(centers = as.factor(centers)) |>
#   pivot_longer(cols = 5:7,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = centers, y = value))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #~ 9: vary clusters for hclust ----
# 
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_mse =  c('intercept', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw')
#   
# )
# 
# method <- c("hclust")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5,
#   centers = 3:9
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate, centers) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca,
#                                          centers = centers)
#                    },
#                    .progress = TRUE)
# 
# results |> head()
# 
# 
# results |>
#   mutate(centers = as.factor(centers)) |>
#   pivot_longer(cols = 5:7,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = centers, y = value))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #~ 10: vary clusters for Leiden ----
# 
# 
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   no_mse =  c('intercept', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw')
#   
# )
# 
# method <- c("leiden")
# pca = c(TRUE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = 1:5,
#   leiden_res = seq(.1, 1.5, .2)
# )
# 
# set.seed(123)
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate, leiden_res) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca,
#                                          leiden_res = leiden_res)
#                    },
#                    .progress = TRUE)
# 
# results |> head()
# 
# 
# results |>
#   mutate(leiden_res = as.factor(leiden_res)) |>
#   pivot_longer(cols = 6:8,
#                names_to = "metric") |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
#   facet_wrap(~ metric, scales = "free_y") +
#   geom_boxplot(aes(x = leiden_res, y = value))








