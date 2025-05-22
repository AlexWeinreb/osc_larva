
# Inits ----
library(tidyverse)


library(getopt)


if(! interactive()){
  spec <- matrix(c(
    'i', 'i', 1, 'integer'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    i = 1
  )
}


set.seed(params$i)


dir_step2 <- "intermediates/2502/250502_step2/"

mat_sc <- qs::qread(file.path(dir_step2, "mat_predictors_scaled.qs"))



# Testing function ----
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




evaluate_clustering <- function(data, method = c("kmeans", "hclust", "leiden"),
                                nb_clust, nb_pcs, verbose = FALSE) {
  
  method <- match.arg(method)
  
  
  if(verbose){
    message("---------------------- ")
    message("   method: ",method)
    message("   nb_clust: ",nb_clust)
    message("   nb_pcs: ",nb_pcs)
  }
  
  
  
  pca_res <- prcomp(data)
  data <- pca_res$x[,seq_len(nb_pcs)]
  
  
  # Cluster assignment
  if (method == "kmeans") {
    km <- stats::kmeans(data, centers = nb_clust, nstart = 25)
    clusters <- km$cluster
    
  } else if (method == "hclust") {
    hc <- fastcluster::hclust(stats::dist(data))
    clusters <- stats::cutree(hc, k = nb_clust)
    
  } else if (method == "leiden") {
    gr <- get_knn_graph(data, k = knn_k)
    part <- leidenbase::leiden_find_partition(
      gr,
      partition_type = "RBConfigurationVertexPartition",
      resolution_parameter = nb_clust
    )
    clusters <- purrr::chuck(part, "membership")
    clusters <- remove_singletons(clusters, SNN = igraph::as_adjacency_matrix(gr))
  }
  
  # Internal metrics
  int_metrics <- clusterCrit::intCriteria(
    as.matrix(data),
    as.integer(clusters),
    c("Calinski_Harabasz", "Davies_Bouldin")
  )
  
  # Silhouette (on sampled subset only)
  
  
  idx <- sample(seq_len(nrow(data)), 4000)
  sil_data <- data[idx, , drop = FALSE]
  sil_clusters <- clusters[idx]
  sil <- tryCatch(cluster::silhouette(sil_clusters, stats::dist(sil_data)),
                  error = \(e) matrix(NA, ncol = 3))
  sil_score <- mean(sil[, 3])
  
  
  data.frame(
    method = method,
    nb_pcs = nb_pcs,
    nb_clust = nb_clust,
    m_silhouette = sil_score,
    m_calinski_harabasz = int_metrics$calinski_harabasz,
    m_davies_bouldin = int_metrics$davies_bouldin
  )
}






#~ run ----

param_grid <- expand_grid(
  method = c("kmeans", "hclust"),
  nb_clust = 2:10,
  nb_pcs = 2:(ncol(mat_sc) - 1),
  replicate = params$i
)

results <- pmap_df(param_grid,
                   \(method, nb_clust, nb_pcs, replicate) {
                     
                     mat_sc1 <- mat_sc[sample(.5*nrow(mat_sc)),]
                     
                     evaluate_clustering(mat_sc1,
                                         method,
                                         nb_clust,
                                         nb_pcs,
                                         verbose = TRUE)
                   })

qs::qsave(results,
          file.path(dir_step2,
                    paste0("compare_clusterings_", params$i, ".qs")))

sessionInfo()

















# ____Old version_____ ----

#~~ Define conditions ----
# 
# all_cols <- colnames(mat_sc)
# 
# subset_list <- list(
#   all =         c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 's1', 's2', 's3', 's4'),
#   coefs =       c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw',           's1', 's2', 's3', 's4'),
#   no_coef =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
#   dtw_t0 =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
#   dtw_ci =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc'),
#   no_dtw =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc'),
#   no_dev =      c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
#   no_mse =      c('intercept',        'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_var_s =    c('intercept', 'mse',          'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_max_s =    c('intercept', 'mse', 'var_s',          'dev_expl', 'amplitude', 'auc', 'dtw'),
#   nodev_dtwt0 = c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
#   no_ampl =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl',              'auc', 'dtw'),
#   no_auc =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude',        'dtw'),
#   no_maxvar_s = c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw'),
#   no_dev_mse =  c('intercept',        'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
#   no_ampl_mse = c('intercept',        'var_s', 'max_s', 'dev_expl',              'auc', 'dtw')
# )
# 
# method <- c("kmeans", "hclust", "leiden")
# pca = c(TRUE, FALSE)
# 
# param_grid <- expand_grid(
#   cols_nm = names(subset_list),
#   method,
#   pca,
#   replicate = params$i,
#   centers = 2:10,
#   leiden_res = seq(.1, 1.5, .2)
# )
# 
# results <- pmap_df(param_grid,
#                    \(cols_nm, method, pca, replicate, centers, leiden_res) {
#                      
#                      mat_sc1 <- mat_sc[sample(.3*nrow(mat_sc)),]
#                      
#                      evaluate_clustering(mat_sc1,
#                                          cols_nm,
#                                          method = method,
#                                          pca = pca,
#                                          centers = centers,
#                                          leiden_res = leiden_res,
#                                          verbose = TRUE)
#                    })
# 
# qs::qsave(results,
#           file.path(dir_step2,
#                     paste0("compare_clusterings", params$i, ".qs")))
# 
# sessionInfo()




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








