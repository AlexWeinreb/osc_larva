


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
  
  mat_sc1 <- mat_sc[ sample(nrow(mat_sc), prop*nrow(mat_sc)),
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
  
  mat_sc1 <- mat_sc[ sample(nrow(mat_sc), prop*nrow(mat_sc)),
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

