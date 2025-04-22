# test clustering of genes to find pulsatile ones


# Inits ----
library(tidyverse)
suppressPackageStartupMessages( library(Seurat) )
library("ElPiGraph.R")


# source("~/project/larval_devt/R/utils_dist_CORT.R")
source("../larval_devt/R/utils_dist_CORT.R")



#~ functions ----


circ_perm <- function(preds){
  
  peak_loc <- apply(preds, 2, which.max)
  
  preds_permuted <- matrix(nrow = nrow(preds), ncol = ncol(preds))
  t_tot <- nrow(preds)
  t_cent <- floor(t_tot/2)
  
  for(i in seq_len(ncol(preds))){
    
    if(peak_loc[[i]] == t_cent){
      
      preds_permuted[, i] <- preds[,i]
      
    } else if(peak_loc[[i]] == t_tot){
      
      preds_permuted[, i] <- c(preds[ (t_cent+1):t_tot, i],
                               preds[ 1:t_cent ,i])
      
    } else if(peak_loc[[i]] < t_cent){
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      a <- peak_loc[[i]]
      b <- t_cent + a
      
      preds_permuted[, i] <- c(preds[ (b+1):t_tot, i],
                               preds[ 1:b, i])
      
    } else{
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      b <- peak_loc[[i]]
      a <- 1 + b - (t_cent+1)
      
      
      preds_permuted[, i] <- c(preds[ (a+1):b, i],
                               preds[ (b+1):t_tot, i],
                               preds[ 1:a, i])
    }
    
  }
  
  colnames(preds_permuted) <- colnames(preds)
  rownames(preds_permuted) <- rownames(preds)
  
  preds_permuted
}


#~ params ----
params <- list(
  batch_rmed_dir = "intermediates/2502/250330_step1",
  out_dir = "intermediates/2502/250421_step2_boot",
  model = "auto",
  i = 1,
  prop_thres = 0.01,
  cnt_thres = 30,
  nb_subsamples_ptDE = 10
)




n_rep_pt_global <- 50



set.seed(123)


cell_types <- list.files(params$batch_rmed_dir,
                         pattern = "_seu\\.qs$") |>
  stringr::str_remove("_seu\\.qs$")











# Cell type ----
cell_type <- cell_types[[params$i]]

cell_type


message(params$i,"/", length(cell_types), ": ", cell_type)

log_content <- paste0(cell_type, "\n")


message("---- load")

gene_expressions <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type, "_gene_expressions.qs"))
)

seu <- qs::qread(
  file.path(params$batch_rmed_dir,
            paste0(cell_type, "_seu.qs"))
)


# save the plot with phases
gg_phase <- FetchData(seu,
                      vars = c("PC_1", "PC_2", "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked, alpha = cell_rho)) +
  ggtitle(cell_type)

gg_phase

# ggsave(paste0(cell_type,"_phase.png"), gg_phase,
#        path = params$out_dir,
#        width = 7, height = 5, units = "in")





#

#~ Pick genes ----
pos_high_genes <- which(
  (gene_expressions$prop_cells >= params$prop_thres) &
    (gene_expressions$nb_cells >= params$cnt_thres)
)

high_genes <- gene_expressions$gene_name[pos_high_genes]

length(high_genes)



# pred_manual <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")
# high_genes <- pred_manual$gene_name



#~ Normalize ----

embding <- seu |>
  FetchData(vars = c("PC_1", "PC_2")) |>
  as.matrix()


stopifnot(all.equal( rownames(embding), rownames(seu[[]]) ))


#~ pseudotime ----


capture.output(
  CurveEPG <- computeElasticPrincipalCircle(X = embding,
                                            NumNodes = 20,
                                            Do_PCA = FALSE,
                                            nReps = n_rep_pt_global,
                                            drawPCAView = FALSE,
                                            drawAccuracyComplexity = FALSE,
                                            drawEnergy = FALSE,
                                            ProbPoint = .6),
  file = nullfile()
)

stopifnot( length(CurveEPG) == n_rep_pt_global+1L )
targt <- CurveEPG[[ n_rep_pt_global+1L ]]
boots <- CurveEPG[ seq_len(n_rep_pt_global) ]


gg_int_cca_circlepath <- PlotPG(embding,
                                TargetPG = targt, BootPG = boots,
                                Do_PCA = FALSE,VizMode = c("Target","Boot"),
                                DimToPlot = 1:2)[[1]] +
  theme_classic() +
  scale_linewidth_manual(values = c(1,2,3))

gg_int_cca_circlepath

# ggsave(paste0(cell_type,"_circlepath.png"), gg_int_cca_circlepath,
#        path = params$out_dir,
#        width = 6, height = 5, units = "in", dpi = 600)


capture.output(
  subgraph <- GetSubGraph(Net = ConstructGraph(targt),
                          Structure = 'circle')[[ "Circle_1" ]],
  file = nullfile()
)


PartStruct <- PartitionData(X = embding,
                            NodePositions = targt$NodePositions)

ProjStruct <- project_point_onto_graph(X = embding,
                                       NodePositions = targt$NodePositions,
                                       Edges = targt$Edges$Edges,
                                       Partition = PartStruct$Partition)


pseudotime <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(subgraph))[[ "Pt" ]]




#~  GAM ----

message("   -- fit GAM")
mat_cnt <- GetAssayData(seu, assay = "RNA", layer = "count")[high_genes,]



# mods <- lapply(high_genes,
#                \(.gene){
#                  dat <- data.frame(
#                    prob = mat_cnt[.gene,] > 0,
#                    pseudotime = pseudotime/max(pseudotime)
#                  )
#                  
#                  mgcv::gam(prob ~ s(pseudotime, k = 6, bs = 'cc'),
#                            data = dat, family = "binomial")
#                }) |>
#   setNames(high_genes)

mods <- lapply(high_genes,
               \(.gene){
                 dat <- data.frame(
                   prob = mat_cnt[.gene,],
                   pseudotime = pseudotime/max(pseudotime)
                 )
                 
                 mgcv::gam(prob ~ s(pseudotime, k = 6, bs = 'cc'),
                           data = dat, family = mgcv::nb(link = "log"))
               }) |>
  setNames(high_genes)

# qs::qsave(mods,
#           "intermediates/2502/250421_test_autoencoder/250421_AMPHso_models.qs")

message("   -- get distance")

len <- 128
all_preds <- vapply(mods,
                    \(.mod) predict(.mod,
                                    type = "response",
                                    newdata = data.frame(pseudotime = (0:(len-1))/len)),
                    FUN.VALUE = double(len))
all_preds_centered <- circ_perm(all_preds)

dim(all_preds_centered)

# qs::qsave(all_preds_centered,
#           "intermediates/2502/250421_test_autoencoder/250421_AMPHso_all_preds.qs")



# tests ----

library(tidyverse)
library(keras)

pred_manual <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

all_preds_centered <- qs::qread("intermediates/2502/250421_test_autoencoder/250421_AMPHso_all_preds.qs")

stopifnot(all(pred_manual$gene_name %in% colnames(all_preds_centered)))

pheatmap::pheatmap(log1p(all_preds_centered),
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   # scale = "column",
                   # clustering_distance_cols = as.dist(dist_mat),
                   annotation_col = pred_manual |>
                     dplyr::select(manual, gene_name) |>
                     column_to_rownames("gene_name"))


#~ Embedding AE ----

# reticulate::use_virtualenv("r-tensorflow")
# keras::install_keras(method = "virtualenv")


# dataset

len <- nrow(all_preds_centered)

genes <- colnames(all_preds_centered)
x_train <- all_preds_centered[, sample(genes, .8*length(genes))] |> t()
x_test <- all_preds_centered[, setdiff(genes, rownames(x_train))] |> t()

x_train <- x_train |> log1p()
x_test <- x_test |> log1p()

stopifnot(all.equal(
  genes |> sort(),
  union(rownames(x_train), rownames(x_test)) |> sort()
))

dim(x_train);dim(x_test)

matplot(x_train[sample(nrow(x_train), 5),] |> t(), type = "l")




# autoencoder
n_bottleneck <- 8

# Define the encoder
encoder <- keras_model_sequential() |>
  layer_conv_1d(filters = 64, kernel_size = 3, activation = 'relu', input_shape = c(len, 1),
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
  layer_conv_1d_transpose(filters = 64, kernel_size = 3,
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





# Train the model
autoencoder |>
  fit(x_train,
      x_train,
      epochs = 50,
      batch_size = 128,
      validation_data = list(x_test, x_test))





# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = autoencoder$input,
                                        outputs = get_layer(autoencoder, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, t(all_preds_centered))

reconstructed <- predict(autoencoder, t(all_preds_centered))[,,1]

rownames(intermediate_output) <- rownames(reconstructed) <- colnames(all_preds_centered)

# opar <- par(no.readonly = TRUE)

par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)
exple <- sample(rownames(x_train), 5)


matplot(t(x_train)[,exple], type = "l")
matplot(t(reconstructed)[,exple], type = "l")
par(opar)
matplot(t(intermediate_output)[,exple], type = "l")


hc <- hclust(dist(apply(intermediate_output, 1, \(x) x/max(x)) |> t()))

pheatmap::pheatmap(log1p(all_preds_centered),
                   cluster_rows = FALSE,
                   cluster_cols = hc,
                   # scale = "column",
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   annotation_col = pred_manual |>
                     dplyr::select(manual, gene_name) |>
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
        value = "cluster") |> View()


res <- pred_manual |>
  left_join(
    tibble(gene_name = names(km$cluster),
           cluster = km$cluster),
    by = "gene_name"
  ) |>
  arrange(desc(manual), cluster)

table(res$manual, res$cluster)


# examples from each cluster

par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)

matplot(all_preds_centered[,sample(names(km$cluster)[km$cluster == 1], 10)], type = "l")
matplot(all_preds_centered[,sample(names(km$cluster)[km$cluster == 2], 10)], type = "l")
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

























