# Run clustering

message("Starting, ", date())

opar <- par(no.readonly = TRUE)


# Inits ----
library(tidyverse)

# load
dir_step2 <- "intermediates/2502/250609_step2/"

dir_clust <- "intermediates/2502/250609_cluster"

predictors <- list.files(dir_step2,
                         pattern = "_descriptors\\.qs$") |>
  map_dfr(~qs::qread(file.path(dir_step2, .x))) |>
  as_tibble() |>
  mutate(cellgene = paste0(cell_type, "|", gene_name))







message("Scaling predictors")

# # we want a low baseline, it's log (bc counts), and can get quite high, basically anything above 1 is probably bad
# hist(log10(1 + predictors$baseline))
# hist(exp(-log1p(predictors$baseline)))
# plot(log10(1 + predictors$baseline), exp(-log1p(predictors$baseline)))
# 
# # we like high peaks, to avoid flat genes. It's log (bc counts). There are outliers.
# hist(log10(1 + predictors$max_peak))
# hist(DescTools::Winsorize(log10(1 + predictors$max_peak)))
# 
# # we want something close to an ideal peak
# hist(predictors$dist_dtw)
# plot(predictors$dist_dtw, exp(-.05 * predictors$dist_dtw))
# 
# 
# # avoid asymmetric curves
# hist(predictors$asymmetry)
# plot(predictors$asymmetry, exp(- .2 * predictors$asymmetry), log = "x")
## > not using: the more highly expressed genes can have quite a bit of asymmetry

# we want a fit where pseudotime is a good predictor of expression
# hist(predictors$dev_expl)



predictors_transformed <- predictors |>
  mutate(low_baseline = exp(-log1p(baseline)),
         peak_amplitude = DescTools::Winsorize(log10(1 + max_peak)),
         shape = exp(-.05 * dist_dtw),
         dev_expl = dev_expl)


mat_pred_trans <- predictors_transformed |>
  select(cellgene, low_baseline, peak_amplitude, shape, dev_expl) |>
  column_to_rownames("cellgene") |>
  as.matrix()



# par(mfrow = c(2,2), mar = c(3, 2, 2, 1) + 0.1)
# 
# for(i in seq_len(ncol(mat_pred_trans))){
#   hist(mat_pred_trans[,i], main = colnames(mat_pred_trans)[[i]])
# }
# par(opar)



message("Scale")


mat_pred <- apply(mat_pred_trans, 2,
      \(col){
        
        bc <- MASS::boxcox(col ~ 1, plotit=FALSE)
        lambda <- bc$x[which.max(bc$y)]
        
        message(lambda)
        
        (col^lambda - 1) / lambda
      })



mat_pred <- scale(mat_pred)


# par(mfrow = c(2,2), mar = c(3, 2, 2, 1) + 0.1)
# 
# for(i in seq_len(ncol(mat_pred))){
#   hist(mat_pred[,i], main = colnames(mat_pred)[[i]])
# }
# par(opar)



qs::qsave(mat_pred, file.path(dir_clust, "250519_mat_predictors.qs"))




message("Clustering Manhattan")

set.seed(123)
hc <- fastcluster::hclust(dist(mat_pred, method = "manhattan"), method = "ward.D2")


message("Done. Saving...")

qs::qsave(hc, file.path(dir_clust, "250519_hclust_manhattan.qs"))



message("Clustering Euclidean")
set.seed(123)
hc <- fastcluster::hclust(dist(mat_pred, method = "euclidean"), method = "ward.D2")


message("Done. Saving...")

qs::qsave(hc, file.path(dir_clust, "250519_hclust_euclidean.qs"))


message("-----------------")

sessionInfo()



