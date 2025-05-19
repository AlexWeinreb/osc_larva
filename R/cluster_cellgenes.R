# Run clustering

message("Starting, ", date())

opar <- par(no.readonly = TRUE)


# Inits ----
library(tidyverse)

# load
dir_step2 <- "intermediates/2502/250519_step2/"

predictors <- list.files(dir_step2,
                         pattern = "_descriptors\\.qs$") |>
  map_dfr(~qs::qread(file.path(dir_step2, .x))) |>
  as_tibble()


message("Selecting cell-genes")
std_cutoff <- .2
hist(predictors$dev_expl); abline(v = std_cutoff, lty = 'dashed', lwd = 2)
table(predictors$dev_expl > std_cutoff)

predictors_sel <- predictors |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(dev_expl > std_cutoff)



message("Scaling predictors")

# # we want a low baseline, it's log (bc counts), and can get quite high, basically anything above 1 is probably bad
# hist(log10(1 + predictors_sel$baseline))
# hist(exp(-log1p(predictors_sel$baseline)))
# plot(log10(1 + predictors_sel$baseline), exp(-log1p(predictors_sel$baseline)))
# 
# # we like high peaks, to avoid flat genes. It's log (bc counts). There are outliers.
# hist(log10(1 + predictors_sel$max_peak))
# hist(DescTools::Winsorize(log10(1 + predictors_sel$max_peak)))
# 
# # we want something close to an ideal peak
# hist(predictors_sel$dist_dtw)
# plot(predictors_sel$dist_dtw, exp(-.05 * predictors_sel$dist_dtw))
# 
# 
# # we don't want asymmetric curves. Looks like above 1 is displaying quite a bit of asymmetry
# hist(predictors_sel$asymmetry)
# plot(predictors_sel$asymmetry, exp(- .2 * predictors_sel$asymmetry), log = "x")



predictors_transformed <- predictors_sel |>
  mutate(low_baseline = exp(-log1p(baseline)),
         peak_amplitude = DescTools::Winsorize(log10(1 + max_peak)),
         shape = exp(-.05 * dist_dtw),
         symmetry = exp(- .2 * asymmetry))


mat_pred_trans <- predictors_transformed |>
  select(cellgene, low_baseline, peak_amplitude, shape, symmetry) |>
  column_to_rownames("cellgene") |>
  as.matrix()


message("Scale")

# par(mfrow = c(2,2), mar = c(3, 2, 2, 1) + 0.1)
# 
# for(i in seq_len(ncol(mat_pred_trans))){
#   hist(mat_pred_trans[,i], main = colnames(mat_pred_trans)[[i]])
# }
# par(opar)


mat_pred <- scale(mat_pred_trans)



# par(mfrow = c(2,2), mar = c(3, 2, 2, 1) + 0.1)
# 
# for(i in seq_len(ncol(mat_pred))){
#   hist(mat_pred[,i], main = colnames(mat_pred)[[i]])
# }
# par(opar)



qs::qsave(mat_pred, file.path(dir_step2, "250519_mat_predictors.qs"))




message("Clustering!")

set.seed(123)
hc <- fastcluster::hclust(dist(mat_pred, method = "euclidean"), method = "ward.D2")

message("Done. Saving...")

qs::qsave(hc, file.path(dir_step2, "250519_hclust.qs"))

message("-----------------")

sessionInfo()



