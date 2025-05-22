# Inits ----

library(tidyverse)




dir_step2 <- "intermediates/2502/250512_step2/"

dir_out <- "intermediates/2502/250512_clustering"

mat_sc <- qs::qread(file.path(dir_step2, "250513_mat_predictors_scaled.qs"))



manual <- readxl::read_excel("data/manual_annotations/250507_manual.xlsx") |>
  select(2:4) |>
  filter(manual %in% c("yes", "no"))


prop <- .3

source("R/compare_kmeans_utils.R")


# Parameters

library(getopt)



if(! interactive()){
  spec <- matrix(c(
    'which', 'f', 1, 'character'
  ), byrow=TRUE, ncol=4)
  
  params <- getopt(spec)
  
} else{
  # Options for interactive
  params <- list(
    which = "nb_clust_pca"
  )
}



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



# run ----
set.seed(123)

if(params$which == "methods"){
  
  
  #~ methods ----
  
  param_grid <- expand_grid(
    cols = names(subset_list)[[1]],
    method = c("pca_km", "cor_pam", "cor_hc"),
    nb_clust = 4:5,
    nb_comp = 3:4,
    replicate = 1:10
  )
  
  #~ run
  
  message("Running! ", date())
  
  results <- pmap_df(param_grid,
                     \(cols, method, nb_clust, nb_comp, replicate) {
                       
                       run_once_algos(cols,
                                      method,
                                      nb_clust,
                                      nb_comp,
                                      replicate)
                     },
                     .progress = TRUE) |>
    as_tibble()
  
  
  
  qs::qsave(results, file.path(dir_out, "250513_methods.qs"))
  
} else if(params$which == "nb_clust_pca"){
  
  #~ nb clusters pca ----
  
  param_grid <- expand_grid(
    cols = names(subset_list)[[1]],
    method = c("pca_km"),
    nb_clust = 3:9,
    nb_comp = 2:6,
    replicate = 1:10
  )
  
  #~ run
  results <- pmap_df(param_grid,
                     \(cols, method, nb_clust, nb_comp, replicate) {
                       
                       run_once_km(cols,
                                   method,
                                   nb_clust,
                                   nb_comp,
                                   replicate)
                     },
                     .progress = TRUE) |>
    as_tibble()
  
  qs::qsave(results, file.path(dir_out, "250513_clust_pca.qs"))
  
  
} else if(params$which == "nb_clust_som"){
  
  #~ nb clusters som ----
  
  param_grid <- expand_grid(
    cols = names(subset_list)[[1]],
    method = c("som_km"),
    nb_clust = 3:9,
    nb_comp = c(3,6,8,10,15,20),
    replicate = 1:10
  )
  
  #~ run
  results <- pmap_df(param_grid,
                     \(cols, method, nb_clust, nb_comp, replicate) {
                       
                       run_once_km(cols,
                                   method,
                                   nb_clust,
                                   nb_comp,
                                   replicate)
                     },
                     .progress = TRUE) |>
    as_tibble()
  
  qs::qsave(results, file.path(dir_out, "250513_clust_som.qs"))
  
  
} else if(params$which == "columns"){
  
  #~ columns ----
  param_grid <- expand_grid(
    cols = names(subset_list),
    method = c("pca_km"),
    nb_clust = 6,
    nb_comp = 3,
    replicate = 1:10
  )
  
  #~ run
  results <- pmap_df(param_grid,
                     \(cols, method, nb_clust, nb_comp, replicate) {
                       
                       run_once_km(cols,
                                   method,
                                   nb_clust,
                                   nb_comp,
                                   replicate)
                     },
                     .progress = TRUE) |>
    as_tibble()
  
  qs::qsave(results, file.path(dir_out, "250513_columns.qs"))
  
} else{
  stop("Which not recognized")
}


message("Done ", date())


message("-------- ")

sessionInfo()



