# results from `compare_clustering_methods.R`

# inits ----

library(tidyverse)

dir_step2 <- "intermediates/2502/250422_step2_nb_centered/"

results <- list.files(dir_step2, pattern = "compare_clusterings") |>
  map_dfr(~{
    qs::qread(file.path(dir_step2, .x)) |>
      add_column(replicate = .x)
  }) |>
  as_tibble() |>
  mutate(replicate = str_extract(replicate,
                                 "compare_clusterings([0-9]).qs",
                                 group = 1) |>
           as.factor())




results


results |>
  mutate(centers = as.factor(centers),
         leiden_res = as.factor(leiden_res),
         param = if_else(method == "leiden",
                         leiden_res,
                         centers)) |>
  pivot_longer(cols = 6:8,
               names_to = "metric") |>
  filter(metric == "silhouette") |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(pca),
             cols = vars(columns),
             scales = "free_y") +
  geom_boxplot(aes(x = param, y = value, fill = method))



results |>
  mutate(centers = as.factor(centers),
         leiden_res = as.factor(leiden_res),
         param = if_else(method == "leiden",
                         leiden_res,
                         centers)) |>
  pivot_longer(cols = 6:8,
               names_to = "metric") |>
  filter(metric == "silhouette",
         method != "leiden") |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(pca),
             cols = vars(columns),
             scales = "free_y") +
  geom_boxplot(aes(x = param, y = value, fill = method))


# compare methods across all conditions
results |>
  mutate(centers = as.factor(centers),
         leiden_res = as.factor(leiden_res),
         param = if_else(method == "leiden",
                         leiden_res,
                         centers)) |>
  pivot_longer(cols = 6:8,
               names_to = "metric") |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_wrap(~metric,
             scales = "free_y") +
  geom_boxplot(aes(x = method, y = value))

results |>
  mutate(centers = as.factor(centers),
         leiden_res = as.factor(leiden_res),
         param = if_else(method == "leiden",
                         leiden_res,
                         centers)) |>
  select(columns, pca, param, replicate, method, silhouette) |>
  summarize(max_sil = max(silhouette, na.rm = TRUE),
            .by = c(columns, pca, replicate, method)) |>
  pivot_wider(id_cols = everything(),
              names_from = "method",
              values_from = "max_sil") |>
  summarize(sil_km_hc = kmeans - hclust,
            sil_km_leid = kmeans - leiden,
            .by = c(columns, pca, replicate)) |>
  pivot_longer(cols = starts_with("sil_"),
               names_to = "comparison") |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  geom_boxplot(aes(x = comparison, y = value))
#> kmeans basically always better than Leiden (use silhouette to compare)
#> # using max silhouette accross cluster sizes


# comparing km to hclust at same cl size
results |>
  filter(method != "leiden") |> select(-leiden_res) |>
  summarize(silhouette = mean(silhouette),
            .by = c(columns, pca, replicate, centers, method)) |>
  mutate(centers = as.factor(centers)) |>
  pivot_wider(id_cols = everything(),
              names_from = "method",
              values_from = "silhouette") |>
  # summarize(sil_km_hc = kmeans - hclust,
  #           .by = c(columns, pca, replicate, centers)) |>
  pivot_longer(cols = c(kmeans, hclust),
               names_to = "metric") |>
  filter(as.numeric(centers) > 3) |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_wrap(~ columns) +
  geom_boxplot(aes(x = interaction(pca, centers), y = value,
                   fill = metric, color = metric))
#> kmeans basically always better than Leiden (use silhouette to compare)



# compare cluster nb for kmeans

results |>
  filter(method == "kmeans") |> select(-leiden_res) |>
  summarize(silhouette = mean(silhouette),
            calinski_harabasz = mean(calinski_harabasz),
            davies_bouldin = mean(davies_bouldin),
            .by = c(columns, pca, replicate, centers)) |>
  mutate(centers = as.factor(centers)) |>
  pivot_longer(cols = 5:7,
               names_to = "metric") |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(metric),
             cols = vars(columns),
             scales = "free_y") +
  geom_boxplot(aes(x = centers, y = value,
                   fill = pca, color = pca))

#> always better after PCA


results |>
  filter(method == "kmeans") |> select(-leiden_res) |>
  summarize(silhouette = mean(silhouette),
            calinski_harabasz = mean(calinski_harabasz),
            davies_bouldin = mean(davies_bouldin),
            .by = c(columns, pca, replicate, centers)) |>
  mutate(centers = as.factor(centers)) |>
  pivot_longer(cols = 5:7,
               names_to = "metric") |>
  filter(pca,
         as.numeric(as.character(centers)) > 3) |>
  ggplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(metric),
             cols = vars(columns),
             scales = "free_y") +
  geom_boxplot(aes(x = centers, y = value))

#> no_dev, no_dev_mse, no_ampl_mse with 4-6 centers best result

# all =         c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci', 's1', 's2', 's3', 's4'),
# no_dev_mse =  c('intercept',        'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
# no_dev =      c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw', 'dtw_ci'),
# no_ampl_mse = c('intercept',        'var_s', 'max_s', 'dev_expl',              'auc', 'dtw')



# only plot manually selected subset of cols



# all =         c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci', 's1', 's2', 's3', 's4'),
# coefs =       c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw',           's1', 's2', 's3', 's4'),
# no_coef =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw', 'dtw_ci'),
# dtw_t0 =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
# dtw_ci =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc',        'dtw_ci'),
# no_dtw =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc'),
# no_dev =      c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw', 'dtw_ci'),
# no_mse =      c('intercept',        'var_s', 'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
# no_var_s =    c('intercept', 'mse',          'max_s', 'dev_expl', 'amplitude', 'auc', 'dtw'),
# no_max_s =    c('intercept', 'mse', 'var_s',          'dev_expl', 'amplitude', 'auc', 'dtw'),
# nodev_dtwt0 = c('intercept', 'mse', 'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
# no_ampl =     c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl',              'auc', 'dtw'),
# no_auc =      c('intercept', 'mse', 'var_s', 'max_s', 'dev_expl', 'amplitude',        'dtw'),
# no_maxvar_s = c('intercept', 'mse',                   'dev_expl', 'amplitude', 'auc', 'dtw'),
# no_dev_mse =  c('intercept',        'var_s', 'max_s',             'amplitude', 'auc', 'dtw'),
# no_ampl_mse = c('intercept',        'var_s', 'max_s', 'dev_expl',              'auc', 'dtw')



## all
cols <- c('all','coefs','no_coef','dtw_t0','dtw_ci','no_dtw','no_dev','no_mse',
          'no_var_s','no_max_s','nodev_dtwt0','no_ampl','no_auc','no_maxvar_s',
          'no_dev_mse','no_ampl_mse')

# remove all, coefs, no_coefs
cols <- c('dtw_t0','dtw_ci','no_dtw','no_dev','no_mse',
          'no_var_s','no_max_s','nodev_dtwt0','no_ampl','no_auc','no_maxvar_s',
          'no_dev_mse','no_ampl_mse')















