
# this is a copy of previous version of `genes_by_cell_type_4_analysis.R`,
# useful to decide the thresholds for step 3 (curve shape)


# Inits ----
library(tidyverse)
library(Seurat)
library(mgcv) # if not loaded, `predict()` may call `predict.lm()` instead of gam

library(wbData)


gids <- wb_load_gene_ids(295) |>
  add_row(X = "a",
          gene_id = "nsIs198",
          symbol = "GFP",
          sequence = "GFP",
          status = "Live",
          biotype = "protein_coding_gene",
          name = "GFP"
  )

source("R/utils_heatmap_processing.R")


dir_step1 <- "intermediates/2502/250330_step1/"
dir_step2 <- "intermediates/2502/250330_step2/"

dir_step3 <- "intermediates/2502/250330_step3_genes_by_celltype/"



ct2tissue <- read_csv("data/cell_type2tissue.csv")




# Load ----

all_raw <- list.files(dir_step2, ".+res_gam.qs$") |>
  map(~ qs::qread( file.path(dir_step2, .x) ) ) |>
  bind_rows()

all_raw |>
  filter(cell_type %in% c("AM_PHso", "intestine")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = cell_type),
             alpha = .2)








# Manually annotate some genes ----
# use the manual annotation to chosse threshold

unique(all_raw$cell_type)

ct <- "AM_PHso"
subseu <- qs::qread( file.path(dir_step1, paste0(ct, "_seu.qs")) )



# select a set of genes to annotate
gois <- all_raw |>
  filter(cell_type == ct) |>
  # filter(dev_expl  > .05, dev_expl < .25,
  #        amplitude > .45, amplitude < .7) |>
  pull(gene_name) |> fct_inorder() |> sample(15) |> sort() |> as.character()

all_raw |>
  filter(cell_type == ct, gene_name %in% gois)
i=0

# look at them
i <- i+1
goi <- gois[[i]]

g1 <- Seurat::FeaturePlot(subseu,
                          features = goi,
                          reduction = "pca",
                          pt.size = 2, #min.cutoff = 0,max.cutoff = 1,
                          alpha = .5) +
  ggtitle(goi, s2i(goi, gids))

row_nb <- which(all_raw$gene_name == goi & all_raw$cell_type == ct)
mod <- all_raw$gam_fit[[ row_nb ]]
all_raw[row_nb,]

g2 <- ggplot() +
  theme_classic() +
  ylab(expression(P* "[" * N[UMI] > 0*"]" )) +
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod$model$pseudotime,
                           count = mod$y),
             alpha = .5) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod$model$pseudotime,
                          prediction = mod$fitted.values)) +
  ggtitle(paste0(all_raw$cell_type[[row_nb]], " -- ", all_raw$gene_name[[row_nb]]),
          paste0(" ampl: ", all_raw$amplitude[[row_nb]] |> round(2),
                 " dev: ", all_raw$dev_expl[[row_nb]] |> round(2),
                 " area: ", all_raw$area_under_curve[[row_nb]] |> round() )
  )


patchwork::wrap_plots(g1, g2)



# # illustrate area under curve
# 
# ggplot() +
#   theme_classic() +
#   ylab(expression(P* "[" * N[UMI] > 0*"]" )) +
#   geom_ribbon(aes(x = pseudotime,
#                   ymin = 0,
#                   ymax = prediction),
#               data = tibble(pseudotime = mod$model$pseudotime,
#                             prediction = mod$fitted.values),
#               fill = "lightyellow",
#               alpha = .6) +
#   geom_point(aes(x = pseudotime,
#                  y = count),
#              data = tibble(pseudotime = mod$model$pseudotime,
#                            count = mod$y),
#              alpha = .5) +
#   geom_line(aes(x = pseudotime,
#                 y = prediction),
#             data = tibble(pseudotime = mod$model$pseudotime,
#                           prediction = mod$fitted.values)) +
#   ggtitle(paste0(all_raw$cell_type[[row_nb]], " -- ", all_raw$gene_name[[row_nb]]),
#           paste0(" ampl: ", all_raw$amplitude[[row_nb]] |> round(2),
#                  " dev: ", all_raw$dev_expl[[row_nb]] |> round(2),
#                  " area: ", all_raw$area_under_curve[[row_nb]] |> round() )
#   )



#~ check manual selection results ----
manual <- readxl::read_excel(file.path(dir_step2, "manual_annotation.xlsx")) |>
  select(-amplitude, -dev_expl)

manual <- left_join(manual, all_raw,
                     by = c("cell_type", "gene_name"))


manual |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = 'none') +
  scale_color_manual(values = c("grey", "firebrick2", "chartreuse4")) +
  # geom_vline(aes(xintercept = .15), color = 'grey') +
  # geom_hline(aes(yintercept = .65), color = 'grey') +
  geom_point(aes(x = dev_expl, y = amplitude, color = manual, size = area_under_curve),
             alpha = .6)


# use logistic regression to pick threshold that matches manual annotation
training <- manual |>
  filter(manual == "no" | manual == "yes") |>
  mutate(manual = case_match(manual,
                             "no"~ "non-peak",
                             "yes" ~ "peak") |>
           factor(levels = c("non-peak", "peak")))

mod <- glmnet::cv.glmnet(x = training |>
                           select(amplitude, dev_expl, area_under_curve) |>
                           as.matrix(),
                         y = training |>
                           pull(manual),
                         type.measure = "class",
                         family = "binomial")

# mod <- glm(manual ~ amplitude + dev_expl + area, data = training, family = "binomial")



stopifnot(all( rownames(coef(mod)) == c("(Intercept)", "amplitude", "dev_expl", "area_under_curve") ))



all_raw$peaky <- predict(mod,
                         newx = all_raw |>
                           select(amplitude, dev_expl, area_under_curve) |>
                           as.matrix(),
                         type = "class",
                         s = "lambda.1se") |>
  as.factor()

manual <- left_join(manual |> select(cell_type, gene_name, manual),
                    all_raw,
                    by = c("cell_type", "gene_name"))

manual |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "firebrick2", "chartreuse4")) +
  geom_point(aes(x = dev_expl, y = amplitude, color = manual, size = area_under_curve, shape = peaky),
              alpha = .6) +
  geom_abline(slope = - coef(mod)[3] / coef(mod)[2],
              intercept = - coef(mod)[1] / coef(mod)[2],
              color = 'grey')


# explore annotation
all_raw |>
  filter(cell_type == "AM_PHso") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = peaky),
             alpha = .2)

all_raw |>
  filter(cell_type == "intestine") |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = peaky),
             alpha = .2)

all_raw |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_bar(aes(x = cell_type, fill = peaky))

all_raw |>
  summarize(nb_peaky = sum(peaky == "peak"),
            nb_genes = n(),
            prop_peaky = nb_peaky / nb_genes,
            .by = cell_type) |>
  left_join(ct2tissue, by = "cell_type") |>
  ggplot() +
  theme_classic() +
  labs(x = "Number of genes tested",
       y = "Proportion of genes selected") +
  aes(x = nb_genes, y = prop_peaky, label = cell_type, color = tissue) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = FALSE)




# heatmap ----

# save a heatmap for each cell type

filtered_data <- all_raw |>
  group_by(cell_type) |>
  nest() |>
  ungroup() |>
  mutate(nb_peaky = map_int(data,
                            ~{
                              sum(.x[["peaky"]] == "peak")
                            })) |>
  filter(nb_peaky > 5)

heatmaps_list <- set_names(filtered_data$data,
                           filtered_data$cell_type) |>
  map(
    \(.dat){
      
      all_preds <- vapply(.dat[["gam_fit"]][.dat$peaky == "peak"],
                          \(.mod) predict(.mod,
                                          type = "response",
                                          newdata = data.frame(pseudotime = (0:199)/200)),
                          FUN.VALUE = double(200))
      
      peak_location <- apply(all_preds, 2, which.max)
      
      all_preds[, order(peak_location)]
      
    },
    .progress = TRUE
  )

# sapply(heatmaps_list, dim)



#~ uniformity ----


filtered_data$uniformity_index <- map_dbl(heatmaps_list,
                                          \(hm) {
                                            diffrange <- hm |>
                                              rowMeans() |>
                                              range() |>
                                              diff()
                                            
                                            1 - diffrange
                                          })



hist(filtered_data$uniformity_index, breaks = 50)






#~ distance ----


filtered_data$similarity_diag <- map_dbl(heatmaps_list,
                                         \(hm){
                                           
                                           n_t <- nrow(hm)
                                           n_g <- ncol(hm)
                                           
                                           sig_cent <- make_ref_sig(n_t)
                                           
                                           null_mat <- sapply(seq_len(n_g),
                                                              \(i){
                                                                circ_perm(sig_cent, floor(n_t * i/n_g))
                                                              })
                                           
                                           
                                           dist <- norm(hm - null_mat, "1")
                                           similarity <- dist/n_t
                                           1 - similarity
                                         })


hist(filtered_data$similarity_diag, breaks = 30)


ggplot(filtered_data) +
  theme_classic() +
  geom_text(aes(x = uniformity_index, y = similarity_diag, label = cell_type))



# dists with pvals

filtered_data$p_val <- map_dbl(heatmaps_list,
                               \(hm){
                                 
                                 n_t <- nrow(hm)
                                 n_g <- ncol(hm)
                                 
                                 sig_cent <- make_ref_sig(n_t)
                                 
                                 null_mat <- sapply(seq_len(n_g),
                                                    \(i){
                                                      circ_perm(sig_cent, floor(n_t * i/n_g))
                                                    })
                                 
                                 
                                 
                                 similarity_perms <- c(
                                   {
                                     dist <- norm(hm - null_mat, "1")
                                     similarity <- dist/n_t
                                     1 - similarity
                                   },
                                   replicate(10000,{
                                     hm_perm <- hm[,sample(ncol(hm))]
                                     dist <- norm(hm_perm - null_mat, "1")
                                     similarity <- dist/n_t
                                     1 - similarity
                                   })
                                 )
                                 
                                 pval <- mean(similarity_perms >= similarity_perms[[1]])
                                 pval
                               },
                               .progress = TRUE)


filtered_data$p_adj <- filtered_data$p_val |>
  p.adjust(method = "holm")
hist(filtered_data$p_val, breaks = 30)


ggplot(filtered_data) +
  theme_classic() +
  aes(x = uniformity_index, y = similarity_diag, label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()


ggplot(filtered_data) +
  theme_classic() +
  aes(x = similarity_diag, y = -log10(p_adj), label = cell_type, shape = p_adj < .05, color = p_adj < .05) +
  geom_point() +
  ggrepel::geom_text_repel()


# # only label every n-th column
# ct <- "BWM"
# hmp_sparsified <- heatmaps_list[[ct]]
# hmp_sparsified <- null_mat
# 
# spar_index <- 7
# colnames_to_sparsify <- setdiff(seq_len(ncol(hmp_sparsified)),
#                                 spar_index * seq_len( ncol(hmp_sparsified) / spar_index ) )
# 
# colnames(hmp_sparsified)[colnames_to_sparsify] <- ""
# head(colnames(hmp_sparsified), 20)
# 
# 
# pheatmap::pheatmap(hmp_sparsified,
#                    cluster_rows = FALSE,
#                    cluster_cols = FALSE,
#                    show_rownames = FALSE,
#                    fontsize = 7,
#                    main = ct)





# Save ----

all_raw |>
  select(- gam_fit) |>
  qs::qsave(file.path(dir_step3, "all_genes.qs"))


filtered_data |>
  select(-data) |>
  qs::qsave(file.path(dir_step3, "cell_types.qs"))

iwalk(heatmaps_list,
      ~ squash::savemat(t(.x)[, nrow(.x):1],
                        filename = file.path(dir_step3,
                                             paste0(.y, "_heatmap.png"))))


iwalk(heatmaps_list,
      ~ qs::qsave(.x,
                  file.path(dir_step3,
                            paste0(.y, "_heatmap.qs")))
)
  
  
  



















