
# this is a copy of previous version of `genes_by_cell_type_4_analysis.R`,
# useful to decide the thresholds for step 3 (curve shape)


# Inits ----
library(tidyverse)
library(Seurat)

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


dir_step1 <- "intermediates/2502/250319_step1"
dir_step2 <- "intermediates/2502/250325_step2_binom/"





all_raw <- list.files(dir_step2, ".+res_gam.qs$") |>
  map(~ qs::qread( file.path(dir_step2, .x) ) ) |>
  bind_rows()

all_raw |>
  filter(cell_type %in% c("AM_PHso", "intestine")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = cell_type),
             alpha = .2)









unique(all_raw$cell_type)
# individual genes prep ----

ct <- "AM_PHso"
subseu <- qs::qread( file.path(dir_step1, paste0(ct, "_seu.qs")) )



# plot genes ----
appg_genes_strict <- c("abu-11", "pqn-54", "pqn-2","abu-15","abu-1","abu-7","abu-8",
                       "abu-6","abu-14","abu-4","pqn-57", "pqn-71","pqn-13")



appg_genes_strict |> intersect(all_raw$gene_name)
goi <- "abu-14"
goi <- "DH11.5"
goi <- sample(ptDE_proc$gene_name, 1)

gois <- all_raw |>
  filter(cell_type == ct) |>
  # filter(dev_expl  > .05, dev_expl < .25,
  #        amplitude > .45, amplitude < .7) |>
  # filter(area_ratio < 75) |>
  # filter(rmse_sub/rmse > 1) |>
  pull(gene_name) |> fct_inorder() |> sample(10) |> sort() |> as.character()

all_raw |>
  filter(cell_type == ct, gene_name %in% gois)
i=0

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
  geom_point(aes(x = pseudotime,
                 y = count),
             data = tibble(pseudotime = mod$model$pseudotime,
                           count = log1p(mod$y)),
             alpha = .5) +
  geom_line(aes(x = pseudotime,
                y = prediction),
            data = tibble(pseudotime = mod$model$pseudotime,
                          prediction = log1p(mod$fitted.values))) +
  ggtitle(paste0(all_raw$cell_type[[row_nb]], " -- ", all_raw$gene_name[[row_nb]]),
          paste0(" ampl: ", all_raw$amplitude[[row_nb]] |> round(2),
                 " dev: ", all_raw$dev_expl[[row_nb]] |> round(2) )
  )


patchwork::wrap_plots(g1, g2)


#~ check manual selection ----
manual <- readxl::read_excel(file.path(dir_step2, "manual_annotation.xlsx"))



manual |>
  ggplot() +
  theme_classic() +
  scale_color_manual(values = c("grey", "firebrick2", "chartreuse4")) +
  # geom_vline(aes(xintercept = .15), color = 'grey') +
  # geom_hline(aes(yintercept = .65), color = 'grey') +
  geom_point(aes(x = dev_expl, y = amplitude, color = manual, shape = cell_type))

# make model
training <- manual |>
  filter(manual == "no" | manual == "yes") |>
  mutate(manual = case_match(manual,
                             "no"~ "non-peak",
                             "yes" ~ "peak") |>
           factor(levels = c("non-peak", "peak")))

mod <- glmnet::cv.glmnet(x = training |>
                           select(amplitude, dev_expl) |>
                           as.matrix(),
                         y = training |>
                           pull(manual),
                         type.measure = "class",
                         family = "binomial")


manual |>
  ggplot() +
  theme_classic() +
  scale_color_manual(values = c("grey", "firebrick2", "chartreuse4")) +
  geom_point(aes(x = dev_expl, y = amplitude, color = manual, shape = cell_type)) +
  geom_abline(slope = - coef(mod)[3] / coef(mod)[2],
              intercept = - coef(mod)[1] / coef(mod)[2],
              color = 'grey')

all_raw$class <- predict(mod,
                         newx = all_raw |>
                           select(amplitude, dev_expl) |>
                           as.matrix(),
                         type = "class",
                         s = "lambda.1se") |>
  as.character()




all_raw |>
  filter(cell_type %in% c("AM_PHso", "intestine")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dev_expl, y = amplitude, color = cell_type, shape = class),
             alpha = .2)


















