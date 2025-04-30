# compare gene lists



library(tidyverse)

# Load ----



manual_amphso <- readxl::read_excel("intermediates/2502/250330_step2/manual_AMPHso.xlsx")

manual_several <- readxl::read_excel("intermediates/2502/250330_step2/manual_annotation.xlsx") |>
  dplyr::select(-amplitude, -dev_expl)


# From clustering
dir_step2_clust <- "intermediates/2502/250422_step2_nb_centered/"

cluster_results <- read_csv(file.path(dir_step2_clust, "250424_cluster_results.csv"))

cluster_results |>
  summarize(n_puls = sum(shape == "pulsatile"),
            n_tot = n(),
            .by = "cell_type") |>
  mutate(prop_puls = round( 100 * n_puls / n_tot )) |>
  arrange(prop_puls)



# from bootstrap
dir_step2_boot <- "intermediates/2502/250424_step2_boot_nb/"

boot_ci <- list.files(dir_step2_boot,
                      pattern = "^250424_bootstraps_ci_.*\\.qs$") |>
  map_dfr(\(.ct){
    
    ct_name <- str_match(.ct, "^250424_bootstraps_ci_(.*)\\.qs$")[,2]
    
    qs::qread(file.path(
      dir_step2_boot,
      paste0("250424_bootstraps_ci_",ct_name,".qs")
    )) |>
      add_column(cell_type = ct_name, .before = 1)
  })

table(boot_ci$cell_type)


#~ process bootstrap ----

threshold_ci <- 50

gg <- boot_ci |>
  slice_sample(n = 5000) |>
  dplyr::arrange(upper, t0) |>
  dplyr::mutate(cellgene = forcats::fct_inorder(paste0(cell_type,"|", gene_name))) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  geom_errorbar(aes(x = cellgene,
                    ymin = lower, ymax = upper),
                alpha = .1,
                linewidth = .1,
                width = 0) +
  geom_point(aes(x = cellgene, y = t0)) +
  geom_hline(yintercept = threshold_ci,
             linetype = 'dotted',
             color = 'red3')

# ggsave("250424_boot_selection.png", path = dir_step2_boot, plot = gg,
#        width = 6.5, height = 10, units = "in")



#~ threshold on CI ----

#subset for convenient illustration

# gg_subset <- 
boot_ci |>
  slice_sample(n = 500) |>
  dplyr::arrange(upper) |>
  dplyr::mutate(cellgene = forcats::fct_inorder(paste0(cell_type,"|", gene_name))) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  geom_errorbar(aes(x = cellgene,
                    ymin = lower, ymax = upper),
                alpha = .1,
                linewidth = .1,
                width = 0) +
  geom_point(aes(x = cellgene, y = t0)) +
  geom_hline(yintercept = threshold_ci,
             linetype = 'dotted',
             color = 'red3')

# ggsave("250424_boot_selection_subset2.png", path = dir_step2_boot, plot = gg_subset,
#        width = 6.5, height = 5, units = "in")


ct <- "hypodermis"
boot_ci |>
  filter(cell_type == ct) |>
  left_join(manual_several |> filter(cell_type == ct),
            by = join_by(gene_name)) |>
  filter(!is.na(manual)) |>
  dplyr::arrange(upper) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  scale_color_manual(values = c("mid" = 'grey', "yes" = "green3", "no" = "red3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper, color = manual)) +
  geom_point(aes(x = gene_name, y = t0, color = manual)) +
  coord_flip() +
  geom_hline(yintercept = threshold_ci,
             linetype = 'dotted',
             color = 'red3')


boot_ci |>
  filter(cell_type == "AM_PHso") |>
  left_join(manual_amphso,
            by = join_by(gene_name)) |>
  filter(!is.na(manual)) |>
  dplyr::arrange(upper) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8)) +
  coord_flip() +
  scale_color_manual(values = c("maybe" = 'grey', "yes" = "green3", "no" = "red3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper, color = manual)) +
  geom_point(aes(x = gene_name, y = t0, color = manual)) +
  geom_hline(yintercept = threshold_ci,
             linetype = 'dashed',
             linewidth = 1.5,
             color = 'orange2')


boot_ci |>
  left_join(cluster_results) |>
  filter(!is.na(shape)) |>
  # mutate(guess = case_when(
  #   str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
  #   str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  # )) |>
  # filter(!is.na(guess)) |>
  slice_sample(n = 500) |>
  dplyr::arrange(upper) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  scale_alpha_manual(values = c("nonpulsatile" = .8, "pulsatile" = 1)) +
  scale_size_manual(values = c("nonpulsatile" = .5, "pulsatile" = 1.5)) +
  scale_color_manual(values = c("nonpulsatile" = 'grey', "pulsatile" = "green3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper),
                alpha = .05) +
  geom_point(aes(x = gene_name, y = t0,
                 color = shape,
                 alpha = shape,
                 size =  shape)) +
  geom_hline(yintercept = threshold_ci,
             linetype = 'dotted')




#~ threshold on t0 ----
threshold_t0 <- 25


#subset for convenient illustration

# gg_subset <- 
boot_ci |>
  slice_sample(n = 500) |>
  dplyr::arrange(t0) |>
  dplyr::mutate(cellgene = forcats::fct_inorder(paste0(cell_type,"|", gene_name))) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  geom_errorbar(aes(x = cellgene,
                    ymin = lower, ymax = upper),
                alpha = .1,
                linewidth = .1,
                width = 0) +
  geom_point(aes(x = cellgene, y = t0)) +
  geom_hline(yintercept = threshold_t0,
             linetype = 'dotted',
             color = 'red3')

# ggsave("250424_boot_selection_subset2.png", path = dir_step2_boot, plot = gg_subset,
#        width = 6.5, height = 5, units = "in")


ct <- "hypodermis"
boot_ci |>
  filter(cell_type == ct) |>
  left_join(manual_several |> filter(cell_type == ct),
            by = join_by(gene_name)) |>
  filter(!is.na(manual)) |>
  dplyr::arrange(t0) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  scale_color_manual(values = c("mid" = 'grey', "yes" = "green3", "no" = "red3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper, color = manual)) +
  geom_point(aes(x = gene_name, y = t0, color = manual)) +
  coord_flip() +
  geom_hline(yintercept = threshold_t0,
             linetype = 'dotted',
             color = 'red3')


boot_ci |>
  filter(cell_type == "AM_PHso") |>
  left_join(manual_amphso,
            by = join_by(gene_name)) |>
  filter(!is.na(manual)) |>
  dplyr::arrange(t0) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8)) +
  coord_flip() +
  scale_color_manual(values = c("maybe" = 'grey', "yes" = "green3", "no" = "red3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper, color = manual)) +
  geom_point(aes(x = gene_name, y = t0, color = manual)) +
  geom_hline(yintercept = threshold_t0,
             linetype = 'dashed',
             linewidth = 1.5,
             color = 'orange2')


boot_ci |>
  left_join(cluster_results) |>
  filter(!is.na(shape)) |>
  # mutate(guess = case_when(
  #   str_detect(gene_name, "col\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "cutl\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "grl\\-[0-9]+") ~ "puls",
  #   str_detect(gene_name, "rps\\-[0-9]+") ~ "nonpuls",
  #   str_detect(gene_name, "rpl\\-[0-9]+") ~ "nonpuls"
  # )) |>
  # filter(!is.na(guess)) |>
  slice_sample(n = 500) |>
  dplyr::arrange(t0) |>
  dplyr::mutate(gene_name = forcats::fct_inorder(gene_name)) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  scale_alpha_manual(values = c("nonpulsatile" = .8, "pulsatile" = 1)) +
  scale_size_manual(values = c("nonpulsatile" = .5, "pulsatile" = 1.5)) +
  scale_color_manual(values = c("nonpulsatile" = 'grey', "pulsatile" = "green3")) +
  geom_errorbar(aes(x = gene_name,
                    ymin = lower, ymax = upper),
                alpha = .05) +
  geom_point(aes(x = gene_name, y = t0,
                 color = shape,
                 alpha = shape,
                 size =  shape)) +
  geom_hline(yintercept = threshold_t0,
             linetype = 'dotted')



#~ do threshold ----


table(ci = boot_ci$upper <= threshold_ci,
      t0 = boot_ci$t0 <= threshold_t0)


results_both <- full_join(
  boot_ci |>
    mutate(dtw = if_else(upper <= threshold_ci & t0 <= threshold_t0,
                         "pulsatile",
                         "nonpulsatile")) |>
    select(cell_type, gene_name, dtw)
  ,
  cluster_results,
  by = c("cell_type", "gene_name")
) |>
  mutate(cellgene = paste0(cell_type, "|", gene_name)) |>
  filter(!is.na(shape),
         !is.na(dtw))



list(dtw = results_both$cellgene[ which(results_both$dtw == "pulsatile") ],
     shape = results_both$cellgene[ which(results_both$shape == "pulsatile") ]) |>
  eulerr::euler() |>
  plot(quantities = TRUE)



res_amph <- results_both |> filter(cell_type == "AM_PHso")

list(dtw = res_amph$cellgene[ which(res_amph$dtw == "pulsatile") ],
     shape = res_amph$cellgene[ which(res_amph$shape == "pulsatile") ]) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

res_amph |>
  filter(dtw != shape) |> View()

mat_preds <- qs::qread(file.path(dir_step2_clust,
                                 "AM_PHso_preds.qs"))

plot(mat_preds[,"nhr-31"], type = "l")










