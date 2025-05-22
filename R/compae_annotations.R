library(Seurat)
library(tidyverse)
library(wbData)

gids_230 <- wb_load_gene_ids(230)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")



osc_raw <- readxl::read_excel("../10x_grl18/data/oscillating/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID),
         gene_name = i2s(gene_id, gids) )




# load datasets ----


# older
oseu <- qs::qread("../10x_grl18/intermediates/240620_reanalysis_panglial_2021/241022_ds_with_angles.qs")
Idents(oseu) <- "tissue"


DimPlot(oseu,
        group.by = "cell_type",
        label = TRUE) +
  NoLegend()

# 2502 reclustering
dir_out <- "intermediates/2502/250328_assembled"
nseu <- qs::qread( file.path(dir_out, "250329_seu_all_herma.qs"))


DimPlot(nseu,
        label = TRUE) +
  NoLegend()


# Compare ----
ometa <- oseu[[]] |>
  rownames_to_column("cell_bc") |>
  select(cell_bc, orig.ident, tissue, cell_type, cell_rho, cell_phase) |>
  filter(! orig.ident %in% c("231019GRL", "240111")) |>
  mutate(
    old_orig.ident = orig.ident,
    orig.ident = case_match(
    orig.ident,
    "210901_glia_WS280_pla_2020-07-30_batch1_CHB3840b" ~ "200730_batch1_CHB3840b",
    "210901_glia_WS280_pla_2020-10-13_batch2_CHB3840b_CEG_fqs" ~ "201013_batch2_CHB3840b_CEG_fqs",
    "210901_glia_WS280_pla_2021-04-13_batch3_CHB3840b" ~ "210413_batch3_CHB3840b",
    "210901_glia_WS280_pla_2021-04-27_batch5_CHB3840b" ~ "210427_batch5_CHB3840b",
    "220210_glia_WS280_pla_2022-02-10_batch6_OH17400/" ~ "220210_OH17400",
    .default = orig.ident
    )) |>
  mutate(cell_bc2 = cell_bc |> str_remove("^s[0-9]_") |> str_remove("_[0-9]+$"),
         orig_bc = paste0(orig.ident, "|", cell_bc2) ) |>
  as_tibble()


nmeta <- nseu[[]] |>
  rownames_to_column("cell_bc") |>
  select(cell_bc, orig.ident, tissue, cell_type, cell_rho, cell_phase) |>
  mutate(cell_bc2 = cell_bc |> str_remove("^s[0-9]+_") |> str_remove("_[0-9]+$"),
         orig_bc = paste0(orig.ident, "|", cell_bc2) ) |>
  as_tibble()

all.equal(
  sort(unique(ometa$orig.ident)),
  sort(unique(nmeta$orig.ident))
  )


list(old = ometa$orig_bc,
     new = nmeta$orig_bc) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


setdiff(ometa$orig_bc, nmeta$orig_bc) |> sample(10)
setdiff(nmeta$orig_bc, ometa$orig_bc) |> sample(10)


#~ compare annotations ----
comb <- full_join(
  ometa |> select(orig_bc, tissue_old = tissue, cell_type_old = cell_type),
  nmeta |> select(orig_bc, tissue_new = tissue, cell_type_new = cell_type),
  by = "orig_bc"
)

comb |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = tissue_old, y = tissue_new),
             alpha = .1) +
  geom_text(
    aes(x = tissue_old, y = tissue_new,
            label = nb_cells),
    data = comb |>
      summarize(nb_cells = n(),
                .by = c(tissue_old, tissue_new)),
    color = 'red3'
  )




comb |>
  filter(tissue_new == "glia" | tissue_old == "glia") |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_jitter(aes(x = cell_type_old, y = cell_type_new),
              alpha = .1) +
  geom_text(
    aes(x = cell_type_old, y = cell_type_new,
        label = nb_cells),
    data = comb |>
      filter(tissue_new == "glia" | tissue_old == "glia") |>
      summarize(nb_cells = n(),
                .by = c(cell_type_old, cell_type_new)),
    color = 'red3'
  )



comb |>
  filter(cell_type_new == "ILso" | cell_type_new == "PHsh" |
           cell_type_old == "ILso" | cell_type_old == "PHsh") |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_jitter(aes(x = cell_type_old, y = cell_type_new),
              alpha = .1) +
  geom_text(
    aes(x = cell_type_old, y = cell_type_new,
        label = nb_cells),
    data = comb |>
      filter(cell_type_new == "ILso" | cell_type_new == "PHsh" |
               cell_type_old == "ILso" | cell_type_old == "PHsh") |>
      summarize(nb_cells = n(),
                .by = c(cell_type_old, cell_type_new)),
    color = 'red3'
  )

#~ mean rho ----
comb <- full_join(
  ometa |> select(orig_bc, tissue_old = tissue, cell_type_old = cell_type, cell_rho_old = cell_rho),
  nmeta |> select(orig_bc, tissue_new = tissue, cell_type_new = cell_type, cell_rho_new = cell_rho),
  by = "orig_bc"
)

comb |>
  filter(cell_type_old == "ILso",
         cell_type_new == "ILso") |>
  pivot_longer(cols = starts_with("cell_rho_"),
               names_prefix = "cell_rho_",
               names_to = "annot",
               values_to = "cell_rho") |>
  ggplot() +
  theme_minimal() +
  geom_density(aes(x = cell_rho, fill = annot),
               alpha = .5)


#~ mean phase ----
comb <- full_join(
  ometa |> select(orig_bc, tissue_old = tissue, cell_type_old = cell_type, cell_phase_old = cell_phase),
  nmeta |> select(orig_bc, tissue_new = tissue, cell_type_new = cell_type, cell_phase_new = cell_phase),
  by = "orig_bc"
)

comb |>
  filter(cell_type_old == "ILso",
         cell_type_new == "ILso") |>
  # pivot_longer(cols = starts_with("cell_rho_"),
  #              names_prefix = "cell_rho_",
  #              names_to = "annot",
  #              values_to = "cell_rho") |>
  ggplot() +
  theme_minimal() +
  geom_point(aes(x = cell_phase_old, y = cell_phase_new),
               alpha = .5)



#~ gene expression ----

comb <- full_join(
  ometa,
  nmeta,
  by = "orig_bc"
)

old_bcs_ilso <- comb |>
  filter(cell_type.x == "ILso",
         cell_type.y == "ILso") |>
  pull(cell_bc.x)

new_bcs_ilso <- comb |>
  filter(cell_type.x == "ILso",
         cell_type.y == "ILso") |>
  pull(cell_bc.y)

genes_common <- intersect(rownames(oseu),
                          s2i(rownames(nseu), gids))

# all(genes_common %in% rownames(oseu))
# all(i2s(genes_common, gids, warn_missing = TRUE) %in% rownames(nseu))


old_sub <- GetAssayData(oseu,
                        assay = "RNA",
                        layer = "counts")[genes_common, old_bcs_ilso]


new_sub <- GetAssayData(nseu,
                        assay = "RNA",
                        layer = "counts")[i2s(genes_common, gids), new_bcs_ilso]

dim(old_sub); dim(new_sub)

all.equal(i2s(rownames(old_sub), gids, warn_missing = TRUE),
          rownames(new_sub))
rownames(old_sub) <- i2s(rownames(old_sub), gids, warn_missing = TRUE)


cbind(old = rowSums(old_sub),
      new = rowSums(new_sub)) |>
  as.data.frame() |>
  rownames_to_column("gene_name") |>
  mutate(is_bulk_osc = gene_name %in% osc_raw$gene_name[osc_raw$OscAmplitude > 1.5]) |>
  # slice_sample(n = 10000) |>
  ggplot() +
  theme_classic() +
  scale_x_continuous(transform = "log1p") +
  scale_y_continuous(transform = "log1p") +
  geom_point(aes(x = old, y = new, color = is_bulk_osc),
             alpha = .2)

patchwork::wrap_plots(
  cbind(old = rowSums(old_sub),
        new = rowSums(new_sub)) |>
    as.data.frame() |>
    rownames_to_column("gene_name") |>
    mutate(is_bulk_osc = gene_name %in% osc_raw$gene_name[osc_raw$OscAmplitude > 1.5]) |>
    # slice_sample(n = 10000) |>
    ggplot() +
    theme_classic() +
    # scale_x_continuous(transform = "log1p") +
    scale_y_continuous(transform = "log1p") +
    geom_hline(aes(yintercept = 1)) +
    geom_boxplot(aes(x = is_bulk_osc, y = new / old, fill = is_bulk_osc)) +
    ggtitle("all genes")
  ,
  cbind(old = rowSums(old_sub),
        new = rowSums(new_sub)) |>
    as.data.frame() |>
    rownames_to_column("gene_name") |>
    filter(old > 1e3 | new > 1e3) |>
    mutate(is_bulk_osc = gene_name %in% osc_raw$gene_name[osc_raw$OscAmplitude > 1.5]) |>
    # slice_sample(n = 10000) |>
    ggplot() +
    theme_classic() +
    # scale_x_continuous(transform = "log1p") +
    scale_y_continuous(transform = "log1p") +
    geom_hline(aes(yintercept = 1)) +
    geom_boxplot(aes(x = is_bulk_osc, y = new / old, fill = is_bulk_osc)) +
    ggtitle("only genes with > 1000 reads")
)



#~~ single cell ----
osc_table <- osc_raw |>
  filter(gene_name %in% rownames(new_sub),
         Class == "Osc") |>
  select(gene_name, gene_id,
         osc_amplitude = OscAmplitude,
         peak_phase_deg = PeakPhase) |>
  group_by(gene_name) |>
  slice_sample(n = 1) |>
  ungroup()

comb <- full_join(
  ometa |> select(orig_bc, tissue_old = tissue, cell_type_old = cell_type, cell_rho_old = cell_rho, cell_bc_old = cell_bc),
  nmeta |> select(orig_bc, tissue_new = tissue, cell_type_new = cell_type, cell_rho_new = cell_rho, cell_bc_new = cell_bc),
  by = "orig_bc"
)

candidates <- comb |>
  filter(cell_type_old == "ILso",
         cell_type_new == "ILso") |>
  filter(cell_rho_old > .7)

dat_old <- enframe(old_sub[,candidates$cell_bc_old[[1]]],
                      name = "gene_name",
                      value = "expression") |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(!is.na(osc_amplitude))

dat_old |>
  ggplot() +
  theme_bw() +
  coord_polar() +
  scale_x_continuous(limits = c(0,360), n.breaks = 15) +
  ylab("expression") +
  geom_segment(aes(x = peak_phase_deg,
                   xend = peak_phase_deg,
                   y = 0,
                   yend = expression),
               linewidth = .25) +
  geom_segment(aes(x = mean_angle %% (360),
                   xend = mean_angle %% (360),
                   y = 0,
                   yend = mean_rho),
               data = dat_old |>
                 summarize(mean_angle = circhelp::weighted_circ_mean(peak_phase_deg*pi/180, expression)*180/pi,
                           mean_rho = max(expression)*circhelp::weighted_circ_rho(peak_phase_deg*pi/180, expression)),
               linewidth = 1,
               color = 'red3')



dat_new <- enframe(new_sub[,candidates$cell_bc_new[[1]]],
                   name = "gene_name",
                   value = "expression") |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(!is.na(osc_amplitude))

dat_new |>
  ggplot() +
  theme_bw() +
  coord_polar() +
  scale_x_continuous(limits = c(0,360), n.breaks = 15) +
  ylab("expression") +
  geom_segment(aes(x = peak_phase_deg,
                   xend = peak_phase_deg,
                   y = 0,
                   yend = expression),
               linewidth = .25) +
  geom_segment(aes(x = mean_angle %% (360),
                   xend = mean_angle %% (360),
                   y = 0,
                   yend = mean_rho),
               data = dat_new |>
                 summarize(mean_angle = circhelp::weighted_circ_mean(peak_phase_deg*pi/180, expression)*180/pi,
                           mean_rho = max(expression)*circhelp::weighted_circ_rho(peak_phase_deg*pi/180, expression)),
               linewidth = 1,
               color = 'red3')






# Reannotate old rho ----

osc_table <- osc_table |> filter(osc_amplitude > 1.5)
mat <- LayerData(oseu, layer = "data", features = osc_table$gene_id)

stopifnot(all(osc_table$gene_id %in% rownames(oseu)))
mat <- mat[osc_table$gene_id,]

# genes_max <- sparseMatrixStats::rowMaxs(mat)
# genes_max[genes_max == 0] <- 1
# 
# mat <- mat / genes_max

stopifnot(identical(
  rownames(mat),
  osc_table$gene_id
))
stopifnot(identical(
  colnames(mat),
  rownames(oseu[[]])
))


# oseu$cell_phase <- angle_from_mat(mat, osc_table$peak_phase_deg)
oseu$cell_rho2 <- rho_from_mat(mat, osc_table$peak_phase_deg)



# new rho ----


osc_table <- osc_table |> filter(osc_amplitude > 1.5)
mat <- LayerData(nseu, layer = "data", features = osc_table$gene_name)


# genes_max <- sparseMatrixStats::rowMaxs(mat)
# genes_max[genes_max == 0] <- 1
# 
# mat <- mat / genes_max

stopifnot(identical(
  rownames(mat),
  osc_table$gene_name
))
stopifnot(identical(
  colnames(mat),
  rownames(nseu[[]])
))


# oseu$cell_phase <- angle_from_mat(mat, osc_table$peak_phase_deg)
nseu$cell_rho2 <- rho_from_mat(mat, osc_table$peak_phase_deg)



#~~ from old scripts ----

wb_ids_osc <- intersect(rownames(oseu),
                        osc_raw$WB_ID)

osc_noperm <- osc_raw  |>
  filter(WB_ID %in% wb_ids_osc) |>
  filter(Class == "Osc",
         OscAmplitude > 1.5)

mat <- GetAssayData(oseu)[osc_noperm$WB_ID[!is.na(osc_noperm$PeakPhase)] |> intersect(rownames(oseu)),]
genes_max <- sparseMatrixStats::rowMaxs(mat)
genes_max[genes_max == 0] <- 1

mat <- mat / genes_max


stopifnot(identical(rownames(mat),
                    osc_noperm$WB_ID))
stopifnot(identical(colnames(mat),
                    rownames(oseu[[]])))


oseu$cell_rho2 <- rho_from_mat(mat, osc_noperm$PeakPhase)







# Redo cell rho ----

ometa <- oseu[[]] |>
  rownames_to_column("cell_bc") |>
  select(cell_bc, orig.ident, tissue, cell_type, cell_rho, cell_phase, cell_rho2) |>
  filter(! orig.ident %in% c("231019GRL", "240111")) |>
  mutate(
    old_orig.ident = orig.ident,
    orig.ident = case_match(
      orig.ident,
      "210901_glia_WS280_pla_2020-07-30_batch1_CHB3840b" ~ "200730_batch1_CHB3840b",
      "210901_glia_WS280_pla_2020-10-13_batch2_CHB3840b_CEG_fqs" ~ "201013_batch2_CHB3840b_CEG_fqs",
      "210901_glia_WS280_pla_2021-04-13_batch3_CHB3840b" ~ "210413_batch3_CHB3840b",
      "210901_glia_WS280_pla_2021-04-27_batch5_CHB3840b" ~ "210427_batch5_CHB3840b",
      "220210_glia_WS280_pla_2022-02-10_batch6_OH17400/" ~ "220210_OH17400",
      .default = orig.ident
    )) |>
  mutate(cell_bc2 = cell_bc |> str_remove("^s[0-9]_") |> str_remove("_[0-9]+$"),
         orig_bc = paste0(orig.ident, "|", cell_bc2) ) |>
  as_tibble()


nmeta <- nseu[[]] |>
  rownames_to_column("cell_bc") |>
  select(cell_bc, orig.ident, tissue, cell_type, cell_rho, cell_phase, cell_rho2) |>
  mutate(cell_bc2 = cell_bc |> str_remove("^s[0-9]+_") |> str_remove("_[0-9]+$"),
         orig_bc = paste0(orig.ident, "|", cell_bc2) ) |>
  as_tibble()


comb <- full_join(
  ometa |> select(orig_bc, tissue_old = tissue, cell_type_old = cell_type, cell_rho_old = cell_rho, cell_rho_old2 = cell_rho2),
  nmeta |> select(orig_bc, tissue_new = tissue, cell_type_new = cell_type, cell_rho_new = cell_rho, cell_rho_new2 = cell_rho2),
  by = "orig_bc"
)

comb |>
  filter(cell_type_old == "ILso",
         cell_type_new == "ILso") |>
  pivot_longer(cols = starts_with("cell_rho_"),
               names_prefix = "cell_rho_",
               names_to = "annot",
               values_to = "cell_rho") |>
  ggplot() +
  theme_minimal() +
  geom_density(aes(x = cell_rho, fill = annot),
               alpha = .5)



