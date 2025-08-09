# Inits ----

opar <- par(no.readonly = TRUE)

library(tidyverse)
library(ggrastr)
library(Seurat)
library(wbData)

source("R/utils.R")
source("R/mean_phase_rho.R")

# SCTransform requires 1.3 GB for this data
# options(future.globals.maxSize = 1.5 * 1024^3)



gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")




osc_raw <- readxl::read_excel("data/msb209498-sup-0003-datasetev1.xlsx",
                              sheet = "Dataset EV1 WBidToGeneNames_Osc",
                              na = "NA") |>
  mutate(gene_id = wb_clean_gene_names(WB_ID),
         gene_name = i2s(gene_id, gids) )







dir_out <- "intermediates/2502/250605_assembled"
# dir_out_individual_cts <- file.path(dir_out, "250330_cell_types")

dir_third_processed <- "intermediates/2502/250529_third_processed/"





# Just herma from third pass ----

#~ load ----


files_list <- list.files(dir_third_processed) |>
  enframe(value = "filename",
          name = NULL) |>
  mutate(filename2 = filename |>
           str_replace("seu_1","g1") |>
           str_replace("seu_2","g2")) |>
  filter(str_detect(filename2, "230[35][0-9]{2}_(g1|g2)")) |>
  separate_wider_regex(filename2,
                       patterns = c(
                         "^230[35][0-9]{2}_",
                         group = "(?:g1|g2)",
                         "_",
                         cell_type = "[a-zA-Z0-9_]+",
                         "\\.qs$"
                       ))


seu_list <- files_list |>
  mutate(stage = case_match(group,
                            "g1" ~ "L2",
                            "g2" ~ "L4")) |>
  pmap(\(filename, group, cell_type, stage){
    seu <- qs::qread(file.path( dir_third_processed, filename ))
    
    DefaultAssay(seu) <- "RNA"
    seu[["SCT"]] <- NULL
    
    seu <- JoinLayers(seu)
    
    seu$stage <- stage
    seu$cell_type <- cell_type
    Idents(seu) <- cell_type
    
    seu
    
  },
  .progress = TRUE)


seu <- merge_fast(seu_list, cell_keys = TRUE)

rm(seu_list)

# qs::qsave(seu, file.path(dir_out,
#                          "250606_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250606_seu_all_herma.qs"))


#~ tissue ----

ct2tissue <- readxl::read_excel("data/cell_annotations.xlsx",
                                sheet = 1L) |>
  select(cell_type, tissue) |>
  column_to_rownames("cell_type")

seu$tissue <- ct2tissue[ seu$cell_type , "tissue" ]

stopifnot( !any(is.na(seu$tissue)) )


#~ norm and dim reduc ----
seu <- SCTransform(seu)

seu <- RunPCA(seu, npcs = 200, verbose = FALSE)

npca <- 40

ElbowPlot(seu, ndims = 200) +
  geom_vline(aes(xintercept = npca))





DimPlot(seu,
        group.by = "stage",
        reduction = "pca",
        pt.size = 2,
        alpha = .05) +
  NoLegend()


seu <- RunUMAP(seu, dims = 1:npca)


DimPlot(seu,
        group.by = "orig.ident",
        reduction = "umap",
        pt.size = 2,
        alpha = .05) +
  NoLegend()


DimPlot(seu,
        group.by = "stage",
        reduction = "umap",
        pt.size = 2,
        alpha = .05)

DimPlot(seu,
        group.by = "cell_type",
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .05) +
  NoLegend()

DimPlot(seu,
        group.by = "tissue",
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .05) +
  NoLegend()




# Export UMAPs ----


# seu2 <- seu
# seu2$tissue[seu2$cell_type == "ILso"] <- "ILso"
# DimPlot(seu2,
#         group.by = "tissue",
#         reduction = "umap",
#         label = FALSE,
#         pt.size = 2,
#         alpha = .1) +
#   NoLegend()
# 
# 
# FetchData(seu2, vars = c("tissue", "umap_1", "umap_2")) |>
#   ggplot() +
#   theme_classic() +
#   geom_point(aes(x = umap_1, y = umap_2, color = tissue),
#              alpha = .1, size = 2,
#              show.legend = FALSE)
# 
# ggsave("UMAP_tissue.png", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)
# ggsave("UMAP_tissue.pdf", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)
# rm(seu2)



DimPlot(seu,
        group.by = "cell_type",
        reduction = "umap",
        label = FALSE,
        pt.size = 2,
        alpha = .1) +
  NoLegend()


FetchData(seu, vars = c("cell_type", "umap_1", "umap_2")) |>
  ggplot() +
  theme_classic() +
  hues::scale_color_iwanthue() +
  geom_point(aes(x = umap_1, y = umap_2, color = cell_type),
             alpha = .1, size = 2,
             show.legend = FALSE)

# ggsave("UMAP_cell_type.png", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)
# ggsave("UMAP_cell_type.pdf", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)

DimPlot(seu,
        group.by = "cell_type",
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .1) +
  NoLegend()



DimPlot(seu,
        group.by = "stage",
        reduction = "umap",
        label = FALSE,
        pt.size = 2,
        alpha = .1)


FetchData(seu, vars = c("stage", "umap_1", "umap_2")) |>
  ggplot() +
  theme_classic() +
  scale_color_brewer(type = "qual") +
  geom_point(aes(x = umap_1, y = umap_2, color = stage),
             alpha = .1, size = 2,
             show.legend = FALSE)

# ggsave("UMAP_stage.png", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)

FetchData(seu, vars = c("stage", "umap_1", "umap_2")) |>
  slice_sample(n = 10) |>
  ggplot() +
  theme_classic() +
  scale_color_brewer(type = "qual") +
  geom_point(aes(x = umap_1, y = umap_2, color = stage),
             alpha = 1, size = 2)

# ggsave("UMAP_stage_legend.pdf", path = "presentations/figures/250610_umap/",
#        width = 60, height = 60, units = "mm",
#        scale = 2)


# colorspace::darken("#beaed4", amount = .2)




# Oscillations ----
# seu <- qs::qread( file.path(dir_out, "250606_seu_all_herma.qs"))

osc_table <- osc_raw |>
  filter(gene_name %in% rownames(seu),
         Class == "Osc",
         OscAmplitude > 1.5) |>
  select(gene_name, gene_id,
         osc_amplitude = OscAmplitude,
         peak_phase_deg = PeakPhase) |>
  group_by(gene_name) |>
  slice_sample(n = 1) |>
  ungroup() |>
  arrange(factor(gene_name, levels = rownames(seu)))

# there are duplicated gene names, removed with slice_sample(); to look at them:
# duplicated_genes <- osc_table$gene_name[duplicated(osc_table$gene_name)]
# map(duplicated_genes,
#     ~ filter(osc_table, gene_name == .x))




#~ average phase and rho ----


mat <- LayerData(seu, assay = "SCT", layer = "data", features = osc_table$gene_name)


# normalize by max across all cells
# (doesn't seem to make a difference in practice, but in principle,
# a gene could "look" high in a cell it's not expressed in if it's much higher
# in other cells)
genes_max <- sparseMatrixStats::rowMaxs(mat)
genes_max[genes_max == 0] <- 1

mat <- mat / genes_max


stopifnot(identical(
  rownames(mat),
  osc_table$gene_name
))
stopifnot(identical(
  colnames(mat),
  rownames(seu[[]])
))


seu$cell_phase <- angle_from_mat(mat, osc_table$peak_phase_deg)
seu$cell_rho <- rho_from_mat(mat, osc_table$peak_phase_deg)



# qs::qsave(seu, file.path(dir_out,
#                          "250606_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250606_seu_all_herma.qs"))






# FetchData(seu,
#           vars = c("umap_1", "umap_2", "cell_phase", "cell_rho")) |>
#   ggplot() +
#   theme_classic() +
#   theme(legend.position = "none") +
#   scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
#                         limits = c(0, 360)) +
#   geom_point(aes(x = umap_1, y = umap_2,
#                  color = cell_phase,
#                  alpha = cell_rho),
#              alpha = .1)




# permutation test on rho

set.seed(123)
pvals_by_cell <- tibble(perm = 0:10000) |>
  rowwise() |>
  mutate(rho_data = map(perm,
                        \(.perm){
                          if(perm == 0){
                            rho_from_mat(mat, osc_table$peak_phase_deg) |>
                              enframe(name = "cell_bc",
                                      value = "rho")
                          } else{
                            rho_from_mat(mat, sample(osc_table$peak_phase_deg)) |>
                              enframe(name = "cell_bc",
                                      value = "rho")
                          }
                        },
                        .progress = TRUE)) |>
  unnest(rho_data) |>
  nest(.by = cell_bc) |>
  mutate(p_val = map_dbl(data,
                         \(dat){
                           mean( dat$rho >= dat$rho[[ 1 ]] )
                         })) |>
  select(-data) |>
  mutate(FDR = p.adjust(p_val, method = "BH")) |>
  column_to_rownames("cell_bc")

# qs::qsave(pvals_by_cell, file.path(dir_out, "250606_permutations_10000.qs"))

pvals_by_cell <- qs::qread(file.path(dir_out, "250606_permutations_10000.qs"))


# hist(pvals_by_cell$p_val, breaks = 30, main = NULL, xlab = "Distribution of p-values")
# hist(pvals_by_cell$FDR, breaks = 50)

table(pvals_by_cell$FDR < .05)
#> FALSE  TRUE 
#>  9606 14451


seu$length_FDR <- pvals_by_cell[rownames(FetchData(seu, vars = "ident")), "FDR"]
seu$length_signif <- (pvals_by_cell[rownames(FetchData(seu, vars = "ident")), "FDR"] < .05)
seu$cell_phase_masked <- if_else(seu$length_signif, seu$cell_phase, NA_real_)




# qs::qsave(seu, file.path(dir_out,
#                          "250606_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250606_seu_all_herma.qs"))



# plot phases

FetchData(seu,
          vars = c("umap_1", "umap_2",
                   "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  aes(x = umap_1, y = umap_2) +
  geom_point(aes(color = cell_phase_masked),
             shape = 16,
             size = 2,
             alpha = .2,
             show.legend = FALSE)

# ggsave("umap_phase.png", path = "presentations/figures/250610_umap/",
#        width = 80, height = 80, units = "mm",
#        scale = 1.5)
# ggsave("umap_phase.pdf", path = "presentations/figures/250610_umap/",
#        width = 80, height = 80, units = "mm",
#        scale = 1.5)







#~ Local coherence index ----

k_neighbors <- 20L


seu <- FindNeighbors(seu,
                     return.neighbor = TRUE,
                     reduction = "pca",
                     dims = 1:npca,
                     k.param = k_neighbors + 1L)


# qs::qsave(seu, file.path(dir_out, "250606_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250606_seu_all_herma.qs"))






cell_phases <- FetchData(seu, vars = c("cell_phase", "cell_rho", "tissue", "cell_type"))


## Normalized

# dotprod_by_cell <- mean_dotprod_norm(cell_phases,
#                                      seu.nn = seu@neighbors$SCT.nn,
#                                      k = k_neighbors)
# dotprod_by_cell <- qs::qread(file.path(dir_out, "250606_dotprod_by_cell.qs"))


## Unnormalized
dotprod_by_cell <- cbind(cell_phases,
                         coherence = mean_dotprod(cell_phases,
                                                  seu.nn = seu@neighbors$SCT.nn,
                                                  k = k_neighbors)) |>
  as_tibble()

# qs::qsave(dotprod_by_cell, file.path(dir_out, "250610_dotprod_by_cell.qs"))
dotprod_by_cell <- qs::qread(file.path(dir_out, "250610_dotprod_by_cell.qs"))



# Plot by cell type and cluster

cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 30) |>
  pull(cell_type)

dotprod_agg_by_ct <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  summarize(mean_coherence = mean(coherence),
            .by = c(tissue, cell_type)) |>
  mutate(tissue = factor(tissue,
                         levels = c("skin", "glia","pharynx","muscle","neuron","reproductive","other"))) |>
  arrange(tissue, desc(mean_coherence)) |>
  mutate(cell_type = fct_inorder(cell_type))

dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  mutate(cell_type = factor(cell_type, levels = levels(dotprod_agg_by_ct$cell_type))) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence") + xlab(NULL) +
  # geom_hline(aes(yintercept = -.2)) + geom_hline(aes(yintercept = .9)) +
  # coord_cartesian(ylim = c(0,.9)) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = coherence, color = tissue),
                               alpha = .2,
                               shape = 16) +
  geom_point(aes(x = cell_type, y = mean_coherence),
             data = dotprod_agg_by_ct)





#~~ perm test ----

# the neighbors etc do not change between permutations: compute once and reuse
precomputed <- local({
  
  
  cell_phases = FetchData(seu,
                          vars = c("cell_phase", "cell_rho",
                                   "tissue", "cell_type"))
  
  
  cell_neighbors = tibble(
    cell = rownames(cell_phases),
    neighbors = map(cell,
                    \(cell) TopNeighbors(seu@neighbors$SCT.nn,
                                         cell = cell,
                                         n = (k_neighbors + 1L) ) |>
                      setdiff(cell),
                    .progress = TRUE)
  ) |>
    unnest(neighbors)
  
  # use a subsample
  cell_same_type = tibble(
    cell = rownames(cell_phases),
    cell_type = cell_phases$cell_type,
    nb_cells_in_type = map_int(cell_type,
                               \(.ct) sum(cell_phases$cell_type == .ct)),
    nb_cells_use = pmin(nb_cells_in_type - 1L, k_neighbors),
    neighbors = pmap(list(.ct = cell_type, .n = nb_cells_use, .cell = cell),
                     \(.ct, .n, .cell) rownames(cell_phases)[cell_phases$cell_type == .ct] |> setdiff(.cell) |> sample(.n),
                     .progress = TRUE)
  ) |>
    select(cell, neighbors) |>
    unnest(neighbors) |>
    filter(cell != neighbors)
  
  
  neigh_cell_indices = match(cell_neighbors$cell, rownames(cell_phases))
  neigh_neigh_indices = match(cell_neighbors$neighbors, rownames(cell_phases))
  
  
  type_cell_indices = match(cell_same_type$cell, rownames(cell_phases))
  type_neigh_indices = match(cell_same_type$neighbors, rownames(cell_phases))
  
  
  list(cell_phases = cell_phases,
       cell_neighbors = cell_neighbors, cell_same_type = cell_same_type,
       neigh_cell_indices = neigh_cell_indices, neigh_neigh_indices = neigh_neigh_indices,
       type_cell_indices = type_cell_indices, type_neigh_indices = type_neigh_indices
  )
})

# qs::qsave(precomputed, file.path(dir_out, "250606_precomputed_perm.qs"))
# precomputed <- qs::qread(file.path(dir_out, "250606_precomputed_perm.qs"))

xx <- run_permutation_test_by_celltype_and_phase(0, precomputed, mat, osc_table$peak_phase_deg)

all.equal(
  dotprod_by_cell |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type)) |> arrange(tissue, cell_type),
  xx |> filter(permutation == 0) |> select(-permutation) |> arrange(tissue, cell_type)
)

plot(
  dotprod_by_cell |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type)) |> arrange(tissue, cell_type) |> pull(mean_coherence),
  xx |> filter(permutation == 0) |> select(-permutation) |> arrange(tissue, cell_type) |> pull(mean_coherence)
); abline(a=0, b=1)


# takes ~1h
set.seed(123)
mean_dotprod_by_celltype_res_perm <- map_dfr(
  0:10000,
  ~ run_permutation_test_by_celltype_and_phase(.x, precomputed, mat, osc_table$peak_phase_deg),
  .progress = TRUE
)
# qs::qsave(mean_dotprod_by_celltype_res_perm,
#           file.path(dir_out, "250606_coherence_perm10000.qs"))
# mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_out, "250606_coherence_perm10000.qs"))

all.equal(
  mean_dotprod_by_celltype_res_perm |> filter(permutation == 0) |> select(-permutation),
  dotprods_normalized(precomputed, mat, osc_table$peak_phase_deg) |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type))
)



library(furrr)
plan(multicore, workers = 4)
precomputed2 <- precomputed[c("cell_phases", "cell_neighbors",
                              "neigh_cell_indices", "neigh_neigh_indices")]
set.seed(123)
mean_dotprod_by_celltype_res_perm <- future_map_dfr(
  0:10000,
  ~ run_permutation_test_unnorm_by_celltype_and_phase(.x, precomputed2, mat, osc_table$peak_phase_deg),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

# qs::qsave(mean_dotprod_by_celltype_res_perm,
#           file.path(dir_out, "250606_coherence_unnorm_perm10000.qs"))
# mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_out, "250606_coherence_unnorm_perm10000.qs"))






# mean_dotprod_by_celltype_res_perm |>
#   arrange(tissue) |> mutate(cell_type = fct_inorder(cell_type)) |>
#   mutate(`permutated` = !(permutation == 0)) |>
#   ggplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c('red3','grey')) +
#   scale_alpha_manual(values = c(1,.01)) +
#   ylab("Mean (in cell type) of mean (with neighbors) dot product") + xlab(NULL) +
#   geom_point(aes(x = cell_type, y = mean_coherence, color = permutated, alpha = permutated),
#              size = 2)


#~~~ res ----
p_vals <- mean_dotprod_by_celltype_res_perm |>
  group_by(tissue, cell_type) |>
  nest() |>
  summarize(p_val = map_dbl(data,
                            \(dat){
                              mean(dat$mean_coherence >= dat$mean_coherence[[1]])
                            }),
            .groups = 'drop') |>
  mutate(p_adj = p.adjust(p_val, method = "holm"))


hist(p_vals$p_val, breaks = 30)
hist(p_vals$p_adj, breaks = 30, xlim = c(0,1))



# remove spaces in names
dotprod_by_cell <- dotprod_by_cell |>
  mutate(cell_type = str_replace_all(cell_type, "_", " "))

p_vals <- p_vals |>
  mutate(cell_type = str_replace_all(cell_type, "_", " "))


# filter on nb of cells

cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 30) |>
  pull(cell_type)


dotprod_agg_by_ct <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(
    p_adj = if_else(is.na(p_adj), 1, p_adj),
    tissue = factor(tissue,
                    levels = c("skin", "glia", "pharynx", "muscle",
                               "neuron", "reproductive", "other"))
  ) |>
  arrange(tissue, desc(mean_coherence)) |>
  mutate(cell_type = fct_inorder(cell_type))

dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  mutate(cell_type = factor(cell_type, levels = levels(dotprod_agg_by_ct$cell_type))) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence") + xlab(NULL) +
  # geom_hline(aes(yintercept = -.4)) + geom_hline(aes(yintercept = 1.4)) +
  # coord_cartesian(ylim = c(-.2,.9)) +
  scale_fill_manual(values = c(`TRUE` = "orange", `FALSE` = "white")) +
  geom_tile(aes(x = cell_type, y = diff(range(dotprod_by_cell$coherence))/2 + min(dotprod_by_cell$coherence),
                fill = p_adj < .05 ),
            height = diff(range(dotprod_by_cell$coherence)),
            alpha = .1,
            data = dotprod_agg_by_ct) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = coherence, color = tissue),
                               alpha = .2,
                               shape = 16) +
  geom_point(aes(x = cell_type, y = mean_coherence),
             data = dotprod_agg_by_ct)

gg_dotprod_by_cell
# gg copy: 1000x550
# 1300 x 450

# ggsave("local_phase_coherence.png", plot = gg_dotprod_by_cell,
#        path = "presentations/figures/local_phase_coherence",
#        width = 200, height = 70, units = "mm",
#        scale = 2)
# ggsave("local_phase_coherence.pdf", plot = gg_dotprod_by_cell,
#        path = "presentations/figures/local_phase_coherence",
#        width = 200, height = 70, units = "mm",
#        scale = 2)



# save plot with ggrastr

dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  mutate(cell_type = factor(cell_type, levels = levels(dotprod_agg_by_ct$cell_type))) |>
  ggplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm")
    ) +
  ylab("Local phase coherence") +
  xlab(NULL) +
  scale_fill_manual(
    values = c(`TRUE` = alpha("orange", alpha = .08), `FALSE` = "white")
    ) +
  geom_tile(
    aes(
      x = cell_type,
      y = diff(range(dotprod_by_cell$coherence))/2 + min(dotprod_by_cell$coherence),
      fill = p_adj < .05
    ),
    height = diff(range(dotprod_by_cell$coherence)),
    data = dotprod_agg_by_ct
  ) +
  ggrastr::geom_quasirandom_rast(
    aes(x = cell_type, y = coherence, color = tissue),
    alpha = .3,
    shape = 16,
    size = 1,
    width = .4,
    raster.dpi = 500
  ) +
  geom_point(
    aes(x = cell_type, y = mean_coherence),
    data = dotprod_agg_by_ct,
    size = 1
  )


# ggsave("local_phase_coherence.pdf",
#        path = "presentations/figures/local_phase_coherence",
#        width = 215, height = 75, units = "mm")



# save for other scripts
local_coherence_by_ct <- dotprod_by_cell |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(p_adj = if_else(is.na(p_adj), 1, p_adj)) |>
  as_tibble()

# qs::qsave(local_coherence_by_ct,
#           file.path(dir_out, "250610_coherence_by_ct.qs"))











# Plot GFP ----

samples_table <- read_tsv("data/samples_table.tsv")

# two ways to plot same thing

# FeaturePlot(seu, features = "nsIs198", pt.size = 2, alpha = .2, cols = c("bisque2", "green4"))

dat <- FetchData(seu, vars = c("orig.ident", "nsIs198", "umap_1", "umap_2", "cell_type"))



ggplot(dat) +
  theme_classic() +
  scale_color_gradient(low = "bisque2", high = "green4") +
  geom_point(aes(x = umap_1, y = umap_2, color = nsIs198),
             alpha = .2, size = 2)

dat2 <- left_join(dat, samples_table,
                  by = c(orig.ident = "sample_name"))

ggplot() +
  theme_void() +
  scale_color_gradient(low = "bisque2", high = "darkolivegreen") +
  scale_fill_gradient(low = "bisque2", high = "dodgerblue2") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "mm")
  ) +
  geom_point_rast(
    aes(x = umap_1, y = umap_2, color = nsIs198),
    alpha = .1, size = .8,
    raster.dpi = 500,
    data = filter(dat2, promoter == "grl-18")
  ) +
  geom_point_rast(
    aes(x = umap_1, y = umap_2, fill = nsIs198),
    alpha = .1, size = 1.5, shape = 21, stroke = NA,
    raster.dpi = 500,
    data = filter(dat2, promoter != "grl-18")
  )

# ggsave("umap_both_sorts_rastr.pdf",
#        path = "presentations/figures/250610_umap",
#        width = 56, height = 56, units = "mm")



# plot separately both sorts
dat2 |>
  filter(promoter == "grl-18") |>
  ggplot() +
  theme_classic() +
  scale_color_gradient(low = "bisque2", high = "green4") +
  geom_point(aes(x = umap_1, y = umap_2, color = nsIs198),
             alpha = .2, size = 2)


dat2 |>
  filter(promoter != "grl-18") |>
  ggplot() +
  theme_classic() +
  scale_color_gradient(low = "bisque2", high = "green4") +
  geom_point(aes(x = umap_1, y = umap_2, color = nsIs198),
             alpha = .2, size = 2)



# Some numbers for manuscript
dat2 |>
  count(promoter,
        nsIs198 > 0)

dat2 |>
  count(promoter,
        nsIs198 > 0,
        cell_type == "ILso")


dat2 |>
  filter(nsIs198 > 0,
         promoter == "grl-18") |>
  pull(cell_type) |>
  fct_lump_n(n = 5) |>
  enframe(value = "cell_type", name = NULL) |>
  count(cell_type) |>
  arrange(desc(n))




# UMAPs per tissue ----
tissues <- seu$tissue |> unique()
i=0

i=i+1
tissue_here <- tissues[[i]]
tissue_here

sub <- subset(seu, tissue == tissue_here)
table(sub$cell_type)
sub <- SCTransform(sub)|>
  RunPCA(npcs = 200, verbose = FALSE)

npca <- 40

ElbowPlot(sub, ndims = 200) +
  geom_vline(aes(xintercept = npca))


sub <- RunUMAP(sub, dims = 1:npca)

DimPlot(sub,
        group.by = "cell_type",
        label = TRUE) + NoLegend()

ggsave(paste0("celltype_",tissue_here,".png"), path = "presentations/figures/250611_umap_tissues/",
       width = 60, height = 60, units = "mm",
       scale = 2)
ggsave(paste0("celltype_",tissue_here,".pdf"), path = "presentations/figures/250611_umap_tissues/",
       width = 60, height = 60, units = "mm",
       scale = 2)


FetchData(sub,
          vars = c("umap_1", "umap_2",
                   "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  aes(x = umap_1, y = umap_2) +
  geom_point(aes(color = cell_phase_masked),
             shape = 16,
             size = 2,
             alpha = .2,
             show.legend = FALSE)

ggsave(paste0("phase_",tissue_here,".png"), path = "presentations/figures/250611_umap_tissues/",
       width = 60, height = 60, units = "mm",
       scale = 2)
ggsave(paste0("phase_",tissue_here,".pdf"), path = "presentations/figures/250611_umap_tissues/",
       width = 60, height = 60, units = "mm",
       scale = 2)


# ____________ ----
# Illustrations ----

#~ Individual cells ILso ----
sub <- qs::qread(
  file.path("intermediates/2502/250522_step1",
            paste0("ILso", "_seu_unsmoothed.qs"))
)


dat <- FetchData(sub,
                 vars = c("PC_1", "PC_2", "cell_phase_masked")) |>
  mutate(selected = FALSE,
         selected = {x <- selected; x[c(4,242)] <- TRUE; x})

dat$PC_1 <- -dat$PC_1


stopifnot(identical(
  rownames(dat)[c(4,242)],
  colnames(GetAssayData(sub, assay = "SCT", layer = "data"))[c(4,242)]
))



# PCA
DimPlot(sub,
        reduction = "pca",
        label = FALSE,
        pt.size = .8,
        alpha = .3,
        cells.highlight = rownames(dat)[c(4,242)],
        sizes.highlight = 3,
        cols.highlight = 'red') +
  NoLegend()


dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = PC_1, y = PC_2),
             alpha = .4, color = "grey") +
  geom_point(aes(x = PC_1, y = PC_2),
             data = dat |> filter(selected),
             size = 3, color = 'red3')

ggsave("phases_ILso_cells_4-242_pca.pdf", path = "presentations/250523_umap/",
       width = 6, height = 6, units = "in")
ggsave("phases_ILso_cells_4-242_pca.png", path = "presentations/250523_umap/",
       width = 6, height = 6, units = "in")


dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50)) +
  geom_point(aes(x = PC_1, y = PC_2, color = cell_phase_masked),
             size = 2,
             alpha = .7)

# ggsave("pca_ILso_color.png", path = "presentations/",
#        width = 6, height = 6, units = "in")


#~ cell, radial ----

# cells 4, 242

cell_nb <- 242
dat_1_cell <- enframe(mat[,cell_nb],
                      name = "gene_name",
                      value = "expression") |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(!is.na(osc_amplitude))



dat_1_cell |>
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
  geom_segment(aes(x = mean_angle,
                   xend = mean_angle,
                   y = 0,
                   yend = mean_rho),
               data = tibble(mean_angle = sub$cell_phase[[cell_nb]],
                             mean_rho = sub$cell_rho[[cell_nb]]),
               linewidth = 1,
               color = 'red3')


# ggsave(paste0("phases_ILso_cell_",cell_nb,".pdf"),
#        path = "presentations/",
#        width = 6, height = 6, units = "in")


dat_1_cell |>
  ggplot() +
  theme_bw() +
  coord_polar() +
  scale_x_continuous(limits = c(0,360), n.breaks = 15) +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50)) +
  ylab("expression") +
  geom_segment(aes(x = peak_phase_deg,
                   xend = peak_phase_deg,
                   y = 0,
                   yend = expression,
                   color = peak_phase_deg),
               linewidth = .25) +
  geom_segment(aes(x = mean_angle %% (360),
                   xend = mean_angle %% (360),
                   y = 0,
                   yend = mean_rho),
               data = tibble(mean_angle = sub$cell_phase[[cell_nb]],
                             mean_rho = sub$cell_rho[[cell_nb]]),
               linewidth = 1,
               color = 'black') #+theme(legend.position = 'none')


# ggsave(paste0("phases_ILso_cell_",cell_nb,"_col.pdf"),
#        path = "presentations/",
#        width = 6, height = 6, units = "in")

tibble(mean_angle = sub$cell_phase[[cell_nb]],
       mean_rho = sub$cell_rho[[cell_nb]],
       FDR = sub$length_FDR[[cell_nb]])









#~ BWM ----
sub <- subset(seu, cell_type == "BWM") |>
  SCTransform() |>
  RunPCA(npcs = 2, verbose = FALSE)

cell_nb <- 14


# note we need the full matrix here to ensure same normalization
mat <- LayerData(seu, assay = "SCT", layer = "data", features = osc_table$gene_name)
genes_max <- sparseMatrixStats::rowMaxs(mat)
genes_max[genes_max == 0] <- 1
mat <- mat / genes_max

mat <- mat[,colnames(sub)]


# PCA

dat <- FetchData(sub,
                 vars = c("PC_1", "PC_2", "cell_phase_masked")) |>
  mutate(selected = FALSE,
         selected = {x <- selected; x[c(cell_nb)] <- TRUE; x})
dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = PC_1, y = PC_2),
             alpha = .4, color = "grey") +
  geom_point(aes(x = PC_1, y = PC_2),
             data = dat |> filter(selected),
             size = 3, color = 'red3')

# ggsave(paste0("phases_BWM_cells_",cell_nb,"_pca.pdf"),
#        path = "presentations/",
#        width = 6, height = 6, units = "in")


# ggsave(paste0("phases_BWM_cells_",cell_nb,"_pca.png"),
#        path = "presentations/",
#        width = 6, height = 6, units = "in")


dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50)) +
  geom_point(aes(x = PC_1, y = PC_2, color = cell_phase_masked),
             size = 2,
             alpha = .7)

# ggsave("pca_BWM_color.png", path = "presentations/",
#        width = 6, height = 6, units = "in")



# 1 cell
cell_nb <- 14
dat_1_cell <- enframe(mat[,cell_nb],
                      name = "gene_name",
                      value = "expression") |>
  left_join(osc_table,
            by = "gene_name") |>
  filter(!is.na(osc_amplitude))


dat_1_cell |>
  ggplot() +
  theme_bw() +
  coord_polar() +
  scale_x_continuous(limits = c(0,360), n.breaks = 15) +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50)) +
  ylab("expression") +
  geom_segment(aes(x = peak_phase_deg,
                   xend = peak_phase_deg,
                   y = 0,
                   yend = expression,
                   color = peak_phase_deg),
               linewidth = .25) +
  geom_segment(aes(x = mean_angle %% (360),
                   xend = mean_angle %% (360),
                   y = 0,
                   yend = mean_rho),
               data = tibble(mean_angle = sub$cell_phase[[cell_nb]],
                             mean_rho = sub$cell_rho[[cell_nb]]),
               linewidth = 1,
               color = 'black') #+theme(legend.position = 'none')


# ggsave(paste0("phases_BWM_cell_",cell_nb,"_col.pdf"),
#        path = "presentations/",
#        width = 6, height = 6, units = "in")




# permutation test result illustr ----

stopifnot(all.equal(
  rownames(mat),
  osc_table$gene_name
))

mat2 <- mat[, cell_nb, drop = FALSE]



empirical <- rho_from_mat(mat2, osc_table$peak_phase_deg)

stopifnot( empirical == sub$cell_rho[cell_nb] )

perms <- replicate(n = 10000,
                   rho_from_mat(mat2, sample(osc_table$peak_phase_deg)))

as_tibble(perms) |>
  ggplot() +
  theme_classic() +
  # scale_x_continuous(limits = c(0, 1)) +
  xlab("Average phase length") +
  geom_histogram(aes(x = value),
                 bins = 50,
                 color = 'white',
                 linewidth = .3) +
  geom_vline(xintercept = empirical,
             color = 'red3',
             linewidth = 1.5)

# ggsave(paste0("perm_BWM_cell_",cell_nb,".pdf"),
#        path = "presentations/",
#        width = 100, height = 95, units = "mm")

sum(c(empirical,perms) >= empirical)
length(perms)




tibble(mean_angle = sub$cell_phase[[cell_nb]],
       mean_rho = sub$cell_rho[[cell_nb]],
       pval = sum(c(empirical,perms) >= empirical) / length(perms),
       FDR = sub$length_FDR[[cell_nb]])


#~ ILso permutations ----

sub <- subset(seu, cell_type == "ILso") |>
  SCTransform() |>
  RunPCA(npcs = 2, verbose = FALSE)

# note we need the full matrix here to ensure same normalization
mat <- LayerData(seu, assay = "SCT", layer = "data", features = osc_table$gene_name)
genes_max <- sparseMatrixStats::rowMaxs(mat)
genes_max[genes_max == 0] <- 1
mat <- mat / genes_max

mat <- mat[,colnames(sub)]



stopifnot(all.equal(
  rownames(mat),
  osc_table$gene_name
))


cell_nb <- 4

mat2 <- mat[, cell_nb, drop = FALSE]



empirical <- rho_from_mat(mat2, osc_table$peak_phase_deg)

stopifnot( empirical == sub$cell_rho[cell_nb] )

perms <- replicate(n = 10000,
                   rho_from_mat(mat2, sample(osc_table$peak_phase_deg)))

as_tibble(perms) |>
  ggplot() +
  theme_classic() +
  # scale_x_continuous(limits = c(0, 1)) +
  xlab("Average phase length") +
  geom_histogram(aes(x = value),
                 bins = 50,
                 color = 'white',
                 linewidth = .3) +
  geom_vline(xintercept = empirical,
             color = 'red3',
             linewidth = 1.5)

# ggsave(paste0("perm_ILso_cell_",cell_nb,".pdf"),
#        path = "presentations/",
#        width = 100, height = 95, units = "mm")

tibble(mean_angle = sub$cell_phase[[cell_nb]],
       mean_rho = sub$cell_rho[[cell_nb]],
       pval = sum(c(empirical,perms) >= empirical) / length(perms),
       FDR = sub$length_FDR[[cell_nb]])




# UMAP phases ----

FetchData(seu,
          vars = c("umap_1", "umap_2",
                   "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  scale_alpha_continuous(limits = c(0,1),
                         trans = scales::transform_exp()) +
  scale_fill_gradient(high = "black", low = "grey90",
                      limits = c(0,1),
                      trans = scales::transform_exp()) +
  aes(x = umap_1, y = umap_2) +
  geom_point(aes(fill = cell_rho),alpha = 0) +
  geom_point(aes(color = cell_phase_masked,
                 alpha = .2*cell_rho),
             shape = 16,
             size = 2,
             show.legend = FALSE)

# ggsave("umap_all.pdf", path = "presentations/",
#        width = 100, height = 70, units = "mm",
#        scale = 2)


# no alpha scale


FetchData(seu,
          vars = c("umap_1", "umap_2",
                   "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  # theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  aes(x = umap_1, y = umap_2) +
  geom_point(aes(color = cell_phase_masked),
             shape = 16,
             size = 2,
             alpha = .2,
             show.legend = FALSE)

# ggsave("umap_all.png", path = "presentations/",
#        width = 140, height = 120, units = "mm",
#        scale = 2)





FetchData(seu,
          vars = c("PC_1", "PC_2", "cell_type", "stage")) |>
  filter(cell_type == "AM_PHso") |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = stage),
             alpha = .2)

sub <- seu |>
  subset(cell_type == "AM_PHso") |>
  SCTransform(verbose = FALSE) |>
  RunPCA(npcs = 2, verbose = FALSE)


FetchData(sub,
          vars = c("PC_1", "PC_2", "cell_type", "stage")) |>
  filter(cell_type == "AM_PHso") |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = stage),
             alpha = .2)



FetchData(sub,
          vars = c("PC_1", "PC_2", "length_signif", "cell_phase")) |>
  mutate(cell_phase_masked = if_else(length_signif,
                                     cell_phase,
                                     NA)) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = cell_phase_masked),
             alpha = .2)

sub <- qs::qread(file.path(dir_out_individual_cts, "AM_PHso.qs"))


# compare older version ----
# the results
seu_new <- qs::qread(file.path(dir_out, "250605_seu_all_herma.qs"))
seu_old <- qs::qread("intermediates/2502/250509_assembled/250509_seu_all_herma.qs" )



# Find correspondences of sample indices by majority vote
correspondence <- inner_join(
  seu_old[[]] |>
    rownames_to_column("cell_bc") |>
    separate_wider_delim(cell_bc,
                         delim = "_",
                         names = c("sample", "cell_bc"),
                         too_many = "merge") |>
    select(sample, cell_bc, cell_type, tissue),
  seu_new[[]] |>
    rownames_to_column("cell_bc") |>
    separate_wider_delim(cell_bc,
                         delim = "_",
                         names = c("sample", "cell_bc"),
                         too_many = "merge") |>
    select(sample, cell_bc, cell_type, tissue),
  by = "cell_bc",
  relationship = "many-to-many"
) |>
  count(sample.x, sample.y, name = "shared_count") |>
  group_by(sample.x) |>
  slice_max(shared_count, n = 1, with_ties = FALSE) |>
  select(s_old = sample.x, s_new = sample.y)


merged <- full_join(
  seu_old[[]] |>
    rownames_to_column("cell_bc") |>
    separate_wider_delim(cell_bc,
                         delim = "_",
                         names = c("sample", "cell_bc"),
                         too_many = "merge") |>
    left_join(correspondence,
              by = c(sample = "s_old")) |>
    select(-sample) |> rename(sample = s_new),
  seu_new[[]] |>
    rownames_to_column("cell_bc") |>
    separate_wider_delim(cell_bc,
                         delim = "_",
                         names = c("sample", "cell_bc"),
                         too_many = "merge"),
  by = c("sample", "cell_bc")
) |>
  rename(cell_type_old = cell_type.x,
         cell_type_new = cell_type.y,
         tissue_old = tissue.x,
         tissue_new = tissue.y)

merged |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = tissue_old, y = tissue_new),
              alpha = .1)

merged |>
  filter(tissue_old == "other") |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = cell_type_old, y = cell_type_new),
              alpha = .1)


merged |>
  filter(tissue_old == "skin" | tissue_new == "skin") |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = cell_type_old, y = cell_type_new),
              alpha = .1)



merged |>
  filter(tissue_old == "glia" | tissue_new == "glia") |>
  ggplot() +
  theme_classic() +
  geom_jitter(aes(x = cell_type_old, y = cell_type_new),
              alpha = .1)



merged |>
  filter(tissue_old == "glia" | tissue_new == "glia",
         cell_type_old == "ADE_PDEso", is.na(cell_type_new)) |> 
  select(starts_with("stage")) |>
  count(stage.x, stage.y)



merged |>
  filter(tissue_old == "glia" | tissue_new == "glia",
         cell_type_old == "socket_s9" | cell_type_old == "socket_s7" |cell_type_new == "glia_sheath_2") |>
  count(stage.x, stage.y, cell_type_old, cell_type_new)





xx <- merged |>
  filter(cell_type_old == "OLso" | cell_type_new == "OLso") |>
  pull(cell_bc)

xx <- merged |>
  filter(cell_type_new == "CEPso") |>
  pull(cell_bc)

sub$tmp <- colnames(sub) %in% xx
DimPlot(sub, group.by = "tmp")


xx <- colnames(sub)[sub$seurat_clusters == 2]

merged |>
  filter(cell_bc %in% colnames(sub)[sub$seurat_clusters == 2]) |>
  count(stage.x, stage.y, cell_type_old, cell_type_new)





# ___________ ----

# ### Assemble everything (male and herma) ----
# 
# samples_table <- read_tsv("data/samples_table.tsv",
#                           col_types = "cffff") |>
#   mutate(file_path = file.path(dir_prefilt,
#                                sample_name |>
#                                  paste0(".qs")))
# 
# 
# 
# seu_list <- samples_table$file_path |>
#   map(qs::qread)
# 
# 
# seu <- merge(seu_list[[1]],
#              seu_list[-1])
# 
# rm(seu_list)
# 
# seu <- JoinLayers(seu)
# 
# 
# 
# #~ annotate metadata ----
# samples_table_lut <- samples_table |>
#   column_to_rownames("sample_name")
# 
# seu$sex <- samples_table_lut[seu$orig.ident, "sex"]
# seu$promoter <- samples_table_lut[seu$orig.ident, "promoter"]
# seu$stage <- samples_table_lut[seu$orig.ident, "stage"]
# 
# 
# 
# # process ----
# seu <- SCTransform(seu,
#                    conserve.memory = TRUE)
# 
# seu <- RunPCA(seu, npcs = 200, verbose = FALSE)
# 
# npca <- 40
# 
# ElbowPlot(seu, ndims = 200) +
#   geom_vline(aes(xintercept = npca))
# 
# # qs::qsave(seu, file.path(dir_out, "250328_seu_pca.qs"))
# 
# DimPlot(seu,
#         group.by = "sex",
#         reduction = "pca",
#         pt.size = 2,
#         alpha = .05) +
#   NoLegend()
# 
# 
# seu <- RunUMAP(seu, dims = 1:npca)
# 
# 
# DimPlot(seu,
#         group.by = "sex",
#         reduction = "umap",
#         pt.size = 2,
#         alpha = .05) +
#   NoLegend()

















