# Inits ----


library(tidyverse)
library(Seurat)
library(wbData)

source("R/utils.R")
source("R/mean_phase_rho.R")

# SCTransform requires 1.3 GB for this data
options(future.globals.maxSize = 1.5 * 1024^3)


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







dir_out <- "intermediates/2502/250328_assembled"
dir_out_individual_cts <- file.path(dir_out, "250330_cell_types")

dir_third_processed <- "intermediates/2502/250313_third_processed"








# Just herma from third pass ----

#~ load ----


files_list <- list.files(dir_third_processed) |>
  enframe(value = "filename",
          name = NULL) |>
  mutate(filename2 = filename |>
           str_replace("seu_1","g1") |>
           str_replace("seu_2","g2")) |>
  filter(str_detect(filename2, "2[35]031[0-9]_(g1|g2)")) |>
  separate_wider_regex(filename2,
                       patterns = c(
                         "^2[35]031[0-9]_",
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




#~ tissue ----

ct2tissue <- read_csv("data/cell_type2tissue.csv") |>
  column_to_rownames("cell_type")

seu$tissue <- ct2tissue[ seu$cell_type , "tissue" ]



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




# Oscillations ----
osc_table <- osc_raw |>
  filter(gene_name %in% rownames(seu),
         Class == "Osc") |>
  select(gene_name, gene_id,
         osc_amplitude = OscAmplitude,
         peak_phase_deg = PeakPhase) |>
  group_by(gene_name) |>
  slice_sample(n = 1) |>
  ungroup()

# there are duplicated gene names, removed with slice_sample(); to look at them:
# duplicated_genes <- osc_table$gene_name[duplicated(osc_table$gene_name)]
# map(duplicated_genes,
#     ~ filter(osc_table, gene_name == .x))




#~ average phase and rho ----


mat <- LayerData(seu, layer = "data", features = osc_table$gene_name)



stopifnot(identical(rownames(mat),
                    osc_table$gene_name))
stopifnot(identical(colnames(mat),
                    rownames(seu[[]])))


seu$cell_phase <- angle_from_mat(mat, osc_table$peak_phase_deg)
seu$cell_rho <- rho_from_mat(mat, osc_table$peak_phase_deg)




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

# qs::qsave(pvals_by_cell, file.path(dir_out, "250329_permutations_10000.qs"))

pvals_by_cell <- qs::qread(file.path(dir_out, "250329_permutations_10000.qs"))


# hist(pvals_by_cell$p_val, breaks = 30, main = NULL, xlab = "Distribution of p-values")
# hist(pvals_by_cell$FDR, breaks = 50)

table(pvals_by_cell$FDR < .05)
#> FALSE  TRUE 
#>  6348 15733 


seu$length_FDR <- pvals_by_cell[rownames(FetchData(seu, vars = "ident")), "FDR"]
seu$length_signif <- (pvals_by_cell[rownames(FetchData(seu, vars = "ident")), "FDR"] < .05)
seu$cell_phase_masked <- if_else(seu$length_signif, seu$cell_phase, NA_real_)




FetchData(seu,
          vars = c("umap_1", "umap_2", "cell_phase_masked", "cell_rho")) |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pals::kovesi.cyclic_mrybm_35_75_c68(50),
                        limits = c(0, 360)) +
  geom_point(aes(x = umap_1, y = umap_2,
                 color = cell_phase_masked,
                 alpha = cell_rho),
             alpha = .1)








#~ Coherence within cell type ----



seu <- FindNeighbors(seu,
                     return.neighbor = TRUE,
                     reduction = "pca",
                     dims = 1:npca,
                     k.param = 21)







cells_phases <- FetchData(seu, vars = c("cell_phase", "cell_rho", "tissue", "cell_type"))





dotprod_by_cell <- cells_phases |>
  mutate(coherence = mean_dotprod(cells_phases,
                                  seu.nn = seu@neighbors$SCT.nn,
                                  k = 20))


# qs::qsave(dotprod_by_cell, file.path(dir_out, "250329_dotprod_by_cell.qs"))
dotprod_by_cell <- qs::qread(file.path(dir_out, "250329_dotprod_by_cell.qs"))




# Plot by cell type and cluster

dotprod_by_cell |>
  arrange(tissue, cell_type) |> mutate(cell_type = fct_inorder(cell_type)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence") + xlab(NULL) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = coherence, color = tissue)) +
  geom_point(aes(x = cell_type, y = mean_coherence),
             data = (
               dotprod_by_cell |>
                 summarize(mean_coherence = mean(coherence),
                           .by = "cell_type")
             ),
             size = 2)


# qs::qsave(seu, file.path(dir_out, "250329_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250329_seu_all_herma.qs"))










#~~ perm test ----

cells_phases <- FetchData(seu, vars = c("cell_phase", "cell_rho", "tissue", "cell_type"))

# the neighbors do not change between permutations: compute once and reuse
cell_neighbors <- tibble(
  cell = rownames(cells_phases),
  neighbors = map(cell,
                  \(cell) TopNeighbors(seu@neighbors$SCT.nn, cell = cell, n = (20+1L) ) |>
                    setdiff(cell),
                  .progress = TRUE)
)



run_permutation_test_by_celltype_rand_phase <- function(.perm){
  
  
  if(.perm > 0){
    
    osc_perm <- osc_table |>
      mutate(peak_phase_deg = sample(peak_phase_deg))
  } else{
    
    osc_perm <- osc_table
  }
  
  cells_phases_perm <- cells_phases
  
  # overwrite with permuted
  cells_phases_perm$cell_phase <- angle_from_mat(mat, osc_perm$peak_phase_deg)
  cells_phases_perm$cell_rho <- rho_from_mat(mat, osc_perm$peak_phase_deg)
  
  cells_phases_xy <- cells_phases_perm |>
    mutate(x = cell_rho * cos(cell_phase *pi/180),
           y = cell_rho * sin(cell_phase *pi/180),
           x = if_else(is.nan(cell_phase), 0, x),
           y = if_else(is.nan(cell_phase), 0, y))
  
  all_cells_neighs_xy <- cell_neighbors |>
    mutate(cell_x = cells_phases_xy[cell, "x"],
           cell_y = cells_phases_xy[cell, "y"]) |>
    unnest(neighbors) |>
    mutate(neigh_x = cells_phases_xy[neighbors, "x"],
           neigh_y = cells_phases_xy[neighbors, "y"])
  
  
  mean_dotprod_by_cell <- tibble(cell = all_cells_neighs_xy$cell,
                                 dotprod = all_cells_neighs_xy$cell_x * all_cells_neighs_xy$neigh_x +
                                   all_cells_neighs_xy$cell_y * all_cells_neighs_xy$neigh_y) |>
    summarize(coherence = mean(dotprod),
              .by = cell)
  
  cells_phases_perm |>
    add_column(coherence = mean_dotprod_by_cell$coherence) |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type)) |>
    add_column(permutation = .perm)
}



set.seed(123)
mean_dotprod_by_celltype_res_perm <- map_dfr(0:10000,
                                             run_permutation_test_by_celltype_rand_phase,
                                             .progress = TRUE)
# qs::qsave(mean_dotprod_by_celltype_res_perm,
#           file.path(dir_out, "250330_coherence_perm10000.qs"))

# mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_out, "250330_coherence_perm10000.qs"))


mean_dotprod_by_celltype_res_perm |>
  arrange(tissue) |> mutate(cell_type = fct_inorder(cell_type)) |>
  mutate(`permutated` = !(permutation == 0)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c('red3','grey')) +
  scale_alpha_manual(values = c(1,.01)) +
  ylab("Mean (in cell type) of mean (with neighbors) dot product") + xlab(NULL) +
  geom_point(aes(x = cell_type, y = mean_coherence, color = permutated, alpha = permutated),
             size = 2)



p_vals <- mean_dotprod_by_celltype_res_perm |>
  group_by(tissue, cell_type) |>
  nest() |>
  summarize(p_val = map_dbl(data,
                            \(dat){
                              mean(dat$mean_coherence >= dat$mean_coherence[[1]])
                            }),
            .groups = 'drop') |>
  mutate(p_adj = p.adjust(p_val, method = "holm"))


hist(p_vals$p_val)


# filter on nb of cells

cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 20) |>
  pull(cell_type)


gg_dotprod_by_cell <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  arrange(tissue, cell_type) |> mutate(cell_type = fct_inorder(cell_type)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence") + xlab(NULL) +
  scale_shape_manual(values = c("-","*")) +
  coord_cartesian(ylim = c(-.1,.7)) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = coherence, color = tissue),
                               alpha = .2) +
  geom_point(aes(x = cell_type, y = mean_coherence, shape = p_adj < .05),
             data = (
               dotprod_by_cell |>
                 filter(cell_type %in% cell_types_to_plot) |>
                 summarize(mean_coherence = mean(coherence),
                           .by = "cell_type") |>
                 left_join(p_vals,
                           by = c("cell_type")) |>
                 mutate(p_adj = if_else(is.na(p_adj), 1, p_adj))
             ),
             size = 8)

gg_dotprod_by_cell
# gg copy: 1000x550





# save for other scripts
local_coherence_by_ct <- dotprod_by_cell |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(p_adj = if_else(is.na(p_adj), 1, p_adj)) |>
  as_tibble()

# qs::qsave(local_coherence_by_ct,
#           file.path(dir_out, "250330_coherence_by_ct.qs"))





# Export cell types ----


# seu <- qs::qread( file.path(dir_out, "250329_seu_all_herma.qs"))



levels(Idents(seu)) |>
  walk(~{
    filename <- file.path(dir_out_individual_cts,
                          paste0(.x, ".qs"))
    
    sub <- subset(seu, idents = .x)
    
    
    if(!file.exists(filename)){
      qs::qsave(sub,
                filename)
      message("saved")
    } else{
      stop("exists!")
    }
  })










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

















