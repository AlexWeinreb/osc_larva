# Inits ----


library(tidyverse)
library(Seurat)
library(wbData)

source("R/utils.R")
source("R/mean_phase_rho.R")

# SCTransform requires 1.3 GB for this data
options(future.globals.maxSize = 1.5 * 1024^3)



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







dir_out <- "intermediates/2502/250509_assembled"
dir_out_individual_cts <- file.path(dir_out, "250330_cell_types")

dir_third_processed <- "intermediates/2502/250313_third_processed"


# the result
# seu <- qs::qread( file.path(dir_out, "250509_seu_all_herma.qs"))

# FeaturePlot(seu, features = "nsIs198", pt.size = 2, alpha = .2, cols = c("bisque2", "green4"))



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

seu$tissue2 <- seu$tissue
seu$tissue2[seu$cell_type == "ILso"] <- "ILso"
DimPlot(seu,
        group.by = "tissue2",
        reduction = "umap",
        label = FALSE,
        pt.size = 2,
        alpha = .1) +
  NoLegend()

# ggsave("UMAP_tissue.png", path = "presentations/",
#        width = 110, height = 120, units = "mm",
#        scale = 2)




# Oscillations ----
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

qs::qsave(pvals_by_cell, file.path(dir_out, "250509_permutations_10000.qs"))

pvals_by_cell <- qs::qread(file.path(dir_out, "250509_permutations_10000.qs"))


# hist(pvals_by_cell$p_val, breaks = 30, main = NULL, xlab = "Distribution of p-values")
# hist(pvals_by_cell$FDR, breaks = 50)

table(pvals_by_cell$FDR < .05)
#> FALSE  TRUE 
#>  9791 12290


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
                 color = cell_phase_masked),
             alpha = .3)









#~ Local coherence index ----



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


# qs::qsave(dotprod_by_cell, file.path(dir_out, "250509_dotprod_by_cell.qs"))
dotprod_by_cell <- qs::qread(file.path(dir_out, "250509_dotprod_by_cell.qs"))



# Plot by cell type and cluster

dotprod_by_cell |>
  arrange(tissue, cell_type) |> mutate(cell_type = fct_inorder(cell_type)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim = c(-.1, .6)) +
  ylab("Local phase coherence") + xlab(NULL) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = coherence, color = tissue)) +
  geom_point(aes(x = cell_type, y = mean_coherence),
             data = (
               dotprod_by_cell |>
                 summarize(mean_coherence = mean(coherence),
                           .by = "cell_type")
             ),
             size = 2)


# qs::qsave(seu, file.path(dir_out, "250509_seu_all_herma.qs"))
# seu <- qs::qread( file.path(dir_out, "250509_seu_all_herma.qs"))






# fancier plot
cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 30) |>
  pull(cell_type) |>
  setdiff("reproductive")


dotprod_agg_by_ct <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  summarize(mean_coherence = mean(coherence),
            .by = c(tissue, cell_type)) |>
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
#           file.path(dir_out, "250509_coherence_perm10000.qs"))

# mean_dotprod_by_celltype_res_perm <- qs::qread(file.path(dir_out, "250509_coherence_perm10000.qs"))


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
hist(p_vals$p_adj)


# filter on nb of cells

cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 30) |>
  pull(cell_type) |>
  setdiff("reproductive")


dotprod_agg_by_ct <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(p_adj = if_else(is.na(p_adj), 1, p_adj),
         signif = cut(p_adj, breaks = c(-Inf, 1e-3,1e-2,5e-2,Inf), labels = c("***","**","*","n.s."))) |>
  arrange(tissue, desc(mean_coherence)) |>
  mutate(cell_type = fct_inorder(cell_type))

gg_dotprod_by_cell <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
   mutate(cell_type = factor(cell_type, levels = levels(dotprod_agg_by_ct$cell_type))) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence") + xlab(NULL) +
  # geom_hline(aes(yintercept = -.2)) + geom_hline(aes(yintercept = .9)) +
  # coord_cartesian(ylim = c(-.2,.9)) +
  geom_tile(aes(x = cell_type, y = .8713012,
                fill = p_adj < .05 ),
            height = 2.212274,
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
#        path = "presentations/",
#        width = 200, height = 70, units = "mm",
#        scale = 2)





# save for other scripts
local_coherence_by_ct <- dotprod_by_cell |>
  summarize(mean_coherence = mean(coherence),
            .by = "cell_type") |>
  left_join(p_vals,
            by = c("cell_type")) |>
  mutate(p_adj = if_else(is.na(p_adj), 1, p_adj)) |>
  as_tibble()

# qs::qsave(local_coherence_by_ct,
#           file.path(dir_out, "250509_coherence_by_ct.qs"))

















# UMAPs per tissue ----
tissues <- seu$tissue |> unique()
i=0

i=i+1
tissue_here <- tissues[[i]]
tissue_here

sub <- subset(seu, tissue == tissue_here)
table(sub$cell_type)
sub <- SCTransform(sub)
sub <- RunPCA(sub, npcs = 200, verbose = FALSE)

npca <- 20

ElbowPlot(sub, ndims = 200) +
  geom_vline(aes(xintercept = npca))


sub <- RunUMAP(sub, dims = 1:npca)

DimPlot(sub,
        group.by = "cell_type",
        label = TRUE) + NoLegend()

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

# ggsave("umap_other.png", path = "presentations/",
#        width = 78, height = 60, units = "mm",
#        scale = 2)



# ____________ ----
# Illustrations ----

#~ Individual cells ----
# seu <- sub
mat <- GetAssayData(seu)

DimPlot(seu,
        reduction = "umap",
        label = FALSE,
        pt.size = .8,
        alpha = .3,
        cells.highlight = colnames(mat)[c(10, 23, 6)],
        sizes.highlight = 3,
        cols.highlight = 'red') +
  NoLegend()


dat <- FetchData(seu,
          vars = c("umap_1", "umap_2")) |>
  mutate(selected = FALSE,
         selected = {x <- selected; x[c(10, 23, 6)] <- TRUE; x})
dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = umap_1, y = umap_2),
             alpha = .1, color = "grey") +
  geom_point(aes(x = umap_1, y = umap_2),
             data = dat |> filter(selected),
             size = 3, color = 'red3')



#~ cell 18 ----
# 10, 23, 6
dat_1_cell <- enframe(mat[,242],
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
  geom_segment(aes(x = mean_angle %% (360),
                   xend = mean_angle %% (360),
                   y = 0,
                   yend = mean_rho),
               data = dat_1_cell |>
                 summarize(mean_angle = circhelp::weighted_circ_mean(peak_phase_deg*pi/180, expression)*180/pi,
                           mean_rho = max(expression)*circhelp::weighted_circ_rho(peak_phase_deg*pi/180, expression)),
               linewidth = 1,
               color = 'red3')



ggsave("phases_ILso_cell_242.pdf", path = "presentations/",
       width = 6, height = 6, units = "in")


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
               data = dat_1_cell |>
                 summarize(mean_angle = circhelp::weighted_circ_mean(peak_phase_deg*pi/180, expression)*180/pi,
                           mean_rho = max(expression)*circhelp::weighted_circ_rho(peak_phase_deg*pi/180, expression)),
               linewidth = 1,
               color = 'black') #+theme(legend.position = 'none')


ggsave("phases_ILso_cell_242_col.pdf", path = "presentations/",
       width = 6, height = 6, units = "in")



# same on PCA instead of UMAP

# seu <- sub
mat <- GetAssayData(seu)

DimPlot(seu,
        reduction = "pca",
        label = FALSE,
        pt.size = .8,
        alpha = .3,
        cells.highlight = colnames(mat)[c(1,242)],
        sizes.highlight = 3,
        cols.highlight = 'red') +
  NoLegend()


dat <- FetchData(seu,
                 vars = c("PC_1", "PC_2")) |>
  mutate(selected = FALSE,
         selected = {x <- selected; x[c(4,242)] <- TRUE; x})
dat |> 
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = PC_1, y = PC_2),
             alpha = .1, color = "grey") +
  geom_point(aes(x = PC_1, y = PC_2),
             data = dat |> filter(selected),
             size = 3, color = 'red3')

# ggsave("phases_ILso_cells_4-242_pca.pdf", path = "presentations/",
#        width = 6, height = 6, units = "in")



# permutation test result illustr ----


mat2 <- mat[,colnames(sub)[[7]], drop = FALSE]

stopifnot(identical(rownames(mat2),
                    osc_table$gene_name))

sub$cell_rho[[7]]

empirical <- rho_from_mat(mat2, osc_table$peak_phase_deg)
perms <- replicate(n = 10000,
                   rho_from_mat(mat2, sample(osc_table$peak_phase_deg)))

as_tibble(perms) |>
  ggplot() +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  xlab("Average phase length") +
  geom_histogram(aes(x = value),
                 bins = 50,
                 color = 'white',
                 linewidth = .3) +
  geom_vline(xintercept = empirical,
             color = 'red3',
             linewidth = 1.5)

ggsave("perm_BWM_cell_7.pdf", path = "presentations/",
       width = 120, height = 75, units = "mm")

sum(perms >= empirical)
length(perms)



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

















