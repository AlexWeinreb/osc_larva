# third pass: annotate neurons

# two subpasses:
# * first, for each seurat_cluster, use majority vote to attribute ident; if needed split clust
# * second, for each label, check if subclusters


library(Seurat)
library(tidyverse)


dir_third <- "intermediates/2502/250311_third"
dir_processed <- "intermediates/2502/250313_third_processed"




# Load ----
cond_here <- "1"
tissue_here <- "neuron"

seu <- qs::qread(file.path( dir_third, paste0("250311_seu_",cond_here,"_", tissue_here,".qs") ))
# seu <- qs::qread(file.path(dir_processed, "230314_seu_neurons.qs"))



DimPlot(seu,
        reduction = "umap",
        group.by = "orig.ident",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend()


DimPlot(seu,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .5) +
  NoLegend()


DimPlot(seu,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .5,group.by = "second_cell_type") +
  NoLegend()

DimPlot(seu,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .5,group.by = "third_neuron_type") +
  NoLegend()




# +++ First part +++ ----


clust_nb <- 0




#~ next subcluster ----
clust_nb <- clust_nb + 1; clust_nb

sub <- subset(seu, idents = clust_nb)
ncol(sub)
# table(sub$second_cell_type)
# 
# table(sub$third_neuron_type)





# recluster
sub <- SCTransform(sub)

nps_max <- pmin(200, ncol(sub) - 2L)
sub <- RunPCA(sub, npcs = nps_max, verbose = FALSE)

npca <- 4

ElbowPlot(sub, ndims = nps_max) +
  geom_vline(aes(xintercept = npca))







sub <- RunUMAP(sub,
               dims = 1:npca,
               n.neighbors = pmin(30, ncol(sub)-2),
               reduction = "pca",
               reduction.name = "umap")


DimPlot(
  sub,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .2
) +
  NoLegend()








#~ save ----

table(sub$second_cell_type)
(maj <- table(sub$second_cell_type) |> sort(decreasing = TRUE) |> head(1) |> names())



seu$third_neuron_type[colnames(seu) %in% colnames(sub)] <- maj



table(seu$third_neuron_type, useNA = 'ifany')

# qs::qsave(seu, file.path(dir_processed, "230318_seu_neurons1.qs"))





#~~ save manual ----
seu$third_neuron_type[colnames(seu) %in% colnames(sub)] <- ""


#~~ to discard ----
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)] <- "unclear"

# also discard the clusters with too few cells:
# seu$third_neuron_type[is.na(seu$third_neuron_type)] <- "unclear"



#~~ save subclusters ----

# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "0"]] <- "RID"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "1"]] <- "RMG"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "2"]] <- "AVB"

# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "0"]] <- "ASER"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "1"]] <- "ASEL"




# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "0"]] <- "ASJ"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "1"]] <- "ASI"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "2"]] <- "ADF"


# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "0"]] <- "RIP"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "1"]] <- "AIM"


# group 1
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "0"]] <- "unclear"
# seu$third_neuron_type[colnames(seu) %in% colnames(sub)[Idents(sub) == "1"]] <- "mechanosensory"




# check ----

genes <- c("F31A9.2","ver-4")

genes <- genes[genes %in% rownames(seu)]

genes

FeaturePlot(sub,
            genes, alpha=.5, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()





#~ markers from whole seurat object ----
# 
# table(sub$doubletFinder, Idents(sub), useNA = 'ifany')
# 
# 
# 
# all_markers <- qs::qread( file.path(dir_third, "250311_marks_1_neuron.qs") )

all_markers |>
  filter(cluster == clust_nb) |>
  filter(p_val_adj < 0.1,
         abs(pct.1 - pct.2) > 0.3,
         pct.1 > .2,
         pct.2 < .2,
         avg_log2FC > 1) |> pull(gene_name) |> paste(collapse = ",") |> message()






#~~ markers subclust ----



sub <- FindNeighbors(sub,
                     dims = 1:npca,
                     reduction = "pca")

sub <- FindClusters(sub,
                    cluster.name = "seurat_clusters",
                    resolution = .1)

DimPlot(
  sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


markers_sub <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)

markers_sub |>
  filter(cluster == "1") |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 1) |> pull(gene_name) |> paste(collapse = ",") |> message()



#~ special cases ----

table(Idents(seu)[seu$second_cell_type == "SMB"]) |> sort()
sub <- subset(seu, idents = c("1","14","43"))



table(Idents(seu)[seu$second_cell_type == "AVA"]) |> sort()
sub <- subset(seu, idents = c("0","24"))

table(Idents(seu)[seu$second_cell_type == "ASJ"]) |> sort()
sub <- subset(seu, idents = c("30","18"))





# +++ Second part +++ ----
seu <- qs::qread(file.path(dir_processed, "230318_seu_neurons1.qs"))

Idents(seu) <- "third_neuron_type"
all_neurs <- levels(Idents(seu))

all_neurs

clust_nb <- 0




#~ next subcluster ----
clust_nb <- clust_nb + 1; neur <- all_neurs[[clust_nb]]; neur

sub <- subset(seu, idents = neur)


table(sub$third_neuron_type)
table(sub$seurat_clusters)[table(sub$seurat_clusters) > 0]





# recluster
sub <- SCTransform(sub)

nps_max <- pmin(200, ncol(sub) - 3L)
sub <- RunPCA(sub, npcs = nps_max, verbose = FALSE)

npca <- 4

ElbowPlot(sub, ndims = nps_max) +
  geom_vline(aes(xintercept = npca))







sub <- RunUMAP(sub,
               dims = 1:npca,
               n.neighbors = pmin(30, ncol(sub)-2),
               reduction = "pca",
               reduction.name = "umap")


DimPlot(
  sub,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .2
) +
  NoLegend()











#~ Check ----


marks <- FindMarkers(seu,
                     ident.1 = "OLQ") |>
  rownames_to_column("gene_name") |>
  as_tibble()




marks |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0)




# +++ save results by neur +++ ----

# seu <- qs::qread(file.path(dir_processed, "230314_seu_neurons2.qs"))

levels(Idents(seu)) |>
  setdiff("unclear") |>
  walk(~{
    filename <- file.path(dir_processed,
                          paste0("230318_g1_", .x, ".qs"))
    
    sub <- subset(seu, idents = .x)
    
    
    if(!file.exists(filename)){
      qs::qsave(sub,
                filename)
      message("saved")
    } else{
      stop("exists!")
    }
  })











