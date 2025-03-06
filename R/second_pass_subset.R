library(tidyverse)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

options(future.globals.maxSize = 1000*1024^2)


dir_first <- "intermediates/2502/250305_per_condition"
dir_second <- "intermediates/2502/250306_per_condition_tissue"


#~ next group ----
group_here <- "4"


# Load first pass results ----
seu <- qs::qread(file.path(dir_first,
                           paste0("250305_seu_merged_",group_here,".qs")))

annot <- readxl::read_excel(file.path(dir_first,
                                      paste0("250305_annot_",group_here,".xlsx"))) |>
  column_to_rownames("cluster")

head(annot)
stopifnot(identical(levels(Idents(seu)), rownames(annot)))

seu$annot_first <- annot$annotation[seu$seurat_clusters]

table(seu$annot_first, useNA = 'ifany')



# Split by tissue ----

tissues <- unique(seu$annot_first) |> setdiff("unclear")
i <- 0

tissues

#~ next tissue ----
i <- i+1
tissue_here <- tissues[[i]]

tissue_here

sub <- subset(seu, annot_first == tissue_here)

sub


#~ Clustering ----
sub <- SCTransform(sub)


sub <- RunPCA(sub, npcs = 200, verbose = FALSE)


npca <- 50

ElbowPlot(sub, ndims = 200) +
  geom_vline(aes(xintercept = npca))


sub <- RunUMAP(sub,
               dims = 1:npca)

DimPlot(sub,
        reduction = "umap",
        group.by = "orig.ident",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend()




sub <- FindNeighbors(sub,
                     dims = 1:npca,
                     verbose = FALSE)

sub <- FindClusters(sub,
                    resolution = .1)


DimPlot(sub,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend()

tissue_here





all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)





file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") )

file.exists(file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") ))
file.exists(file.path( dir_second, paste0("250306_marks_",group_here,"_", tissue_here,".qs") ))


qs::qsave(sub,
          file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") ))
qs::qsave(all_markers,
          file.path( dir_second, paste0("250306_marks_",group_here,"_", tissue_here,".qs") ))






## With integration ----

table(sub$annot_first)
table(sub$orig.ident)


# underrepresented_batches <- names(table(sub$orig.ident))[table(sub$orig.ident) < 3]
# sub <- subset(sub, orig.ident %in% setdiff(unique(sub$orig.ident), underrepresented_batches))

sub <- split(sub, f = sub$orig.ident)



sub <- NormalizeData(sub)
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub)
sub <- RunPCA(sub, npcs = 200, verbose = FALSE)


npca <- 50

ElbowPlot(sub, ndims = 200) +
  geom_vline(aes(xintercept = npca))


# sub <- IntegrateLayers(object = sub,
#                        method = CCAIntegration)
# sub0 <- sub

sub <- IntegrateLayers(object = sub,
                       method = HarmonyIntegration)



sub <- RunUMAP(sub,
               dims = 1:npca,
               reduction = "harmony",
               reduction.name = "umap")



DimPlot(
  sub,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


sub <- FindNeighbors(sub,
                     dims = 1:npca,
                     reduction = "harmony")

sub <- FindClusters(sub,
                    cluster.name = "seurat_clusters",
                    resolution = .8)

DimPlot(
  sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()




sub <- JoinLayers(sub)


all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)





file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") )

file.exists(file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") ))
file.exists(file.path( dir_second, paste0("250306_marks_",group_here,"_", tissue_here,".qs") ))


qs::qsave(sub,
          file.path( dir_second, paste0("250306_seu_",group_here,"_", tissue_here,".qs") ))
qs::qsave(all_markers,
          file.path( dir_second, paste0("250306_marks_",group_here,"_", tissue_here,".qs") ))












