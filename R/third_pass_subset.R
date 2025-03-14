library(tidyverse)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")



dir_second <- "intermediates/2502/250306_per_condition_tissue"
dir_third <- "intermediates/2502/250311_third"




annots <- readxl::read_excel(file.path( dir_second, "250306_annot_tissues.xlsx" ),
                             col_types = "text")

seu_objects <- list.files(dir_second, pattern = "^250306_seu_.+qs$") |>
  enframe(value = "filename",
        name = NULL) |>
  separate_wider_regex(cols = filename,
                       patterns = c("^250306_", "seu_",
                                    group = "[1-4]",
                                    "_",
                                    annot_first = "[a-z]+",
                                    "\\.qs$"),
                       cols_remove = FALSE)

stopifnot(identical(
  seu_objects |>
    select(group, annot_first) |>
    arrange(group, annot_first),
  annots |>
    select(group, annot_first) |>
    distinct() |>
    arrange(group, annot_first)
))



group_here <- 0



# Next condition ----
group_here <- group_here+1

group_here


all_seu <- seu_objects |>
  filter(group == group_here) |>
  rename(.group = group,
         .annot_first = annot_first,
         .filename = filename) |>
  pmap(
    \(.group, .annot_first, .filename){
      
      seu <- file.path(dir_second, .filename) |>
        qs::qread()
      
      lut <- annots |>
        filter(group == .group,
               annot_first == .annot_first) |>
        select(cluster, annot_second, cell_type) |>
        column_to_rownames("cluster")
      
      seu$second_tissue <- lut$annot_second[ Idents(seu) ]
      seu$second_cell_type <- lut$cell_type[ Idents(seu) ]
      
      seu
    }
  )

# can't merge when too few cells in a batch
all_seu <- map(all_seu,
                \(sub){
                  
                  underrepresented_batches <- names(table(sub$orig.ident))[table(sub$orig.ident) < 3]
                  sub <- subset(sub, orig.ident %in% setdiff(unique(sub$orig.ident), underrepresented_batches))
                  sub
                })

seu <- merge(all_seu[[1]],
             all_seu[-1])
Idents(seu) <- "second_tissue"

DefaultAssay(seu) <- "RNA"
seu[["SCT"]] <- NULL

seu <- JoinLayers(seu)

# qs::qsave(seu, file.path(dir_third, "250311_group_4_seu.qs"))

table(Idents(seu), useNA = "ifany")

tissues <- levels(Idents(seu)) |> setdiff("unclear")

tissues


i <- 0


#~ next tissue ----
i <- i+1
tissue_here <- tissues[[i]]

tissue_here

sub <- subset(seu, second_tissue == tissue_here)


sub



# No integration ----

table(sub$annot_first)
table(sub$orig.ident)


sub <- SCTransform(sub)

nps_max <- pmin(200, ncol(sub) - 10L)
sub <- RunPCA(sub, npcs = nps_max, verbose = FALSE)

npca <- 60

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






sub <- FindNeighbors(sub,
                     dims = 1:npca,
                     reduction = "pca")

#~ clust ----
sub <- FindClusters(sub,
                    cluster.name = "seurat_clusters",
                    resolution = 10)

DimPlot(
  sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


#~ save ----
file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") )

file.exists(file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") ))
file.exists(file.path( dir_third, paste0("250311_marks_",group_here,"_", tissue_here,".qs") ))


qs::qsave(sub,
          file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") ))

if(length(levels(Idents(sub))) > 1){
  all_markers <- FindAllMarkers(sub,
                                logfc.threshold = 1,
                                min.pct = .3) |>
    as_tibble() |>
    rename(gene_name = gene)
  
  qs::qsave(all_markers,
            file.path( dir_third, paste0("250311_marks_",group_here,"_", tissue_here,".qs") ))
} else{
  message("Single cluster, skipping markers")
}









# With integration ----
sub <- subset(seu, second_tissue == tissue_here)

table(sub$annot_first)
table(sub$orig.ident)


underrepresented_batches <- names(table(sub$orig.ident))[table(sub$orig.ident) < 3]
sub <- subset(sub, orig.ident %in% setdiff(unique(sub$orig.ident), underrepresented_batches) )


table(sub$orig.ident)


sub <- split(sub, f = sub$orig.ident)




#~ Normalization ----
sub <- SCTransform(sub)

nps_max <- pmin(200, ncol(sub) - 10L)
sub <- RunPCA(sub, npcs = nps_max, verbose = FALSE)

npca <- 30

ElbowPlot(sub, ndims = nps_max) +
  geom_vline(aes(xintercept = npca))





sub <- IntegrateLayers(object = sub,
                       method = HarmonyIntegration)


sub <- RunUMAP(sub,
               dims = 1:npca,
               n.neighbors = pmin(30, ncol(sub)-2),
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

#~ clust ----
sub <- FindClusters(sub,
                    cluster.name = "seurat_clusters",
                    resolution = .4)

DimPlot(
  sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 2,
  alpha = .2
) +
  NoLegend()


sub <- PrepSCTFindMarkers(sub)



#~ save ----
file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") )

file.exists(file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") ))
file.exists(file.path( dir_third, paste0("250311_marks_",group_here,"_", tissue_here,".qs") ))


qs::qsave(sub,
          file.path( dir_third, paste0("250311_seu_",group_here,"_", tissue_here,".qs") ))

if(length(levels(Idents(sub))) > 1){
  all_markers <- FindAllMarkers(sub,
                                logfc.threshold = 1,
                                min.pct = .3) |>
    as_tibble() |>
    rename(gene_name = gene)
  
  qs::qsave(all_markers,
            file.path( dir_third, paste0("250311_marks_",group_here,"_", tissue_here,".qs") ))
} else{
  message("Single cluster, skipping markers")
}



















# check ----

genes <- genes[genes %in% rownames(sub)]

genes

FeaturePlot(sub,
            genes, alpha=.7, pt.size = 1.5,
            combine = FALSE) |># lapply(\(.x) .x + NoLegend()) |>
  patchwork::wrap_plots()



all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_name = gene)
all_markers |>
  filter(pct.1 > .3,
         pct.2 < .1,
         p_val_adj < .1) |> View()










#~ Compare ----

integ
nonint

marks_nonint <- FindAllMarkers(nonint,
                               logfc.threshold = 1,
                               min.pct = .4) |>
  as_tibble() |>
  rename(gene_name = gene)

marks_integ <- FindAllMarkers(integ,
                              logfc.threshold = 1,
                              min.pct = .4) |>
  as_tibble() |>
  rename(gene_name = gene)



pb_integ <- AggregateExpression(integ, assays = "RNA")[["RNA"]]
pb_nonint <- AggregateExpression(nonint, assays = "RNA")[["RNA"]]

cc <- cor(as.matrix(pb_integ), as.matrix(pb_nonint))
cc

printMat::matimage(cc)



DimPlot(integ, label = TRUE) + NoLegend()
DimPlot(nonint, label = TRUE) + NoLegend()

idents <- full_join(
  integ[[]] |>
    rownames_to_column("cell_bc") |>
    select(cell_bc, seurat_clusters) |>
    rename(clust_integ = seurat_clusters)
  ,
  nonint[[]] |>
    rownames_to_column("cell_bc") |>
    select(cell_bc, seurat_clusters) |>
    rename(clust_nonint = seurat_clusters)
  ,
  by = "cell_bc"
)

table(idents$clust_integ, idents$clust_nonint) |> printMat::matimage()

table(idents$clust_integ, idents$clust_nonint) |>
  pheatmap::pheatmap(cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     scale = "column")



xx <- marks_nonint |> filter(cluster == "7", pct.2 < .05) |> pull(gene_name) |> head(6)
xx <- marks_integ |> filter(cluster == "17", pct.2 < .05) |> pull(gene_name) |> head(6)




FeaturePlot(integ,
            features = xx,
            alpha=.7, pt.size = 1)

FeaturePlot(nonint,
            features = xx,
            alpha=.7, pt.size = 1)



