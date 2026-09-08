

library(tidyverse)
library(Seurat)
library(wbData)

gids <- wb_load_gene_ids(295) |>
  tibble::add_row(X = NA, gene_id = "nsIs198",
                  symbol = "nsIs198", sequence = "nsIs198",
                  status = "Live", biotype = "protein_coding_gene",
                  name = "nsIs198")

source("R/utils.R")


dir_save <- "intermediates/2502/250305_per_condition"


rev_markers <- readxl::read_excel("/gpfs/ycga/work/hammarlund/aw853/references/data/tissue_markers_250305.xlsx",
                                  sheet = "selected") |>
  mutate(gene_id = s2i(gene, gids, warn_missing = TRUE) )


stopifnot(!any(duplicated(rev_markers$gene)))



## Load ----
subname <- "4"


seu <- qs::qread(file.path(dir_save,
                           paste0("250305_seu_merged_",subname,".qs")))

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
        alpha = .2) +
  NoLegend()



# reverse approach ----


aggregate_tib(seu,
              gene_names = rev_markers$gene[1:57]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)

aggregate_tib(seu,
              gene_names = rev_markers$gene[58:117]) |>
  ggplot() +
  theme_classic() +
  scale_size_continuous(limits = c(.01,1)) +
  scale_color_gradient2(low = "blue3",
                        mid = "red3",
                        high = "yellow2",
                        midpoint = 2) +
  geom_point(aes(x = cell_id, y = gene_name,
                 size = prop, color = log10(count)),
             alpha = .8)




# plot annotation ----

annot <- readxl::read_excel(file.path(dir_save,
                                      paste0("250305_annot_",subname,".xlsx"))) |>
  column_to_rownames("cluster") |>
  mutate(cell_type = if_else(is.na(cell_type), annotation, cell_type))

head(annot)

seu$annot_first <- annot$annotation[seu$seurat_clusters]
seu$annot_more <- annot$cell_type[seu$seurat_clusters]

DimPlot(seu,
        reduction = "umap",
        group.by = "annot_first",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend() +
  ggtitle(NULL)


DimPlot(seu,
        reduction = "umap",
        group.by = "annot_more",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend() +
  ggtitle(NULL)


seu$tmp <- seu$seurat_clusters == "8"

DimPlot(seu,
        reduction = "umap",
        group.by = "tmp",
        label = TRUE,
        pt.size = 2,
        alpha = .2) +
  NoLegend() +
  ggtitle(NULL)


# check single cluster ----

marks <- FindMarkers(seu,
                     ident.1 = "29") |>
  rownames_to_column("gene") |>
  as_tibble()
marks |>
  mutate(delta = pct.1 - pct.2) |>
  arrange(desc(delta)) |>
  filter(p_val_adj < .1,
         pct.1 > .3,
         pct.2 < .1) |> #pull(gene) |> paste0(collapse = ",") |> message()
  slice_head(n = 10)

FeaturePlot(seu,
            features = c("C32D5.8", "T01B4.3"),
            pt.size = 1.5,
            alpha = .3)

FeaturePlot(seu,
            features = c("hsp-16.2","hsp-16.11"),
            pt.size = 2,
            alpha = .2)



# compare to subset of other clusters
marks <- FindMarkers(seu,
                     ident.1 = "4",
                     ident.2 = c("1","3","9")) |>
  rownames_to_column("gene") |>
  as_tibble()





# recluster subset ----

sub <- subset(seu, idents = "24")

sub <- SCTransform(sub)

sub <- RunPCA(sub, npcs = 100, verbose = FALSE)

ElbowPlot(sub, ndims = 200)


sub <- FindNeighbors(sub,
                     dims = 1:30,
                     verbose = FALSE)

sub <- FindClusters(sub,
                    resolution = .8)

sub <- RunUMAP(sub,
               dims = 1:30)

DimPlot(sub,
        reduction = "umap",
        group.by = "orig.ident",
        label = TRUE,
        pt.size = 2,
        alpha = .5) +
  NoLegend()

DimPlot(sub,
        reduction = "umap",
        label = TRUE,
        pt.size = 2,
        alpha = .5) +
  NoLegend()




sub <- PrepSCTFindMarkers(sub)



all_markers <- FindAllMarkers(sub,
                              logfc.threshold = 1,
                              min.pct = .3) |>
  as_tibble() |>
  rename(gene_id = gene) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))


all_markers |>
  filter(cluster == "3") |>
  filter(p_val_adj < 0.05,
         abs(pct.1 - pct.2) > 0.3,
         avg_log2FC > 0) |>
  pull(gene_name) |> paste0(collapse = " ")


FeaturePlot(sub,
            features = c("egl-20", "col-118", "C39E9.8","W03F9.11") |> s2i(gids))








