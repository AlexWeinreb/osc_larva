
# Inits ----

library(tidyverse)
library(Seurat)

library(wbData)


gids <- wb_load_gene_coords(295)


data_dir <- "intermediates/2502/250304_doubletFinder"
dir_out <- "intermediates/2502/250304_filt_ds"



sample_table <- tibble(
  sample_filename = list.files(data_dir, pattern = "^250304_.+_dblts\\.qs$")
) |>
  separate_wider_regex(sample_filename,
                       patterns = c(
                         "^250304_",
                         sample_id = "[[:alnum:]_]+",
                         "_dblts\\.qs"
                       ),
                       cols_remove = FALSE)


datasets <- set_names(sample_table$sample_filename,
                      sample_table$sample_id) |>
  map(~ qs::qread(file.path(data_dir, .x)) )





# QC ----

map(datasets,
    ~ rownames(GetAssayData(.x,
                            assay="RNA",
                            slot="counts"))) %>%
  set_names(paste0("s_",seq_along(.))) |>
  UpSetR::fromList() |>
  UpSetR::upset(nset = 10)

common_genes <- map(datasets,
                    ~ rownames(GetAssayData(.x,
                                            assay="RNA",
                                            slot="counts"))) |> 
  reduce(intersect)



mt_genes <- gids |>
  filter(name %in% common_genes) |>
  filter(chr == "MtDNA") |>
  pull(name)
mt_genes

gids$name[ str_detect(gids$name, pattern = "^rps|^rpl") ]
ribo_genes <- gids |>
  filter(name %in% common_genes) |>
  filter(str_detect(name, pattern = "^rps|^rpl")) |>
  pull(name)
ribo_genes


for(i in 1:length(datasets)){
  datasets[[i]][["percent.mt"]] <- PercentageFeatureSet(datasets[[i]],
                                                        features = mt_genes)
  datasets[[i]][["percent_ribo"]] <- PercentageFeatureSet(datasets[[i]],
                                                          features = ribo_genes)
  
}



#~ Some global QC ----

xx <- imap(datasets,
           ~ ggplot(slice_sample(.x@meta.data, n=10000)) +
             geom_point(aes(x=nFeature_RNA, y=percent.mt), alpha=.2) +
             theme_classic() +
             scale_x_log10() +
             geom_hline(aes(yintercept=25), color="red") +
             geom_vline(aes(xintercept=100), color="red") +
             xlab("Number of genes") +
             ylab("% genes from mitochondria") +
             ggtitle(.y))
patchwork::wrap_plots(xx)



xx <- imap(datasets,
           ~ ggplot(slice_sample(.x@meta.data, n=10000)) +
             geom_point(aes(x=nCount_RNA, y=percent.mt), alpha=.2) +
             theme_classic() +
             scale_x_log10() +
             geom_hline(aes(yintercept=25), color="red") +
             xlab("RNA count") +
             ylab("% genes from mitochondria") +
             ggtitle(.y))
patchwork::wrap_plots(xx)



xx <- imap(datasets,
           ~ ggplot(slice_sample(.x@meta.data, n=1000)) +
             geom_point(aes(x=nCount_RNA, y=nFeature_RNA), alpha=.2) +
             theme_classic() +
             scale_x_log10() +
             scale_y_log10() +
             geom_hline(aes(yintercept=100), color="red") +
             xlab("RNA count") +
             ylab("Nb genes") +
             ggtitle(.y))
patchwork::wrap_plots(xx)


xx <- imap(datasets,
           ~ ggplot(slice_sample(.x@meta.data, n=1000)) +
             geom_point(aes(x=nCount_RNA, y=percent_ribo), alpha=.2) +
             theme_classic() +
             scale_x_log10() +
             xlab("RNA count") +
             ylab("% ribosomal genes") +
             ggtitle(.y))
patchwork::wrap_plots(xx)





# qs::qsave(datasets, "intermediates/2502_prefilt/250213_datasets.qs")





#~ Per-sample QC ----

# filtered_datasets <- vector(mode = "list", length = length(datasets))
# names(filtered_datasets) <- names(datasets)
# i <- 0

# Run those steps for each sample
i <- i+1
i
ds <- datasets[[i]]
nds <- names(datasets)[[i]]
nds

thrs_genes <- 70
thrs_mito <- 15
patchwork::wrap_plots(
  
  ggplot(slice_sample(ds@meta.data, n=10000)) +
    geom_point(aes(x=nFeature_RNA, y=percent.mt), alpha=.2) +
    theme_classic() +
    scale_x_log10() +
    geom_hline(aes(yintercept=thrs_mito), color="red") +
    geom_vline(aes(xintercept=thrs_genes), color="red") +
    xlab("Number of genes") +
    ylab("% genes from mitochondria"),
  
  ggplot(slice_sample(ds@meta.data, n=10000)) +
    geom_point(aes(x=nFeature_RNA, y=percent_ribo), alpha=.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(aes(xintercept=thrs_genes), color="red") +
    xlab("Number of genes") +
    ylab("% genes from ribosomes"),
  
  ggplot(slice_sample(ds@meta.data, n=10000)) +
    geom_point(aes(x=nCount_RNA, y=nFeature_RNA), alpha=.2) +
    theme_classic() +
    scale_x_log10() +
    geom_hline(aes(yintercept=thrs_genes), color="red") +
    xlab("RNA count") +
    ylab("Nb genes")
  
)

ncol(ds)

ds <- subset(ds, percent.mt < thrs_mito & nFeature_RNA > thrs_genes)
ds

ds <- NormalizeData(ds)
ds <- FindVariableFeatures(ds)
ds <- ScaleData(ds)
ds <- RunPCA(ds, npcs = 100, verbose = FALSE)
ElbowPlot(ds, ndims = 100)

# Clustering
ds <- FindNeighbors(ds, dims = 1:30)
ds <- FindClusters(ds,
                   resolution = 2)

table(Idents(ds))

# ds <- RunUMAP(ds, dims = 1:30)
# DimPlot(ds, label = TRUE) + NoLegend()


as_tibble(ds@meta.data) %>%
  ggplot() +
  geom_violin(aes(x=seurat_clusters, y=nCount_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(x=seurat_clusters, y=nCount_RNA),
                               alpha = .5) +
  scale_y_log10() +
  theme_classic() +
  ylab("Number of UMIs")


as_tibble(ds@meta.data) %>%
  ggplot() +
  geom_violin(aes(x=seurat_clusters, y=nFeature_RNA)) +
  ggbeeswarm::geom_quasirandom(aes(x=seurat_clusters, y=nFeature_RNA),
                               alpha = .5) +
  scale_y_log10() +
  theme_classic() +
  ylab("Number of features")


as_tibble(ds@meta.data) %>%
  ggplot() +
  geom_boxplot(aes(x=seurat_clusters, y=percent.mt)) +
  # geom_violin(aes(x=seurat_clusters, y=percent.mt)) +
  # ggbeeswarm::geom_quasirandom(aes(x=seurat_clusters, y=percent.mt),
  #                              alpha = .5) +
  theme_classic() +
  ylab("Percent mitochondrial genes")



as_tibble(ds@meta.data) %>%
  ggplot() +
  geom_boxplot(aes(x=seurat_clusters, y=percent_ribo)) +
  theme_classic() +
  ylab("Percent ribosomal genes")


# look for clusters that do not have a clear identity
xx <- FindAllMarkers(ds,logfc.threshold = 1)
# xx |>
#   group_by(cluster) |>
#   slice_head(n=20) |>
#   mutate(gene_name = i2s(gene, gids)) |>
#   View()

ds$seurat_clusters |> table()
# more convenient to copy-paste
cl <- 0



cl <- cl+6
cl:(cl+5)
xx |>
  filter(cluster %in% cl:(cl+5),
         avg_log2FC > .1,
         pct.1 > .1,
         p_val_adj < .05) |>
  group_by(cluster) |>
  slice_head(n=20) |>
  pull(gene) |>
  paste0(collapse = ",") |> message()

xx |>
  filter(cluster %in% cl:(cl+5),
         avg_log2FC > .1,
         pct.1 > .1,
         p_val_adj < .05) |>
  group_by(cluster) |> summarize(n = n())





# individual cluster

xx |>
  filter(cluster == "13",
         avg_log2FC > 1.5,
         pct.1 > .1,
         p_val_adj < .1) |>
  pull(gene) |>
  head(150) |>
  paste0(collapse = ",") |>
  message()

xx |>
  filter(cluster == "8",
         avg_log2FC >= 0,
         pct.1 >= pct.2) |>
  mutate(gene = i2s(gene, gids)) |>
  arrange(desc(pct.1))

i; nds
#1: 200730_batch1_CHB3840b,           all good
#2: 201013_batch2_CHB3840b_CEG_fqs    all good
#3: 201013_batch2_CHB3841_CEG_fqs     all good
#4: 210413_batch3_CHB3840b            all_good
#5: 210413_batch3_CHB3841           cl 0 bad
#6: 210420_batch4_CHB3841           cl 3 bad
#7: 210420_batch4_hmn                 all good
#8: 210427_batch5_CHB3840b          cl3 bad
#9: 210427_batch5_CHB3841           cl3 bad
#10: 220210_OH17400                 all ok
#11: 230414                         removing cl1 (maybe reproductive)
#12: 230421_AM                      cl1 bad
#13: 230421_PM                      all ok
#14: 230505_AM                      all ok (cl2 few markers, but may be pharynx)
#15: 230505_PM                      all ok
#16: 231019GRL                      all ok
#17: 240111                         all ok (cl17 repro)





# ds <- subset(ds,
#              idents = setdiff(
#                levels(Idents(ds)),
#                c("1")
#                ))

i
nds
file.exists(file.path( dir_out, paste0(nds,".qs") ))
qs::qsave(ds, file.path( dir_out, paste0(nds,".qs") ))

# rm(datasets)

#~~ end ----













