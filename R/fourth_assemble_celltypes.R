
# Inits ----
library(tidyverse)
library(Seurat)







dir_third_processed <- "intermediates/2502/250313_third_processed"


# Load ----

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

files_list |>
  count(cell_type) |>
  arrange(desc(n)) |>
  count(n)


files_list_by_type <- files_list |>
  summarize(both_names = list(filename),
            n = n(),
            .by = cell_type) |>
  filter(n == 2) |>
  select(-n)


i <- 0

#~ Next cell type ----

i <- i+1; i


.ct <- files_list_by_type$cell_type[[i]]
.files <- files_list_by_type$both_names[[i]]


.ct

stopifnot(length(.files) == 2L)

seu <- merge(
  qs::qread(file.path( dir_third_processed, .files[[1]] )),
  qs::qread(file.path( dir_third_processed, .files[[2]] ))
)

seu


# reprocess
seu <- SCTransform(seu)

nps_max <- pmin(200, ncol(seu) - 2L)
seu <- RunPCA(seu, npcs = nps_max, verbose = FALSE)

npca <- 30

ElbowPlot(seu, ndims = nps_max) +
  geom_vline(aes(xintercept = npca))







seu <- RunUMAP(seu,
               dims = 1:npca,
               n.neighbors = pmin(30, ncol(seu)-2),
               reduction = "pca",
               reduction.name = "umap")


DimPlot(
  seu,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 2,
  alpha = .2
) +
  NoLegend()





DimPlot(seu,
        reduction = "pca",
        group.by = "orig.ident",
        pt.size = 2,
        alpha = .2) +
  NoLegend()














