# Manually annotate genes to help clustering

library(Seurat)
library(tidyverse)


dir_step1 <- "intermediates/2502/250502_step1/"

dir_step2 <- "intermediates/2502/250502_step2/"
all_tests <- qs::qread(file.path(dir_step2, "mat_predictors_scaled.qs")) |>
  as.data.frame() |>
  rownames_to_column("cell_gene") |>
  separate_wider_delim(cell_gene,
                              delim = "|",
                              names = c("cell_type", "gene_name")) |>
  select(-s1,-s2,-s3,-s4)

cell_types <- list.files(dir_step1,
                         pattern = "_seu\\.qs$") |>
  str_remove("_seu\\.qs$")


# pick cell type
ct <- sample(cell_types, 1)
ct
subseu <- qs::qread( file.path(dir_step1, paste0(ct, "_seu.qs")) )


# pick genes
gois <- all_tests |>
  filter(cell_type == ct) |>
  filter(str_detect(gene_name, "\\-")) |>
  pull(gene_name) |> fct_inorder() |> sample(15) |> sort() |> as.character()

all_tests |>
  filter(cell_type == ct, gene_name %in% gois) |>
  select(1:2)
i=0

# look at them
i <- i+1
goi <- gois[[i]]

FetchData(subseu,
                vars = c("PC_1", "PC_2", goi)) |>
  set_names(c("PC_1", "PC_2", "count")) |>
  mutate(count = if_else(count < 1, NA, count )) |>
  ggplot() +
  theme_classic() +
  scale_color_gradient2(low = "grey50", mid = 'grey', high = "blue") +
  geom_point(aes(x = PC_1, y = PC_2,
                 color = count),
             size = 3,
             # shape = 16,
             alpha = .2) +
  ggtitle(bquote(italic(.(goi))))









# Clean duplicates
manual_annot <- readxl::read_excel("data/manual_annotations/250507_manual.xlsx") |>
  select(2:4)


manual_annot |> select(cell_type, gene_name) |> anyDuplicated()

manual_annot |>
  summarize(annots = list(manual),
            .by = c(cell_type, gene_name)) |>
  filter(lengths(annots) > 1)
  # unnest(annots)

manual_annot |>
  filter(manual == "yes" | manual == "no")




