# Unspliced-fraction in each cell type to distinguish cells from nuclei

# Inits ----
library(tidyverse)

dir_velocyto <- "intermediates/2502/250825_velocyto"

# Load ----

seu <- qs::qread("intermediates/2502/250605_assembled/250606_seu_all_herma.qs")

emat_tot <- qs::qread(file.path(dir_velocyto, "emat_tot.qs"))
nmat_tot <- qs::qread(file.path(dir_velocyto, "nmat_tot.qs"))


stopifnot(identical(
  colnames(emat_tot),
  colnames(nmat_tot)
))


# Filter ----
# Drop barcodes that didn't map to a Seurat cell (NA colnames from the join),
# then intersect with the cells actually present in `seu`.

valid <- !is.na(colnames(emat_tot)) & (colnames(emat_tot) %in% colnames(seu))

emat_filt <- emat_tot[, valid, drop = FALSE]
nmat_filt <- nmat_tot[, valid, drop = FALSE]


stopifnot(identical(
  colnames(emat_filt),
  colnames(nmat_filt)
))


# Fraction unspliced ----
# fraction per cell


per_cell <- tibble(
  bc_seu         = colnames(emat_filt),
  spliced        = Matrix::colSums(emat_filt),
  unspliced      = Matrix::colSums(nmat_filt),
  total          = spliced + unspliced,
  frac_unspliced = unspliced / total
) |>
  left_join(
    seu[[]] |> rownames_to_column("bc_seu") |> select(bc_seu, cell_type, percent.mt),
    by = "bc_seu"
  )

stopifnot(all( is.finite(per_cell$frac_unspliced) ))


# Plot ----

per_cell |>
  filter(n() >= 50, .by = cell_type) |>
  mutate(cell_type = fct_reorder(cell_type, frac_unspliced, .fun = median, .desc = TRUE)) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("black", "coral3")) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = frac_unspliced,
                                   color = cell_type == "hypodermis"),
                               alpha = .6,
                               show.legend = FALSE)




per_cell |>
  summarize(
    n             = n(),
    median_frac   = 100 * median(frac_unspliced),
    median_total  = median(total),
    median_mito   = median(percent.mt),
    .by = cell_type
  ) |>
  ggplot() +
  theme_classic() +
  scale_color_manual(values = c("black", "coral")) +
  aes(x = median_mito, y = median_total,
      color = cell_type == "hypodermis") +
  geom_point(show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(label = cell_type),
                           size = 3, max.overlaps = 10,
                           show.legend = FALSE)



