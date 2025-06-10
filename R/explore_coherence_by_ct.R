


cells_phases_xy <- cells_phases |>
  mutate(x = cell_rho * cos(cell_phase *pi/180),
         y = cell_rho * sin(cell_phase *pi/180),
         x = if_else(is.nan(cell_phase), 0, x),
         y = if_else(is.nan(cell_phase), 0, y))


cell_neighbors_xy <- tibble(
  cell = rownames(cells_phases),
  neighbors = map(cell,
                  \(cell) TopNeighbors(seu@neighbors$SCT.nn, cell = cell, n = (20+1L) ) |>
                    setdiff(cell),
                  .progress = TRUE)
) |>
  mutate(cell_x = cells_phases_xy[cell, "x"],
         cell_y = cells_phases_xy[cell, "y"]) |>
  unnest(neighbors) |>
  mutate(neigh_x = cells_phases_xy[neighbors, "x"],
         neigh_y = cells_phases_xy[neighbors, "y"])



cell_type_xy <- tibble(
  cell = rownames(cells_phases),
  neighbors = map(cells_phases$cell_type,
                  \(ct) rownames(cells_phases)[cells_phases$cell_type == ct],
                  .progress = TRUE)
) |>
  mutate(cell_x = cells_phases_xy[cell, "x"],
         cell_y = cells_phases_xy[cell, "y"]) |>
  unnest(neighbors) |>
  filter(cell != neighbors) |>
  mutate(neigh_x = cells_phases_xy[neighbors, "x"],
         neigh_y = cells_phases_xy[neighbors, "y"])




cell_dotp_neighbors <- tibble(cell = cell_neighbors_xy$cell,
                              dotprod = cell_neighbors_xy$cell_x * cell_neighbors_xy$neigh_x +
                                cell_neighbors_xy$cell_y * cell_neighbors_xy$neigh_y) |>
  summarize(mean_dotprod = mean(dotprod),
            .by = cell)


cell_dotp_ct <- tibble(cell = cell_type_xy$cell,
                              dotprod = cell_type_xy$cell_x * cell_type_xy$neigh_x +
                         cell_type_xy$cell_y * cell_type_xy$neigh_y) |>
  summarize(mean_dotprod = mean(dotprod),
            .by = cell)

stopifnot(identical( cell_dotp_neighbors$cell, cell_dotp_ct$cell ))


stopifnot(identical( cell_dotp_neighbors$cell, rownames(cells_phases) ))


left_join(cell_dotp_neighbors |> rename(dotp_neigh = mean_dotprod),
          cell_dotp_ct |> rename(dotp_ct = mean_dotprod),
          by = "cell") |>
  left_join(cells_phases |> rownames_to_column("cell"),
            by = "cell") |>
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~tissue) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_point(aes(x = dotp_ct, y = dotp_neigh, color = cell_type),
             alpha = .3)


dotprod_by_cell <- left_join(cell_dotp_neighbors |> rename(dotp_neigh = mean_dotprod),
                             cell_dotp_ct |> rename(dotp_ct = mean_dotprod),
                             by = "cell") |>
  left_join(cells_phases |> rownames_to_column("cell"),
            by = "cell") |>
  mutate(norm_dotp = ((dotp_neigh+1) / (dotp_ct+1)))




# fancier plot
cell_types_to_plot <- dotprod_by_cell |>
  summarize(nb_cells = n(),
            .by = cell_type) |>
  filter(nb_cells >= 30) |>
  pull(cell_type) |>
  setdiff("reproductive")


dotprod_agg_by_ct <- dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  summarize(mean_coherence = mean(norm_dotp),
            .by = c(tissue, cell_type)) |>
  arrange(tissue, desc(mean_coherence)) |>
  mutate(cell_type = fct_inorder(cell_type))

dotprod_by_cell |>
  filter(cell_type %in% cell_types_to_plot) |>
  mutate(cell_type = factor(cell_type, levels = levels(dotprod_agg_by_ct$cell_type))) |>
  ggplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Local phase coherence norm") + xlab(NULL) +
  scale_y_continuous(labels = scales::label_number(.1)) +
  # geom_hline(aes(yintercept = -.2)) + geom_hline(aes(yintercept = .9)) +
  # coord_cartesian(ylim = c(0,.9)) +
  ggbeeswarm::geom_quasirandom(aes(x = cell_type, y = norm_dotp, color = tissue),
                               alpha = .2,
                               shape = 16) +
  geom_point(aes(x = cell_type, y = mean_coherence),
             data = dotprod_agg_by_ct)
















