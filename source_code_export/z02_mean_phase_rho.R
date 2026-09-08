# Efficiently compute rho (length of averaged vector) from angles in degree
rho_from_mat <- function(mat_expression, phases){
  
  phases_rad <- phases * pi/180
  
  sum_ws <- colSums(mat_expression)
  
  sines <- as.numeric( sin(phases_rad) %*% mat_expression ) / sum_ws
  cosines <- as.numeric( cos(phases_rad) %*% mat_expression ) / sum_ws
  
  rho <- sqrt( sines^2 + cosines^2 )
  rho[sum_ws == 0] <- 0
  
  rho * log10(sum_ws + 1)
  
}



# Efficiently compute angle (in degrees)
angle_from_mat <- function(mat_expression, phases){
  
  phases_rad <- phases * pi/180
  sum_ws <- colSums(mat_expression)
  
  sines <- as.numeric( sin(phases_rad) %*% mat_expression ) / sum_ws
  cosines <- as.numeric( cos(phases_rad) %*% mat_expression ) / sum_ws
  
  res_rad <- atan2( sines , cosines ) %% (2*pi)
  
  res_rad * 180 / pi
  
}



tau <- function(x){
  
  if(all(x == 0)) return(1)
  
  scaled <- x / max(x)
  N <- length(x)
  sum(1-scaled) / (N - 1)
}






#' Mean dotproduct of each cell's phase with its neighbors
#'
#' @param cell_phases a data.frame with rownames and two columns called `cell_phase` and `cell_rho`
#' @param seu.nn a seurat "nearest network" object, e.g. `seu@neighbors$RNA.nn`
#' @param k the number of nearest neighbors to use
#'
#' @return
#'
#' @examples
mean_dotprod <- function(cell_phases, seu.nn, k = 20){
  
  
  cells_phases_xy <- cell_phases |>
    mutate(x = cell_rho * cos(cell_phase *pi/180),
           y = cell_rho * sin(cell_phase *pi/180),
           x = if_else(is.nan(cell_phase), 0, x),
           y = if_else(is.nan(cell_phase), 0, y))
  
  
  cell_neighbors <- tibble(
    cell = rownames(cell_phases),
    neighbors = map(cell,
                    \(cell) TopNeighbors(seu.nn, cell = cell, n = (k+1L) ) |>
                      setdiff(cell),
                    .progress = TRUE)
  )
  
  
  
  
  all_cells_neighs_xy <- cell_neighbors |>
    mutate(cell_x = cells_phases_xy[cell, "x"],
           cell_y = cells_phases_xy[cell, "y"]) |>
    unnest(neighbors) |>
    mutate(neigh_x = cells_phases_xy[neighbors, "x"],
           neigh_y = cells_phases_xy[neighbors, "y"])
  
  
  mean_dotprod_by_cell <- tibble(cell = all_cells_neighs_xy$cell,
                                 dotprod = all_cells_neighs_xy$cell_x * all_cells_neighs_xy$neigh_x +
                                   all_cells_neighs_xy$cell_y * all_cells_neighs_xy$neigh_y) |>
    summarize(mean_dotprod = mean(dotprod),
              .by = cell)
  
  stopifnot(identical( mean_dotprod_by_cell$cell, rownames(cell_phases) ))
  
  mean_dotprod_by_cell$mean_dotprod
}







#' Mean dotproduct of each cell's phase with its neighbors
#'
#' @param cells a data.frame with rownames and three columns called `cell_phase`, `cell_rho`, and `cell_type`
#' @param seu.nn a seurat "nearest network" object, e.g. `seu@neighbors$RNA.nn`
#' @param k the number of nearest neighbors to use
#'
#' @return
#'
#' @examples
mean_dotprod_norm <- function(cells, seu.nn, k = 20){
  
  stopifnot(nchar(rownames(cells)[[1]]) > 1)
  stopifnot(all(
    c("cell_phase", "cell_rho", "cell_type") %in% colnames(cells)
  ))
  
  cells_phases_xy <- cells |>
    mutate(x = cell_rho * cos(cell_phase *pi/180),
           y = cell_rho * sin(cell_phase *pi/180),
           x = if_else(is.nan(cell_phase), 0, x),
           y = if_else(is.nan(cell_phase), 0, y))
  
  
  cell_neighbors_xy <- tibble(
    cell = rownames(cells),
    neighbors = map(cell,
                    \(cell) TopNeighbors(seu@neighbors$SCT.nn, cell = cell, n = (k+1L) ) |>
                      setdiff(cell),
                    .progress = TRUE)
  ) |>
    mutate(cell_x = cells_phases_xy[cell, "x"],
           cell_y = cells_phases_xy[cell, "y"]) |>
    unnest(neighbors) |>
    mutate(neigh_x = cells_phases_xy[neighbors, "x"],
           neigh_y = cells_phases_xy[neighbors, "y"])
  
  
  
  cell_type_xy <- tibble(
    cell = rownames(cells),
    neighbors = map(cells$cell_type,
                    \(ct) rownames(cells)[cells$cell_type == ct],
                    .progress = TRUE)
  ) |>
    mutate(cell_x = cells_phases_xy[cell, "x"],
           cell_y = cells_phases_xy[cell, "y"]) |>
    unnest(neighbors) |>
    filter(cell != neighbors) |>
    mutate(neigh_x = cells_phases_xy[neighbors, "x"],
           neigh_y = cells_phases_xy[neighbors, "y"])
  
  
  
  # dot products
  
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
  
  stopifnot(identical( cell_dotp_neighbors$cell, rownames(cells) ))
  
  # # Plot neighbors vs cell type
  # left_join(cell_dotp_neighbors |> rename(dotp_neigh = mean_dotprod),
  #           cell_dotp_ct |> rename(dotp_ct = mean_dotprod),
  #           by = "cell") |>
  #   left_join(cells |> rownames_to_column("cell"),
  #             by = "cell") |>
  #   ggplot() +
  #   theme_classic() +
  #   theme(legend.position = "none") +
  #   facet_wrap(~tissue) +
  #   geom_hline(aes(yintercept = 0), linetype = "dashed") +
  #   geom_vline(aes(xintercept = 0), linetype = "dashed") +
  #   geom_point(aes(x = dotp_ct, y = dotp_neigh, color = cell_type),
  #              alpha = .3)
  
  
  left_join(cell_dotp_neighbors |> rename(dotp_neigh = mean_dotprod),
            cell_dotp_ct |> rename(dotp_ct = mean_dotprod),
            by = "cell") |>
    left_join(cells |> rownames_to_column("cell"),
              by = "cell") |>
    mutate(coherence = log((dotp_neigh+1) / (dotp_ct+1)))
}






#' Title
#'
#' @param .perm permutation number; if 0 not permuted, otherwise permutation
#' @param cell_phases a data.frame with rownames and column called `cell_type`
#' @param mat expression matrix (scaled)
#' @param osc_table data.frame with column `peak_phase_deg` (bulk reference) in the same order as rows of `mat`
#' @param cell_neighbors data.frame with columns `cell` and `neighbors` giving neighbor cell names
#' @param cell_same_type data.frame with columns `cell` and `neighbors` giving cell names to normalize to
#'
#' @return
#' @export
#'
#' @examples
dotprods_normalized <- function(precomputed, mat, phase){
  
  
  cells <- precomputed$cell_phases
  
  # overwrite with permuted
  cells$cell_phase <- angle_from_mat(mat, phase)
  cells$cell_rho <- rho_from_mat(mat, phase)
  
  
  
  cells_phases_xy <- cells |>
    mutate(x = cell_rho * cos(cell_phase *pi/180),
           y = cell_rho * sin(cell_phase *pi/180),
           x = if_else(is.nan(cell_phase), 0, x),
           y = if_else(is.nan(cell_phase), 0, y))
  
  
  cell_neighbors <- precomputed$cell_neighbors
  cell_neighbors$cell_x <- cells_phases_xy$x[precomputed$neigh_cell_indices]
  cell_neighbors$cell_y <- cells_phases_xy$y[precomputed$neigh_cell_indices]
  cell_neighbors$neigh_x <- cells_phases_xy$x[precomputed$neigh_neigh_indices]
  cell_neighbors$neigh_y <- cells_phases_xy$y[precomputed$neigh_neigh_indices]
  
  
  
  
  cell_same_type <- precomputed$cell_same_type
  cell_same_type$cell_x <- cells_phases_xy$x[precomputed$type_cell_indices]
  cell_same_type$cell_y <- cells_phases_xy$y[precomputed$type_cell_indices]
  cell_same_type$neigh_x <- cells_phases_xy$x[precomputed$type_neigh_indices]
  cell_same_type$neigh_y <- cells_phases_xy$y[precomputed$type_neigh_indices]
  
  
  
  # dot products
  
  cell_dotp_neighbors <- cell_neighbors |>
    collapse::fmutate(dotprod = cell_x * neigh_x + cell_y * neigh_y) |>
    collapse::fgroup_by(cell) |>
    collapse::fsummarize(dotp_neigh = collapse::fmean(dotprod))
  
  
  cell_dotp_ct <- cell_same_type |>
    collapse::fmutate(dotprod = cell_x * neigh_x + cell_y * neigh_y) |>
    collapse::fgroup_by(cell) |>
    collapse::fsummarize(dotp_ct = collapse::fmean(dotprod))
  
  
  
  left_join(cell_dotp_neighbors,
            cell_dotp_ct,
            by = "cell") |>
    left_join(cells |> rownames_to_column("cell"),
              by = "cell") |>
    mutate(coherence = log((dotp_neigh+1) / (dotp_ct+1)))
  
}





run_permutation_test_by_celltype_and_phase <- function(.perm, precomputed, mat, phase){
  
  
  if(.perm == 0L){
    phase <- phase
  } else{
    phase <- sample(phase)
  }
  
  dotprods_normalized(precomputed, mat, phase) |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type)) |>
    add_column(permutation = .perm,
               .before = 1)
}








# Unnormalized ----
#' Title
#'
#' @param .perm permutation number; if 0 not permuted, otherwise permutation
#' @param cell_phases a data.frame with rownames and column called `cell_type`
#' @param mat expression matrix (scaled)
#' @param osc_table data.frame with column `peak_phase_deg` (bulk reference) in the same order as rows of `mat`
#' @param cell_neighbors data.frame with columns `cell` and `neighbors` giving neighbor cell names
#' @param cell_same_type data.frame with columns `cell` and `neighbors` giving cell names to normalize to
#'
#' @return
#' @export
#'
#' @examples
dotprods_unnormalized <- function(precomputed, mat, phase){
  
  
  cells <- precomputed$cell_phases
  
  # overwrite with permuted
  cells$cell_phase <- angle_from_mat(mat, phase)
  cells$cell_rho <- rho_from_mat(mat, phase)
  
  
  
  cells_phases_xy <- cells |>
    mutate(x = cell_rho * cos(cell_phase *pi/180),
           y = cell_rho * sin(cell_phase *pi/180),
           x = if_else(is.nan(cell_phase), 0, x),
           y = if_else(is.nan(cell_phase), 0, y))
  
  
  cell_neighbors <- precomputed$cell_neighbors
  cell_neighbors$cell_x <- cells_phases_xy$x[precomputed$neigh_cell_indices]
  cell_neighbors$cell_y <- cells_phases_xy$y[precomputed$neigh_cell_indices]
  cell_neighbors$neigh_x <- cells_phases_xy$x[precomputed$neigh_neigh_indices]
  cell_neighbors$neigh_y <- cells_phases_xy$y[precomputed$neigh_neigh_indices]
  
  
  
  
  # dot products
  
  cell_dotp_neighbors <- cell_neighbors |>
    collapse::fmutate(dotprod = cell_x * neigh_x + cell_y * neigh_y) |>
    collapse::fgroup_by(cell) |>
    collapse::fsummarize(dotp_neigh = collapse::fmean(dotprod))
  
  
  
  cell_dotp_neighbors |>
    left_join(cells |> rownames_to_column("cell"),
              by = "cell") |>
    mutate(coherence = dotp_neigh)
  
}


run_permutation_test_unnorm_by_celltype_and_phase <- function(.perm, precomputed, mat, phase){
  
  
  if(.perm == 0L){
    phase <- phase
  } else{
    phase <- sample(phase)
  }
  
  dotprods_unnormalized(precomputed, mat, phase) |>
    summarize(mean_coherence = mean(coherence),
              .by = c(tissue, cell_type)) |>
    add_column(permutation = .perm,
               .before = 1)
}

