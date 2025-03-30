# Efficiently compute rho (length of averaged vector) from angles in degree
rho_from_mat <- function(mat_expression, phases){
  
  phases_rad <- phases * pi/180
  
  sum_ws <- colSums(mat_expression)
  
  sines <- as.numeric( sin(phases_rad) %*% mat_expression ) / sum_ws
  cosines <- as.numeric( cos(phases_rad) %*% mat_expression ) / sum_ws
  
  rho <- sqrt( sines^2 + cosines^2 )
  rho[sum_ws == 0] <- 0
  rho
  
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



















