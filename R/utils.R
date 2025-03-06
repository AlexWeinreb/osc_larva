aggregate_tib <- function(seu, gene_names, convert_id = FALSE){
  if(convert_id){
    gene_ids <- s2i(gene_names, gids, warn_missing = TRUE)
  } else{
    gene_ids <- gene_names
  }
  
  not_found <- which( ! gene_ids %in% rownames(seu) )
  
  if( length(not_found) > 0){
    message("Not found: ", gene_names[not_found])
    gene_ids <- gene_ids[-not_found]
    gene_names <- gene_names[-not_found]
  }
  
  all_dat <- FetchData(seu, vars = gene_ids)
  cell_ids <- Idents(seu)
  
  lapply(unique(cell_ids) |> rlang::set_names(),
         \(cid){
           dat_cell <- all_dat[cell_ids == cid, ]
           
           tibble::tibble(cell_id = cid,
                          gene_id = colnames(all_dat),
                          gene_name = gene_names,
                          count = colSums(dat_cell),
                          prop = colMeans(dat_cell > 0))
         }) |>
    purrr::list_rbind() |>
    mutate(gene_name = fct_inorder(gene_name))
  
}


AggregateProportion <- function(obj){
  
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  
  idents_onehot <- Matrix::sparse.model.matrix(~0+ident,
                                               data = data.frame( ident = Idents(obj)))
  colnames(idents_onehot) <- sub("^ident", "", colnames(idents_onehot))
  
  count_per_ident <- (counts > 0) %*% idents_onehot
  
  nb_cells <- colSums(idents_onehot)
  
  prop_per_ident <- apply(count_per_ident, 1, \(x) {x / nb_cells} ) |> t()
  prop_per_ident
}
