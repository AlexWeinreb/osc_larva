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



plotsum <- function(genes, seu, assay = "RNA"){
  
  genes <- genes |> intersect(rownames(seu))
  
  if(length(genes) == 0L) stop("No gene found")
  
  mat <- GetAssayData(seu, assay = assay)[genes,]
  colsums <- colMeans(mat)
  umap <- FetchData(seu, vars = c("umap_1", "umap_2"))
  umap$expr_genes <- colsums[rownames(umap)]
  umap$expr_genes[umap$expr_genes == 0] <- NA_real_
  
  ggplot(umap) +
    theme_minimal() +
    scale_color_gradient(low = 'darkblue', high = 'red3') +
    geom_point(aes(x = umap_1, y = umap_2,
                   color = expr_genes),
               size = 2,
               alpha = .2)
}


plotsum2 <- function(genes1, genes2, seu, assay = "RNA"){
  
  genes1 <- genes1 |> intersect(rownames(seu))
  genes2 <- genes2 |> intersect(rownames(seu))
  
  if(length(genes1) == 0L | length(genes2) == 0L) stop("No gene found")
  
  mat1 <- GetAssayData(seu, assay = assay)[genes1,]
  mat2 <- GetAssayData(seu, assay = assay)[genes2,]
  
  colsums1 <- colMeans(mat1)
  colsums2 <- colMeans(mat2)
  
  umap <- FetchData(seu, vars = c("umap_1", "umap_2"))
  
  umap$expr_genes1 <- colsums1[rownames(umap)]
  umap$expr_genes2 <- colsums2[rownames(umap)]
  umap$blend <- umap$expr_genes1 - umap$expr_genes2
  
  ggplot(umap) +
    theme_minimal() +
    scale_color_gradient2(low = scales::muted('red'),
                          high = scales::muted('blue'),
                          mid = 'lightgrey') +
    geom_point(aes(x = umap_1, y = umap_2,
                   color = blend),
               size = 2,
               alpha = .8)
}


