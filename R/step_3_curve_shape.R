# (on cluster)
# load the fitted object from step 2, analyze shape of curve


# Inits ----

message("Starting ", date())


if(Sys.getenv("CLUSTER") == "mccleary" || Sys.getenv("SLURM_CLUSTER_NAME") == "mccleary" ){
  
  message("cluster")
  in_pt_res_dir <- "intermediates/2502/250319_step2"
  out_pt_res_dir <- "intermediates/2502/250320_step3"
  options(wb_dir_cache = "/gpfs/ycga/work/hammarlund/aw853/references/WS295")
  
} else if(Sys.getenv("COMPUTERNAME") == "DESKTOP-C8IEHJQ"){
  
  message("personal computer")
  pt_res_dir <- "intermediates/240726_genes_by_celltype/241028_second_pass_cluster/"
  
} else{
  stop("computer not recognized")
}



# if mgcv not loaded, might call predict.lm instead of predict.gam for binomial fits
library(mgcv)

library(wbData)


gids <- rbind(
  wb_load_gene_ids(295),
  c(X = "a",
    gene_id = "GFP",
    symbol = "nsIs198",
    sequence = "nsIs198",
    status = "Live",
    biotype = "protein_coding_gene",
    name = "nsIs198"
  )
)







#' circular permutation to put peak in middle
#' @param preds matrix where each column is one gene and each row one timepoint, obtained
#' by predicting from GAM fit at regular intervals
#'
#' @description find max of each gene's time series, circular permute so the max is in the middle
circ_perm <- function(preds){
  
  peak_loc <- apply(preds, 2, which.max)
  
  preds_permuted <- matrix(nrow = nrow(preds), ncol = ncol(preds))
  t_tot <- nrow(preds)
  t_cent <- floor(t_tot/2)
  
  for(i in seq_len(ncol(preds))){
    
    if(peak_loc[[i]] == t_cent){
      
      preds_permuted[, i] <- preds[,i]
      
    } else if(peak_loc[[i]] == t_tot){
      
      preds_permuted[, i] <- c(preds[ (t_cent+1):t_tot, i],
                               preds[ 1:t_cent ,i])
      
    } else if(peak_loc[[i]] < t_cent){
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      a <- peak_loc[[i]]
      b <- t_cent + a
      
      preds_permuted[, i] <- c(preds[ (b+1):t_tot, i],
                               preds[ 1:b, i])
      
    } else{
      
      #  1     a      t_c     b      t_t
      #  |_____|_______|______|______|
      
      
      b <- peak_loc[[i]]
      a <- 1 + b - (t_cent+1)
      
      
      preds_permuted[, i] <- c(preds[ (a+1):b, i],
                               preds[ (b+1):t_tot, i],
                               preds[ 1:a, i])
    }
    
  }
  
  colnames(preds_permuted) <- colnames(preds)
  rownames(preds_permuted) <- rownames(preds)
  
  preds_permuted
}


#' Count number of peaks
#' 
#' @param y the prediction in which to find peaks
#' @param plot if FALSE, simply return nb of peaks, if TRUE, plot it
count_big_peaks <- function(y, plot = FALSE){
  
  #y[ !is.finite(y) ] <- 0
  
  x_peak <- which(scorepeak::detect_localmaxima(y, w = 9)) |>
    setdiff(c(1, length(y)))
  
  big <- which(y[x_peak] > (max(y) / 3) )
  
  x_big_peak <- x_peak[big]
  y_big_peak <- y[x_big_peak]
  
  
  if(plot == TRUE){
    plot(y, type = 'l')
    points(x_peak, y[x_peak], cex = 1.5, col = 'darkred')
    points(x_big_peak, y_big_peak, cex = 1.5, col = 'violet', pch = 19)
    abline(h = max(y) / 3, lty = 'dashed', col = 'grey80')
  }
  
  length(big)
}





cell_types <- list.files(in_pt_res_dir,
                         pattern = "^([a-zA-Z_ ]+)_res_pseudotimeDE.qs$") |>
  stringr::str_match("^([a-zA-Z_ ]+)_res_pseudotimeDE.qs$") |>
  (\(x) x[,2] )()


stopifnot(length(cell_types) > 1)

message("Number of cell types: ", length(cell_types))




# Loop ----
for(.ct in cell_types){
  
  message("Loading ", .ct)
  
  #~ load ----
  res <- qs::qread(file.path(in_pt_res_dir, paste0(.ct, "_res_pseudotimeDE.qs"))) |>
    dplyr::rename(gene_name = gene) |>
    dplyr::mutate(gene_id = s2i(gene_name, gids, warn_missing = TRUE),
           p_adj = p.adjust(para.pv, method = "BH"),
           .after = "gene_name")
  
  #~ extract gam.fit ----
  genes_to_test <- which(res$para.pv < 0.2)
  
  if(length(genes_to_test) == 0){
     message("    !! No gene significant in ", .ct,", skipping cell type !!")
     next
  }
  
  
  all_genes_id <- res$gene_id[genes_to_test]
  all_preds <- sapply(res$gam.fit[genes_to_test],
                           \(.mod) predict(.mod,
                                           type = "response",
                                           newdata = data.frame(pseudotime = (0:199)/200)))
  colnames(all_preds) <- all_genes_id
  
  #~ center ----
  all_preds_perm <- circ_perm(all_preds)
  
  
  #~ get curve properties ----
  # nb of peaks
  all_nb_peaks <- apply(all_preds_perm, 2, count_big_peaks)
  
  # area under curve normalized to amplitude
  all_amplitudes <- apply(all_preds_perm, 2, \(x) diff(range(x)))
  
  all_area_under_curve <- apply(all_preds_perm, 2,
                                \(.y){
                                  .y0 <- .y - min(.y)
                                  pracma::trapz(seq_along(.y0), .y0)
                                })
  
  all_area_ratio <- all_area_under_curve / all_amplitudes
  
  # height of fitted peak normalized to max cell count
  all_peak_max <- matrixStats::colMaxs(all_preds_perm)
  
  all_cell_max <- sapply(res$gam.fit[genes_to_test],
                         \(.mod) max(.mod$model$expv) )
  
  all_peak_norm_height <- all_peak_max / all_cell_max
  
  # model fit quality
  all_r2 <- sapply(res$gam.fit[genes_to_test],
                   \(.mod) summary(.mod)[["r.sq"]] )
  all_dev_expl <- sapply(res$gam.fit[genes_to_test],
                         \(.mod) summary(.mod)[["dev.expl"]] )
  
  #~ export ----
  
  # all res without gam.fit (heavy storage)
  res_exp <- res
  
  res_exp$cell_type <- .ct
  
  res_exp <- res_exp |>
    merge(data.frame(gene_id = all_genes_id,
                     nb_peaks = all_nb_peaks),
          by = "gene_id",
          all = TRUE) |>
    merge(data.frame(gene_id = all_genes_id,
                     area_ratio = all_area_ratio),
          by = "gene_id",
          all = TRUE) |>
    merge(data.frame(gene_id = all_genes_id,
                     peak_norm_height = all_peak_norm_height),
          by = "gene_id",
          all = TRUE) |>
    merge(data.frame(gene_id = all_genes_id,
                     r2 = all_r2),
          by = "gene_id",
          all = TRUE) |>
    merge(data.frame(gene_id = all_genes_id,
                     dev_expl = all_dev_expl),
          by = "gene_id",
          all = TRUE)
  
  # reorder
  res_exp <- res_exp[c("gene_id", "gene_name", "cell_type", "nb_peaks", "area_ratio",
                       "peak_norm_height", "r2", "dev_expl",
                       "para.pv",
                       "p_adj", "test.statistics", "zinf", "aic",
                       "expv.mean", "notes" )]
  message("   Saving")
  qs::qsave(res_exp,
            file.path(out_pt_res_dir, paste0(.ct, "_res_ptDE_preproc.qs")))
  
  
}


message("----------------    sessionInfo    --------------------------")

sessionInfo()



