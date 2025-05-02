
# circular permutation of signal `sig` so that the center of `sig` ends up in `i`
circ_perm_i <- function(sig, i){
  
  n_t <- length(sig)
  m <- floor(n_t / 2)
  
  if(i == m){
    
    sig
    
  } else if(i < m){
    
    c(sig[(m - i + 1):(n_t)],
      sig[(1):(m - i)])
    
  } else{
    
    c(sig[(n_t - i + m + 1):(n_t)],
      sig[(1):(n_t - i + m)])
    
  }
}



# center signal so the peak ends up in the middle
# note might be off by one
circ_perm_mat <- function(preds){
  
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
