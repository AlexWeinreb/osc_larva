

# circular permutation of signal `sig` so that the center of `sig` ends up in `i`
circ_perm <- function(sig, i){
  
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

# create a reference signal of length n and range 0-1
make_ref_sig <- function(n){
  sig <- dnorm(seq_len(n),
               mean = n/2,
               sd = n/6)
  sig <- sig - min(sig)
  sig <- sig / max(sig)
}

# a norm where we take average of summed column differences
mynorm <- function(mat){
  
  colSums(abs(mat)) |> mean()
}
