

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



#' Align two angular signals
#'
#' @param tb reference (in degrees)
#' @param tp experimental signal (in degrees)
#'
#' @returns a list with `invert` and `shift`, to transform the experimental signal into the reference
#'
#' @examples
#' opar <- par(no.readonly = T)
#' n <- 100
#' tb <- runif(n, 1, 360)
#' shift <- runif(1, 0, 360)
#' inv <- sample(c(T,F), 1)
#' tb_inv <- if(inv){
#'   360 - tb
#' } else{
#'    tb
#' }
#'     
#' tp <- ( tb_inv + shift ) %% 360 + rnorm(n, 0, 50)
#' inv; shift
#' res <- align_circular(tb, tp)
#'    
#' tp_aligned <- if (res$invert) {
#'   ((360 - tp) - res$shift) %% 360
#' } else {
#'     (tp - res$shift) %% 360
#' }
#'         
#' par(mfrow = c(1,2))
#' plot(tb, tp)
#' plot(tb, tp_aligned)
#' par(opar)
align_circular <- function(tb, tp) {
  
  shifts <- 0:359
  
  
  transformed <- list(
    original = sapply(shifts, function(s) (tp - s) %% 360),
    inverted = sapply(shifts, function(s) ((360 - tp) - s) %% 360)
  )
  
  
  
  rmses <- lapply(transformed, function(mat) {
    apply(mat, 2, function(x) {
      diff <- (x - tb + 180) %% 360 - 180
      sqrt(mean(diff^2))
    })
  })
  
  
  best_index <- which.min(c(rmses$original, rmses$inverted))
  invert <- best_index > 360
  shift <- shifts[(best_index - 1) %% 360 + 1]
  
  list(invert = invert, shift = shift)
  
}

## Test code




