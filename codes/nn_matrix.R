library(dplyr)
library(spNNGP)

get_NN_ind <- function(ind, ind_distM_i, M){
  l <- ifelse(ind < M, ind, M)
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}

get_nn <- function(coords, m){
  n <- nrow(coords)
  nn_data <- spConjNNGP(
    rep(0, n) ~ 1, coords = coords,
    n.neighbors = m,
    theta.alpha = c("phi" = 5, "alpha" = 0.5),
    sigma.sq.IG = c(2, 1),
    cov.model = "exponential",
    return.neighbor.info = T, fit.rep = F, 
    verbose = F)
  ord <- nn_data$neighbor.info$ord
  nn_idx <- sapply(1:(n - 1), get_NN_ind, nn_data$neighbor.info$n.indx[-1], m) %>%
    t()
  
  return(list(ord = ord, nn_idx = nn_idx))
}
