library(dplyr)
library(spNNGP)

# Distance matrix for location i and its neighbors
i_dist <- function(i, neighbor_idx, coords_ord){
  idx <- c(i, neighbor_idx[[i - 1]])
  d <- dist(coords_ord[idx, ])
  return(d)
}

get_dist_matrix <- function(neighbor_dist, neighbor_idx){
  n <- length(neighbor_dist)
  D <- matrix(0, nrow = n, ncol = n)
  
  for (i in 2:length(neighbor_idx)) {
    l <- length(neighbor_idx[[i - 1]]) # number of neighbors
    d <- neighbor_dist[[i]][1:l]
    D[i, neighbor_idx[[i - 1]]] <- d
    # D[neighbor_idx[[i - 1]], i] <- d
  }
  return(D)
}

get_NN_ind <- function(ind, ind_distM_i, M){
  if (ind < M ) {
    l <- ind
  } else{
    l <- M
  }
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}

get_nn_matrix <- function(coords, m){
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
  coords_ord <- coords[nn_data$neighbor.info$ord,]
  nn_idx <- sapply(1:(n - 1), get_NN_ind, nn_data$neighbor.info$n.indx[-1], m) %>%
    t()
  neighbor_dist <- sapply(
    2:n, i_dist, nn_data$neighbor.info$n.indx[-1], coords_ord)
  distance_matrix <- get_dist_matrix(neighbor_dist, nn_data$neighbor.info$n.indx[-1])
  
  return(list(
    ord = ord, coords_ord = coords_ord, nn_idx = nn_idx, 
    distance_matrix = distance_matrix
  ))
}

m <- 3
nn_data <- get_nn_matrix(coords, m)
nn_data
