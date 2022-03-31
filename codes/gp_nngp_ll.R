library(dplyr)
library(spNNGP)
library(mvtnorm)

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

# Return the parameters of conditional distribution p(y_star | y) from a MVN
cond_mvn <- function(y_new, y, p, coords){
  d <- dist(coords) %>%
    as.matrix()
  d_star <- rbind(y_new, coords) %>%
    dist() %>%
    as.matrix()
  d_star <- d_star[-1, 1]
  sigma <- p[1]
  phi <- p[2]
  tau <- p[3]
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(length(y))
  sigma_star <- sigma^2 * exp(-phi^2 * d_star)
  sigma_star_star <- sigma^2
  ynew_mean <- sigma_star %*% solve(Sigma, y)
  ynew_sigma <- sigma_star_star - sigma_star %*% solve(Sigma, sigma_star)
  y_param <- c(ynew_mean, ynew_sigma)
  return(y_param)
}

# GP likelihood
gp_ll = function(p, coords, y) {
  sigma <- p[1] %>%
    as.numeric()
  phi <- p[2] %>%
    as.numeric()
  tau <- p[3] %>%
    as.numeric()
  d <- dist(coords) %>%
    as.matrix()
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(nrow(coords))
  mu <- rep(0, nrow(coords))
  ll <- dmvnorm(y, mu, Sigma, log = TRUE)
  return(-ll)
}

# NNGP likelihood
nngp_ll = function(p, coords, y, m) {
  sigma <- p[1] %>%
    as.numeric()
  phi <- p[2] %>%
    as.numeric()
  tau <- p[3] %>%
    as.numeric()
  nn_data <- get_nn(coords, m)
  cond_dens <- dnorm(y[1], 0, sqrt(sigma^2 + tau^2))
  # Get the conditional distributions for NNGP
  for (i in 1:nrow(nn_data$nn_idx)) {
    nn <- nn_data$ord[nn_data$nn_idx[i,]]
    parms <- cond_mvn(
      coords[i + 1,], y[nn], c(sigma, phi, tau), matrix(coords[nn,], nrow = length(nn)))
    cond_dens[i + 1] <- dnorm(y[i + 1], parms[1], sqrt(parms[2]))
  }
  ll <- prod(cond_dens) %>%
    log()
  return(-ll)
}