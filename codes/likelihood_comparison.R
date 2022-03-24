# Comparison of GP and NNGP likelihood
#
# The basic idea is that the NNGP likelihood should be a cheaper approximation
# of the GP likelihood.
# 
# In order to assess, we write the joint distribution of NNGP as the product of 
# conditionals, that is
#
# p(y_1) * p(y_2 | N(y_2)) * p(y_3 | N(y_3)) * p(y_n | N(y_n)),
#
# where N(y_k) is the neighbor set of y_k.

library(dplyr)
library(mvtnorm)
source("./nngp-notes/codes/nn_matrix.R")

n <- 5
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

ord <- order(coords[,1])
coords <- coords[ord,]

sigma_sq <- 2
phi <- 3
tau_sq <- 0.9

d <- dist(coords) %>%
  as.matrix()

w <- rmvnorm(1, rep(0, n), sigma_sq * exp(-phi * d)) %>%
  c()
eps <- rnorm(n, sqrt(tau_sq))
y <- w + eps

# Predicting from a GP
gp_pred <- function(y_new, y, sigma, coords){
  d2 <- rbind(y_new, coords) %>%
    dist() %>%
    as.matrix()
  sigma_star <- sigma_sq * exp(-phi * d2[-1, 1])
  sigma_star_star <- sigma_sq
  ynew_mean <- sigma_star %*% solve(sigma, y)
  ynew_sigma <- sigma_star_star - sigma_star %*% solve(sigma, sigma_star)
  y_param <- c(ynew_mean, ynew_sigma)
  return(y_param)
}

sigma <- sigma_sq * exp(-phi * d) + tau_sq * diag(n)

# ------------------------------------------
# Full GP
# -----------------------------------------
sapply(2:n, function(i) gp_pred(coords[i,], y[1:i], sigma[1:i, 1:i], coords[1:i,])) %>%
  t()

# ------------------------------------------
# NNGP
# -----------------------------------------
nn_data <- get_nn(coords, m)
for (i in 1:nrow(nn_data$nn_idx)) {
  nn <- nn_data$ord[nn_data$nn_idx[i,]]
  print(gp_pred(coords[i + 1,], y[nn], sigma[nn, nn], coords[nn,]))
}

# GP likelihood
gp_ll = function(p, x, y) {
  sigma_sq = exp(p[1])
  phi = exp(p[2])
  tau_sq = exp(p[3])
  d <- dist(x) %>%
    as.matrix()
  sigma <- sigma_sq * exp(-phi * d) + tau_sq * diag(nrow(x))
  mu <- rep(0, nrow(x))
  ll <- dmvnorm(y, mu, sigma, log = TRUE)
  return(-ll)
}

# NNGP likelihood
nngp_ll = function(p, x, y, m) {
  sigma_sq = exp(p[1])
  phi = exp(p[2])
  tau_sq = exp(p[3])
  d <- dist(x) %>%
    as.matrix()
  sigma <- sigma_sq * exp(-phi * d) + tau_sq * diag(nrow(x))
  nn_data <- get_nn(coords, m)
  cond_dens <- NULL
  cond_dens[1] <- dnorm(y[1], 0, sqrt(sigma[1, 1]))
  for (i in 1:nrow(nn_data$nn_idx)) {
    nn <- nn_data$ord[nn_data$nn_idx[i,]]
    parms <- gp_pred(coords[i + 1,], y[nn], sigma[nn, nn], coords[nn,])
    cond_dens[i + 1] <- dnorm(y[i + 1], parms[1], sqrt(parms[2]))
  }
  ll <- prod(cond_dens) %>%
    log()
  return(-ll)
}

parms <- log(c(2, 3, 0.9))
gp_ll(parms, coords, y)
nngp_ll(parms, coords, y, m)
