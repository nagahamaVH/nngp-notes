library(dplyr)
library(spNNGP)
library(mvtnorm)

# Rprof()

# -----------------------------------------------------------------------------
# Functions to find the neighrest neighbors
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Functions to prediction of MVN
# -----------------------------------------------------------------------------
# Return the conditional mean and variance from a MVN, i.e the parameters of
# p(idx | given_idx)
cond_mvn <- function(idx, given_idx, y, sigma){
  Sigma <- sigma[given_idx, given_idx]
  sigma_star <- sigma[idx, given_idx]
  sigma_star_star <- sigma[idx, idx]
  ynew_mean <- sigma_star %*% solve(Sigma, y)
  ynew_sigma <- sigma_star_star - sigma_star %*% solve(Sigma, sigma_star)
  y_param <- c(ynew_mean, ynew_sigma)
  return(y_param)
}

# Predict a new observation given pair of coordinates
pred_cond_mvn <- function(x_star, x, y, p) {
  d <- dist(x) %>%
    as.matrix()
  d_star <- rbind(x_star, x) %>%
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

# -----------------------------------------------------------------------------
# Gaussian response
# -----------------------------------------------------------------------------
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
  d <- dist(coords) %>%
    as.matrix()
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(length(y))
  nn_data <- get_nn(coords, m)
  cond_dens <- dnorm(y[1], 0, sqrt(Sigma[1, 1]), log = T)
  # Get the conditional distributions for NNGP
  for (i in 1:nrow(nn_data$nn_idx)) {
    nn <- nn_data$ord[nn_data$nn_idx[i,]]
    parms <- cond_mvn(idx = i + 1, given_idx = nn, y = y[nn], Sigma)
    cond_dens[i + 1] <- dnorm(y[i + 1], parms[1], sqrt(parms[2]), log = T)
  }
  ll <- sum(cond_dens)
  return(-ll)
}

# -----------------------------------------------------------------------------
# Poisson response
# -----------------------------------------------------------------------------
# Joint distribution
joint_pois <- function(y, w, coords, p){
  sigma <- p[1] %>%
    as.numeric()
  phi <- p[2] %>%
    as.numeric()
  d <- dist(coords) %>%
    as.matrix()
  Sigma <- sigma^2 * exp(-phi^2 * d)
  dens <- dpois(y, exp(w)) * dmvnorm(w, sigma = Sigma)
  return(dens)
}

# define before proxy marginal
sigma <- p[1] %>%
  as.numeric()
phi <- p[2] %>%
  as.numeric()
d <- dist(coords) %>%
  as.matrix()
Sigma <- sigma^2 * exp(-phi^2 * d)
Sigma_inv <- solve(Sigma)

# Proxy of marginal using laplace approx.
marginal_proxy <- function(w0, y, Sigma, Sigma_inv){
  n <- length(y)

  ll <- expression(y * w0 - exp(w0) - 1 / 2 * (w0 * Sigma_inv * w0))
  hessian_expr <- D(ll, "w0") %>%
    D(., "w0")
  h <- eval(hessian_expr)
  
  joint <- dpois(y, exp(w0)) * dmvnorm(w0, sigma = Sigma)
  proxy <- joint * (2 * pi)^(n / 2) * h
  return(proxy)
}

gp <- optim(
  init_parms, marginal_proxy, y = y, Sigma = Sigma, Sigma_inv = Sigma_inv, 
  method = "L-BFGS-B")


# dx2x <- deriv(~ x^2, "x") ; dx2x
# mode(dx2x)
# 
# trig.exp <- expression(sin(cos(x + y^2)))
# ( D.sc <- D(trig.exp, "x") )
# all.equal(D(trig.exp[[1]], "x"), D.sc)
