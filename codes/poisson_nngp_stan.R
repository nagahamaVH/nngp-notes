# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())

n <- 500
max_s <- 3

set.seed(163)
coords <- cbind(runif(n, max = max_s), runif(n, max = max_s))

ord <- order(coords[,1])
coords <- coords[ord,]

sigma <- 2.6
l <- .4
beta <- c(-.5, 1.8)

D <- dist(coords) |>
  as.matrix()
C <- sigma^2 * exp(-D / (2 * l^2))
x <- mvtnorm::rmvnorm(1, sigma = C, method = "chol", checkSymmetry = F) |>
  c()
Z <- cbind(1, rnorm(n, 2))
eta <- exp(tcrossprod(Z, t(beta)) + x)
y <- rpois(n, eta)

m <- 10

get_NN_ind <- function(ind, ind_distM_i, M){
  l <- ifelse(ind < M, ind, M)
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}

get_nn <- function(coords, m){
  n <- dim(coords)[1]
  nn_data <- spNNGP::spConjNNGP(
    rep(0, n) ~ 1, coords = coords,
    n.neighbors = m,
    theta.alpha = c("phi" = 5, "alpha" = 0.5),
    sigma.sq.IG = c(2, 1),
    cov.model = "exponential",
    return.neighbor.info = T, fit.rep = F, 
    verbose = F)
  ord <- nn_data$neighbor.info$ord
  nn_idx <- sapply(1:(n - 1), get_NN_ind, nn_data$neighbor.info$n.indx[-1], m) |>
    t()
  
  return(list(ord = ord, nn_idx = nn_idx))
}

nn <- get_nn(coords, m)

stan_data <- list(
  n = n,
  y = y,
  Z = Z,
  p = dim(Z)[2],
  D = D,
  nn = nn$nn_idx,
  m = m)

hist(y)

# Stan model
model <- "./codes/poisson_nngp.stan"

stan_fit <- stan(
  file = model,
  data = stan_data,
  chains = 3,
  iter = 800,
  seed = 171)

stan_fit |>
  mcmc_trace(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens_overlay(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens(pars = c("sigma", "l"), regex_pars = "beta")

summary(stan_fit, pars = c("beta", "sigma", "l"))$summary

