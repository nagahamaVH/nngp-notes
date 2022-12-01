# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(rstan)
library(bayesplot)
source("./codes/NNMatrix.R")

options(mc.cores = parallel::detectCores())

n <- 80
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

m <- 3

nn <- NNMatrix(coords, m)

stan_data <- list(
  n = n,
  y = y,
  Z = Z,
  p = dim(Z)[2],
  nn = nn$NN_ind,
  D = nn$NN_dist,
  D_m = nn$NN_distM,
  m = m)

hist(y)

# Stan model
model <- "./codes/poisson_nngp.stan"

stan_fit <- stan(
  file = model,
  data = stan_data,
  chains = 3,
  iter = 300,
  seed = 171)

stan_fit |>
  mcmc_trace(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens_overlay(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens(pars = c("sigma", "l"), regex_pars = "beta")

summary(stan_fit, pars = c("beta", "sigma", "l"))$summary

