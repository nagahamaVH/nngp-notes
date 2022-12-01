# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(rstan)
library(bayesplot)

# ------------------- Setup ---------------------------------------------------
options(mc.cores = parallel::detectCores())

model_name <- "poisson_gp"

data_board <- pins::board_folder("./data", versioned = T)
model_board <- pins::board_folder("models", versioned = T)
# -----------------------------------------------------------------------------

n <- 200
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

stan_data <- list(
  n = n,
  y = y,
  Z = Z,
  p = dim(Z)[2],
  D = D)

# ------------------------ Stan parameters ------------------------------------
n_chain <- 3
n_it <- 600
model_file <- "./codes/poisson_gp.stan"
# -----------------------------------------------------------------------------

stan_fit <- stan(
  file = model_file,
  data = stan_data,
  chains = n_chain,
  iter = n_it,
  seed = 171)

stan_fit |>
  mcmc_trace(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens_overlay(pars = c("sigma", "l"), regex_pars = "beta")

stan_fit |>
  mcmc_dens(pars = c("sigma", "l"), regex_pars = "beta")

summary(stan_fit, pars = c("beta", "sigma", "l"))$summary

# Save data
data_meta <- list(
  cov_kernel = "se",
  n = n,
  max_s = max_s,
  cov_par = list(
    sigma = sigma,
    l = l
  ),
  beta = beta
)

pins::pin_write(
  data_board, stan_data, model_name, type = "rds", metadata = data_meta)

# Save model
data_version <- data_board |>
  pins::pin_versions(model_name) |>
  tail(1) %>%
  dplyr::pull(version)

model_meta <- list(
  model = model_name,
  file = model_file,
  data = data_version,
  n_chain = n_chain,
  n_it = n_it
)

pins::pin_write(
  model_board, stan_fit, model_name, type = "rds", metadata = model_meta)
