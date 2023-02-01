# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(rstan)
library(bayesplot)
source("./codes/NNMatrix.R")

# ------------------- Setup ---------------------------------------------------
options(mc.cores = parallel::detectCores())

model_name <- "poisson_nngp_zhang"

data_board <- pins::board_folder("./data", versioned = T)
model_board <- pins::board_folder("models", versioned = T)
# -----------------------------------------------------------------------------

stan_fit <- pins::pin_read(model_board, model_name)
data <- pins::pin_read(data_board, model_name)
# -----------------------------------------------------------------------------

np <- nuts_params(stan_fit)

stan_fit |>
  mcmc_trace(
    pars = c("l", "sigma"), regex_pars = "beta", np = np)

stan_fit |>
  mcmc_dens_overlay(pars = c("l", "sigma"), regex_pars = "beta")

# Should be R < 1.05 
stan_fit |>
  rhat() |>
  mcmc_rhat()

# Desired ESS >= 0.5
stan_fit |>
  neff_ratio() |>
  mcmc_neff()

mcmc_pairs(
  stan_fit,
  regex_pars = "beta",
  pars = c("l", "sigma"),
  transformations = list(l = "log", sigma = "log"),
  off_diag_fun = "hex",
  np = np
)

stan_fit |>
  mcmc_dens(pars = c("l", "sigma"), regex_pars = "beta")

stan_fit |>
  mcmc_dens(pars = c("beta[1]", "w[1]", "w[2]", "w[3]"))

summary(stan_fit, pars = c("beta", "sigma", "l"))$summary
