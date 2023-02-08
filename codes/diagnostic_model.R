# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(ggplot2)
library(posterior)
library(bayesplot)

# ------------------- Setup ---------------------------------------------------
options(mc.cores = parallel::detectCores())

model_name <- "poisson_nngp_zhang"

data_board <- pins::board_folder("./data", versioned = T)
model_board <- pins::board_folder("models", versioned = T)
# -----------------------------------------------------------------------------
pins::pin_versions(model_board, model_name)$version[1]

meta <- pins::pin_meta(model_board, pins::pin_versions(model_board, model_name)$version[1])

stan_fit <- pins::pin_read(model_board, model_name, version = pins::pin_versions(model_board, model_name)$version[1])
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
  mcmc_dens(pars = c("l", "sigma", "w[1]", "w[2]"), regex_pars = "beta")

summary(stan_fit, pars = c("beta", "sigma", "l"))$summary

stan_fit |>
  mcmc_intervals(
    regex_pars = "w",
    prob = 0.8,
    prob_outer = 0.99,
    point_est = "mean",
    inner_size = 2,
    point_size = 2
  )

post <- summarise_draws(stan_fit, default_summary_measures())
post <- post |>
  dplyr::filter(
    stringr::str_detect(post$variable, "^w\\[")
  ) |>
  dplyr::mutate(
    w = x
  )

ggplot(post, aes(x = variable, y = median)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  geom_point() +
  geom_point(aes(x = variable, y = w), shape = "x", colour = "red", size = 3) +
  coord_flip() +
  theme_minimal()
