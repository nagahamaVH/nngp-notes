library(dplyr)
library(mvtnorm)
library(rstan)
library(bayesplot)

n <- 120
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1.7
phi_true <- 3.2

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
y <- rpois(n, exp(w))

# ---------------------------------------------------------------
# STAN
# ---------------------------------------------------------------
model <- '
functions{
  matrix GP(matrix x, real sigma, real phi) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sigma^2;
      for (j in (i + 1):N) {
        K[i, j] = sigma^2 * exp(-phi^2 * x[i,j]);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sigma^2;
    return K;
  }
}
data{
  int<lower = 1> N;
  int<lower = 0> y[N];
  matrix[N, N] DMat; // Distance matrix
}
parameters{
  vector[N] k;
  
  // GP standard deviation parameters
  real<lower=0> sigma;
  // GP length-scale parameters
  real<lower=0> phi;
}
model{
  matrix[N,N] SIGMA;
  vector[N] mu;
  SIGMA = GP(DMat, sigma, phi);
  k ~ multi_normal(rep_vector(0,N), SIGMA);
  
  // Priors
  phi ~ uniform(0, 10);
  sigma ~ uniform(0, 10);
  
  y ~ poisson_log(k);
}
'

stan_data <- list(
  N = n,
  y = y,
  DMat = d)

# Compute the distance matrix
stan_fit <- stan(
  model_code = model,
  data = stan_data,
  chains = 3,
  iter = 5000,
  seed = 171)

stan_fit %>%
  mcmc_trace(pars = c("sigma", "phi"))

stan_fit %>%
  mcmc_dens(pars = c("sigma", "phi"))

print(stan_fit, pars = c("sigma", "phi"))
