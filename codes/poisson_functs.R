library(dplyr)
library(tidyr)
library(mvtnorm)
library(msm)

# -----------------------------------------------------------------------------
# Poisson response
# -----------------------------------------------------------------------------

# f(w | y, theta) \propto f(y | w, theta) * f(w | theta)
posterior <- function(y, w, sigma_w, a = 0, mu = rep(0, length(w)), log = F) {
  ll <- sum(dpois(y, a + exp(w), log = T)) + dmvnorm(w, mu, sigma_w, log = T)
  if (!log) {
    ll <- exp(ll)
  }
  return(ll)
}

gradient <- function(w, y) {
  y - exp(w)
}

hessian <- function(w) {
  -exp(w)
}

estim_gaussian_proxy <- function(
    w0, y, sigma_w, control = list(it = 100, tol = 1e-4)) {
  for (i in 1:control$it) {
    hess_m <- diag(hessian(w0), length(w0))
    sigma <- chol2inv(chol(chol2inv(chol(sigma_w)) - hess_m))
    mu <- as.vector(sigma %*% (gradient(w0, y) - hess_m %*% w0))
    if (max(abs(w0 - mu)) < control$tol) break
    if (i == control$it) stop("Max iteraction number reached")
    w0 <- mu
  }

  parms <- list(
    mu = mu, sigma = sigma)
  return(parms)
}

# Laplace approximation
# f(theta | y) \approx f(y | w_hat, theta) * f(w_hat | theta) / MVN(w_hat, sigma)
laplace <- function(y, coords, p_init, w_init, a = 0, mu_w = rep(0, length(y))) {
  sigma <- exp(p_init[1])
  phi <- exp(p_init[2])
  d <- dist(coords) %>%
    as.matrix()
  sigma_w <- sigma^2 * exp(-phi^2 * d)
  print(exp(p_init))
  # Gaussian proxy
  gauss <- estim_gaussian_proxy(w_init, y, sigma_w)
  w_hat <- gauss$mu
  sigma_hat <- gauss$sigma

  denom <- -(length(y) / 2) * log(2 * pi) - 0.5 * log(det(sigma_hat) + sqrt(.Machine$double.xmin))
  ll <- posterior(y, w_hat, sigma_w, a, mu_w, log = T) - denom
  return(ll)
}
