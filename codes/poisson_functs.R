library(dplyr)
library(tidyr)
library(purrr)
library(mvtnorm)
library(msm)

# -----------------------------------------------------------------------------
# Poisson response
# -----------------------------------------------------------------------------

# f(w | y, theta) \propto f(y | w, theta) * f(w | theta)
posterior <- function(y, w, sigma_w, mu = rep(0, length(w)), log = F) {
  ll <- sum(dpois(y, exp(w), log = T)) + dmvnorm(w, mu, sigma_w, log = T)
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

# Gaussian proxy of f(y | w) using Taylor expansion
proxy_fyw0 <- function(y, w0){
  b <- gradient(w0, y) - hessian(w0) * w0
  c <- -hessian(w0)
  c <- ifelse(c < 0, sqrt(.Machine$double.xmin), c)
  mu <- b / c

  parms <- list(b = b, c = c, mu = mu)
  return(parms)
}

# Gaussian proxy of f(w | y) using Taylor expansion
# proxy_fwy0 <- function(mu_w, sigma_w, b_y, c_y){
#   prec_m <- solve(sigma_w)
#   b <- prec_m %*% mu_w + b_y
#   c <- prec_m + diag(c_y, length(mu_w))
# 
#   mu <- solve(c, b) %>%
#     as.vector()
#   parms <- list(mu = mu, sigma = solve(c))
#   return(parms)
# }

# MODIFY TO ANY MU_W
# proxy_fwy02 <- function(sigma_w, y, w0){
#   hess_m <- diag(hessian(w0), length(w0))
#   sigma <- solve(solve(sigma_w) - hess_m)
#   mu <- sigma %*% (gradient(w0, y) - hess_m %*% w0)
#   parms <- list(mu = as.vector(mu), sigma = sigma)
#   return(parms)
# }

estim_gaussian_proxy <- function(
    w0, y, sigma_w, mu_w = rep(0, length(y)), control = list(it = 100, tol = 1e-4)) {
  for (i in 1:control$it) {
    # Taylor expansion of f(y | w) to approximate to a normal distribution
    # fit_fyw0 <- proxy_fyw0(y, w0)

    # Calculating the parameters. Note that f(w | y) is already Gaussian
    # fit_fwy0 <- proxy_fwy0(mu_w, sigma_w, fit_fyw0$b, fit_fyw0$c)

    hess_m <- diag(hessian(w0), length(w0))
    sigma <- solve(solve(sigma_w) - hess_m)
    mu <- as.vector(sigma %*% (gradient(w0, y) - hess_m %*% w0))

    # if (max(abs(w0 - fit_fwy0$mu)) < control$tol) break
    if (max(abs(w0 - mu)) < control$tol) break
    if (i == control$it) stop("Max iteraction number reached")
    # w0 <- fit_fwy0$mu
    w0 <- mu
  }

  # parms <- list(
    # mu = fit_fwy0$mu, sigma = fit_fwy0$sigma)
  parms <- list(
    mu = mu, sigma = sigma)
  return(parms)
}

# Laplace approximation
# f(theta | y) \approx f(y | w_hat, theta) * f(w_hat | theta) / MVN(w_hat, sigma)
laplace <- function(y, coords, p_init, w_hat, sigma_hat, mu_w = rep(0, length(y))) {
  sigma <- exp(p_init[1])
  phi <- exp(p_init[2])
  d <- dist(coords) %>%
    as.matrix()
  sigma_w <- sigma^2 * exp(-phi^2 * d)
  
  # Covariance matrix of Gaussian proxy
  hess_m <- diag(hessian(w_hat), length(w_hat))
  sigma_proxy <- solve(solve(sigma_w) - hess_m)

  denom <- -0.5 * log(det(sigma_proxy) + sqrt(.Machine$double.xmin))
  
  ll <- posterior(y, w_hat, sigma_w, mu_w, log = T) - denom

  return(ll)
}
