library(dplyr)
library(tidyr)
library(purrr)
library(mvtnorm)
library(msm)

# -----------------------------------------------------------------------------
# Poisson response
# -----------------------------------------------------------------------------

# f(w | y, theta) \propto f(y | w, theta) * f(w | theta)
posterior <- function(y, w, mu = rep(0, length(w)), sigma, log = F) {
  ll <- sum(dpois(y, exp(w), log = T)) + dmvnorm(w, mu, sigma, log = T)
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

# # Checking for n = 2
# w0 <- rep(0, n)
# for (i in 1:100) {
#   proxy <- proxy_fyw0(y, w0)
#   e <- abs(w0 - proxy$mu) %>%
#     sum() / length(w0)
#   if (e < 1e-4) {
#     print("Converged!")
#     break
#   }
#   w0 <- proxy$mu
# }
# 
# # Exact
# w_grid <- sapply(w, function(w) seq(w - w/2, w + w/2, length.out = 50)) %>%
#   as_tibble() %>%
#   expand.grid()
# dens <- apply(w_grid, MARGIN = 1, function(x){
#   dpois(y, exp(x), log = T) %>%
#     sum()
# })
# df <- bind_cols(
#   w_grid,
#   d = dens
# )
# 
# df <- bind_cols(
#   w_grid,
#   exact = apply(w_grid, MARGIN = 1, function(x){
#     dpois(y, as.numeric(exp(x)), log = T) %>%
#       sum()
#   }),
#   proxy = apply(w_grid, MARGIN = 1, dmvnorm, mean = proxy$mu, sigma = diag(1 / proxy$c), log = T)
# )
# 
# # Visualising the approximation
# library(plotly)
# plot_ly() %>%
#   add_trace(data = df, x = ~V1, y = ~V2, z = ~exact, type = "mesh3d") %>%
#   add_trace(data = df, x = ~V1, y = ~V2, z = ~proxy, type = "mesh3d")

# Gaussian proxy of f(w | y) using Taylor expansion
proxy_fwy0 <- function(mu_w, sigma_w, b_y, c_y){
  prec_m <- solve(sigma_w)
  b <- prec_m %*% mu_w + b_y
  c <- prec_m + diag(c_y, length(mu_w))

  mu <- solve(c, b) %>%
    as.vector()
  parms <- list(mu = mu, sigma = solve(c))
  return(parms)
}

estim_gaussian_proxy <- function(
    w0, y, sigma_w, mu_w = rep(0, length(y)), control = list(it = 100, tol = 1e-4)) {
  
  # Taylor expansion of f(y | w) to approximate to a normal distribution
  fit_fyw0 <- proxy_fyw0(y, w0)
  
  # Calculating the parameters. Note that f(w | y) is already Gaussian
  fit_fwy0 <- proxy_fwy0(mu_w, sigma_w, fit_fyw0$b, fit_fyw0$c)
  
  parms <- list(
    mu = fit_fwy0$mu, sigma = fit_fwy0$sigma)
  
  return(parms)
}

# Laplace approximation
# f(theta | y) \approx f(y | w_hat, theta) * f(w_hat | theta) / MVN(w_hat, sigma)
laplace <- function(y, coords, p_init, w_hat) {
  sigma <- exp(p_init[1])
  phi <- exp(p_init[2])
  mu <- rep(0, length(y))
  d <- dist(coords) %>%
    as.matrix()
  sigma_w <- sigma^2 * exp(-phi^2 * d)
  
  # Transform singular matrix to non-singular. What about SVD?
  # while (det(sigma_w) == 0) {
  #   sigma_w <- sigma_w + diag(1e-3, nrow(d))
  # }
    
  # Covariance matrix of Gaussian proxy
  prec_proxy <- solve(sigma_w) + diag(-hessian(w_hat), length(w_hat))
  sigma_proxy <- solve(prec_proxy)

  denom <- -0.5 * log(det(sigma_proxy) + sqrt(.Machine$double.xmin))
  ll <- posterior(y, w_hat, mu, sigma_w, log = T) - denom

  return(ll)
}
