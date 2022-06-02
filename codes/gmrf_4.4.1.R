library(dplyr)
library(ggplot2)
library(cowplot)

# f(b | y, theta) \propto f(y | b, theta) * f(b | theta)
posterior <- function(b, y, theta, log = F) {
  mu <- theta[1]
  sigma <- sqrt(theta[2])
  ll <- dpois(y, exp(b), log = T) + dnorm(b, mu, sigma, log = T)
  if (!log) {
    ll <- exp(ll)
  }
  return(ll)
}

gradient <- function(b, y, sigma_inv) {
  y - exp(b) - sigma_inv * b
}

hessian <- function(b, sigma_inv) {
  -(exp(b) + sigma_inv) 
}

# Return the mean and variance of Gaussian approximation for a given value b
gaussian_proxy <- function(b0, y, sigma_inv) {
  b <- gradient(b0, y, sigma_inv) - hessian(b0, sigma_inv) * b0
  c <- -hessian(b0, sigma_inv)
  mu <- b / c
  sigma_sq <- 1 / c
  parms <- list(mu = as.vector(mu), sigma_sq = sigma_sq)
  return(parms)
}

# Fisher scoring (Iteratively reweighted least squares) to match the mode of
# posterior distribution with the Gaussian approximation
estim_gaussian_proxy <- function(b0, y, sigma_inv, control = list(it = 100, tol = 1e-4)) {
  for (i in 1:control$it) {
    proxy <- gaussian_proxy(b0, y, sigma_inv)
    mu <- proxy$mu
    if (abs(b0 - mu) < control$tol) break
    b0 <- mu
  }
  parms <- list(
    mu = as.vector(mu), sigma_sq = proxy$sigma_sq, n_iter = i,
    convergence = ifelse(i < control$it, T, F))
  return(parms)
}

# Simulating data
# y | b ~ Poisson(mu) where b = log mu
# b ~ N(0, 1/0.001)
n <- 1
theta <- c(0, 1/0.001)
y <- 3

fit_gaussian <- estim_gaussian_proxy(0, y, solve(theta[2]))
mu_proxy <- fit_gaussian$mu
sd_proxy <- sqrt(fit_gaussian$sigma_sq) %>%
  as.vector()

# Step size
h <- 0.01

b_grid <- seq(mu_proxy - 5 * sd_proxy, mu_proxy + 5 * sd_proxy, by = h)
unnorm_exact <- posterior(b_grid, y = y, theta = theta)

b0_grid <- c(0, 0.5, 1, 1.5)

plots <- list()
for (i in 1:length(b0_grid)) {
  proxy_parms <- gaussian_proxy(b0_grid[i], y, solve(theta[2]))
  
  df <- tibble(
    b = b_grid,
    exact = unnorm_exact / sum(unnorm_exact * h),
    proxy = dnorm(b_grid, proxy_parms$mu, sqrt(proxy_parms$sigma_sq))
  )

  p <- ggplot(df, aes(x = b, y = exact)) +
    geom_line() +
    geom_line(aes(y = proxy), linetype = 2) +
    geom_point(aes_string(x = b0_grid[i], y = 0)) +
    labs(y = "Density") +
    xlim(c(-2, 3)) +
    theme_light()
  plots[[i]] <- p
}

plot_grid(plotlist = plots, nrow = 2)
