# Comparison of GP and NNGP likelihood for Poisson response
#
# The basic idea is that the NNGP likelihood should be a cheaper approximation
# of the GP likelihood.
# In order to assess, we write the joint distribution of NNGP as the product of 
# conditionals, that is
#
# p(y_1) * p(y_2 | N(y_2)) * p(y_3 | N(y_3)) * p(y_n | N(y_n)),
#
# where N(y_k) is the neighbor set of y_k.

library(dplyr)
library(purrr)
library(ggplot2)
library(mvtnorm)
source("./nngp-notes/codes/stat_utils.R")
source("./nngp-notes/codes/poisson_functs.R")

n <- 120

set.seed(126)
coords <- cbind(runif(n), runif(n))
# coords <- cbind(rep(1:10, each = 10), rep(1:10, 10))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3
mu_w <- 0
# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, mean = rep(mu_w, n), sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
# C <- sigma^2 * exp(-phi^2 * d)
# w <- c(t(chol(C)) %*% rnorm(n))
y <- rpois(n, 2 + exp(w))

hist(y)

table(y == 0) %>%
  prop.table()

# GP
jesus <- function() {
  n_it <- 100
  p_init <- c(log(1), log(1))
  w_init <- log(y + .5)
  for (i in 1:n_it) {
    print(i)
    sigma <- exp(p_init[1])
    phi <- exp(p_init[2])
    sigma_w <- sigma^2 * exp(-phi^2 * d)
    
    # Gaussian proxy
    gauss <- estim_gaussian_proxy(w_init, y, sigma_w)
    w_hat <- gauss$mu
    sigma_hat <- gauss$sigma
    
    # Laplace
    fitted_parms <- optim(
      p_init, laplace, y = y, coords = coords, w_hat = w_hat, 
      sigma_hat = sigma_hat, a = 2, mu_w = rep(mu_w, n), method = "BFGS", hessian = T, 
      control = list(fnscale = -1))
    if (fitted_parms$convergence != 0) stop("Convergence error Laplace")
    
    if (max(abs(p_init - fitted_parms$par)) < 1e-4) break
    if (i == n_it) stop("Max iteraction number reached")
    p_init <- fitted_parms$par
    w_init <- w_hat
    print(exp(fitted_parms$par))
  }
  return(fitted_parms)
}

# debugonce(jesus)
fitted_parms <- jesus()

sigma_hat <- deltamethod(
  list(~exp(x1), ~exp(x2)), fitted_parms$par, -solve(fitted_parms$hessian), ses = F)

ci_parms <- confint2(exp(fitted_parms$par), sigma_hat) %>%
  mutate(
    par = c("sigma", "phi"),
    true = c(sigma_true, phi_true)
  )

ggplot(data = ci_parms, aes(x = "")) +
  facet_wrap(~par, scales = "free") +
  geom_linerange(aes(ymin = lb, ymax = ub)) +
  geom_point(aes(y = mu)) +
  geom_hline(aes(yintercept = true), lty = 2, col = "red") +
  labs(y = "Estimate", x = "")

