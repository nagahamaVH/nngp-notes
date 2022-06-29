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
library(ggplot2)
library(mvtnorm)
source("./nngp-notes/codes/stat_utils.R")
source("./nngp-notes/codes/poisson_functs.R")

n <- 100

set.seed(126)
coords <- cbind(runif(n), runif(n))
# coords <- cbind(rep(1:10, each = 10), rep(1:10, 10))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3

# Generate data
d <- dist(coords) %>%
  as.matrix()
sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
w <- rmvnorm(1, mean = rep(0, n), sigma = sigma_w) %>%
  c()
# C <- sigma^2 * exp(-phi^2 * d)
# w <- c(t(chol(C)) %*% rnorm(n))
y <- rpois(n, exp(w))

hist(y)

table(y == 0) %>%
  prop.table()

# Decay
plot(d[,1], sigma_w[,1])

# Laplace approximation
p_init <- c(log(5), log(5))
w_init <- log(y + .5)
mu_w <- 0
fitted_parms <- optim(
  p_init, laplace, y = y, coords = coords, w_init = w_init, mu_w = rep(mu_w, n), 
  method = "Nelder-Mead", hessian = T, control = list(fnscale = -1))

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
