# Comparison of GP and NNGP likelihood
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
source("./nngp-notes/codes/nn_matrix.R")

n <- 5
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

sigma <- 3
phi <- 4
tau <- 0.6

d <- dist(coords) %>%
  as.matrix()

w <- rmvnorm(1, sigma = sigma^2 * exp(-phi^2 * d)) %>%
  c()
eps <- rnorm(n, 0, tau)
y <- w + eps

# Conditional distribution p(y_star | y) from a MVN
gp_pred <- function(y_new, y, p, coords){
  d <- dist(coords) %>%
    as.matrix()
  d_star <- rbind(y_new, coords) %>%
    dist() %>%
    as.matrix()
  d_star <- d_star[-1, 1]
  sigma <- p[1]
  phi <- p[2]
  tau <- p[3]
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(length(y))
  sigma_star <- sigma^2 * exp(-phi^2 * d_star)
  sigma_star_star <- sigma^2
  ynew_mean <- sigma_star %*% solve(Sigma, y)
  ynew_sigma <- sigma_star_star - sigma_star %*% solve(Sigma, sigma_star)
  y_param <- c(ynew_mean, ynew_sigma)
  return(y_param)
}

# GP likelihood
gp_ll = function(p, x, y) {
  sigma <- p[1] %>%
    as.numeric()
  phi <- p[2] %>%
    as.numeric()
  tau <- p[3] %>%
    as.numeric()
  d <- dist(x) %>%
    as.matrix()
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(nrow(x))
  mu <- rep(0, nrow(x))
  ll <- dmvnorm(y, mu, Sigma, log = TRUE)
  return(-ll)
}

# NNGP likelihood
nngp_ll = function(p, x, y, m) {
  sigma <- p[1] %>%
    as.numeric()
  phi <- p[2] %>%
    as.numeric()
  tau <- p[3] %>%
    as.numeric()
  nn_data <- get_nn(x, m)
  cond_dens <- dnorm(y[1], 0, sqrt(sigma^2 + tau^2))
  # Get the conditional distributions for NNGP
  for (i in 1:nrow(nn_data$nn_idx)) {
    nn <- nn_data$ord[nn_data$nn_idx[i,]]
    parms <- gp_pred(
      x[i + 1,], y[nn], c(sigma, phi, tau), matrix(x[nn,], nrow = length(nn)))
    cond_dens[i + 1] <- dnorm(y[i + 1], parms[1], sqrt(parms[2]))
  }
  ll <- prod(cond_dens) %>%
    log()
  return(-ll)
}

parms_grid <- seq(0.2, sigma + 2, length.out = 20) %>%
  tibble(sigma = ., phi, tau)

gp_uni <- apply(parms_grid, MARGIN = 1, FUN = gp_ll, coords, y)
nngp_uni <- apply(parms_grid, MARGIN = 1, FUN = nngp_ll, coords, y, m)

parms_grid <- parms_grid %>%
  bind_cols("gp" = gp_uni, "nngp" = nngp_uni)

ggplot(parms_grid, aes(x = sigma, y = gp, colour = "GP")) +
  geom_line() +
  geom_line(aes(y = nngp, colour = "NNGP")) +
  geom_point(data = parms_grid[which.min(parms_grid$gp),]) +
  geom_point(
    data = parms_grid[which.min(parms_grid$nngp),],
    aes(y = nngp, colour = "NNGP")) +
  labs(x = expression(sigma), y = "log-likelihood", colour = "") +
  theme_light()

