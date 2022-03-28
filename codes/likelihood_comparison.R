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

n <- 1000
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

sigma_sq <- 2
phi_sq <- 3
tau_sq <- 0.9

d <- dist(coords) %>%
  as.matrix()

w <- rmvnorm(1, rep(0, n), sigma_sq * exp(-phi_sq * d)) %>%
  c()
eps <- rnorm(n, sqrt(tau_sq))
y <- w + eps

# Conditional distribution p(y_star | y) from a MVN
gp_pred <- function(y_new, y, p, coords){
  d <- dist(coords) %>%
    as.matrix()
  d_star <- rbind(y_new, coords) %>%
    dist() %>%
    as.matrix()
  d_star <- d_star[-1, 1]
  sigma_sq <- p[1]
  phi_sq <- p[2]
  tau_sq <- p[3]
  sigma <- sigma_sq * exp(-phi_sq * d) + tau_sq * diag(length(y))
  sigma_star <- sigma_sq * exp(-phi_sq * d_star)
  sigma_star_star <- sigma_sq
  ynew_mean <- sigma_star %*% solve(sigma, y)
  ynew_sigma <- sigma_star_star - sigma_star %*% solve(sigma, sigma_star)
  y_param <- c(ynew_mean, ynew_sigma)
  return(y_param)
}

# GP likelihood
gp_ll = function(p, x, y) {
  sigma_sq <- p[1] %>%
    as.numeric()
  phi_sq <- p[2] %>%
    as.numeric()
  tau_sq <- p[3] %>%
    as.numeric()
  d <- dist(x) %>%
    as.matrix()
  sigma <- sigma_sq * exp(-phi_sq * d) + tau_sq * diag(nrow(x))
  mu <- rep(0, nrow(x))
  ll <- dmvnorm(y, mu, sigma, log = TRUE)
  return(-ll)
}

# NNGP likelihood
nngp_ll = function(p, x, y, m) {
  sigma_sq <- p[1] %>%
    as.numeric()
  phi_sq <- p[2] %>%
    as.numeric()
  tau_sq <- p[3] %>%
    as.numeric()
  nn_data <- get_nn(x, m)
  cond_dens <- NULL
  cond_dens[1] <- dnorm(y[1], 0, sqrt(sigma_sq + tau_sq))
  # Get the conditional distributions for NNGP
  for (i in 1:nrow(nn_data$nn_idx)) {
    nn <- nn_data$ord[nn_data$nn_idx[i,]]
    parms <- gp_pred(
      x[i + 1,], y[nn], c(sigma_sq, phi_sq, tau_sq), matrix(x[nn,], nrow = length(nn)))
    cond_dens[i + 1] <- dnorm(y[i + 1], parms[1], sqrt(parms[2]))
  }
  ll <- prod(cond_dens) %>%
    log()
  return(-ll)
}

parms_grid <- seq(0.2, sigma_sq + 2, length.out = 20) %>%
  tibble(sigma_sq = ., phi_sq, tau_sq)

gp_uni <- apply(parms_grid, MARGIN = 1, FUN = gp_ll, coords, y)
nngp_uni <- apply(parms_grid, MARGIN = 1, FUN = nngp_ll, coords, y, m)

parms_grid <- parms_grid %>%
  bind_cols("gp" = gp_uni, "nngp" = nngp_uni)

ggplot(parms_grid, aes(x = sigma_sq, y = gp, colour = "GP")) +
  geom_line() +
  geom_line(aes(y = nngp, colour = "NNGP")) +
  geom_vline(xintercept = parms_grid$sigma_sq[which.min(parms_grid$gp)], colour = "red") +
  geom_vline(xintercept = parms_grid$sigma_sq[which.min(parms_grid$nngp)], colour = "blue") +
  labs(x = expression(sigma^2), y = "log-likelihood") +
  theme_light()

