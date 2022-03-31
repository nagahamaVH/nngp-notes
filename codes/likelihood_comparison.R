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
library(cowplot)
library(mvtnorm)
source("./nngp-notes/codes/gp_nngp_ll.R")

n <- 50
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma <- 3
phi <- 4
tau <- 0.6

d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma^2 * exp(-phi^2 * d)) %>%
  c()
eps <- rnorm(n, 0, tau)
y <- w + eps

# -----------------------
# Univariate plot
# -----------------------
parms_grid <- seq(sigma - .5, sigma + .5, length.out = 100) %>%
  tibble(sigma = ., phi, tau)

gp <- apply(parms_grid, MARGIN = 1, FUN = gp_ll, coords, y)
nngp <- apply(parms_grid, MARGIN = 1, FUN = nngp_ll, coords, y, m)

parms_grid <- parms_grid %>%
  bind_cols("gp" = gp, "nngp" = nngp)

ggplot(parms_grid, aes(x = sigma, y = gp, colour = "GP")) +
  geom_vline(xintercept = sigma, linetype = 2) +
  geom_line() +
  geom_line(aes(y = nngp, colour = "NNGP")) +
  geom_point(data = parms_grid[which.min(parms_grid$gp),]) +
  geom_point(
    data = parms_grid[which.min(parms_grid$nngp),],
    aes(y = nngp, colour = "NNGP")) +
  labs(
    x = expression(sigma), y = "log-likelihood", colour = "",
    title = paste0("n = ", n, ", m = ", m)) +
  theme_light()

# -----------------------
# Bivariate plot
# -----------------------
parms_grid2 <- tibble(
  sigma = seq(sigma - 1, sigma + 2, length.out = 10), 
  phi = seq(phi - 1, phi + 2, length.out = 10), 
  tau) %>%
  expand.grid()

gp <- apply(parms_grid2, MARGIN = 1, FUN = gp_ll, coords, y)
nngp <- apply(parms_grid2, MARGIN = 1, FUN = nngp_ll, coords, y, m)

parms_grid2 <- parms_grid2 %>%
  bind_cols("gp" = gp, "nngp" = nngp)

p1 <- ggplot(parms_grid2, aes(x = phi, y = sigma, colour = gp)) +
  geom_point(size = 8) +
  geom_point(
    aes(x = phi, y = sigma), data = tibble(phi, sigma), 
    color = "black", size = 8, shape = 18) +
  geom_point(
    aes(x = phi, y = sigma), data = parms_grid2[which.min(parms_grid2$gp),], 
    color = "black", size = 10, shape = 1) +
  labs(
    x = expression(Phi), y = expression(sigma), colour = "log-likelihood",
    title = "GP") +
  scale_colour_viridis_b(direction = -1) +
  theme_light()

p2 <- ggplot(parms_grid2, aes(x = phi, y = sigma, colour = nngp)) +
  geom_point(size = 8) +
  geom_point(
    aes(x = phi, y = sigma), data = tibble(phi, sigma), 
    color = "black", size = 8, shape = 18) +
  geom_point(
    aes(x = phi, y = sigma), data = parms_grid2[which.min(parms_grid2$nngp),], 
    color = "black", size = 10, shape = 1) +
  labs(
    x = expression(Phi), y = expression(sigma), colour = "log-likelihood",
    title = "NNGP") +
  scale_colour_viridis_b(direction = -1) +
  theme_light()

legend <- get_legend(p1)

plot_grid(
  p1 + theme(legend.position = "none"), 
  p2 + theme(legend.position = "none"),
  legend, ncol = 3, rel_widths = c(1, 1, .4))

# -----------------------
# Predicted surface
# -----------------------
df <- tibble(
  c1 = coords[,1],
  c2 = coords[,2],
  y = y
)

gp_pred <- df %>%
  select(c1, c2) %>%
  apply(MARGIN = 1, cond_mvn, y = y, p = c(sigma, phi, tau), coords = coords) %>%
  t()

df <- df %>%
  bind_cols(gp = gp_pred[,1])

p3 <- ggplot(df, aes(x = c1, y = c2, z = y)) +
  geom_density_2d_filled() +
  theme_light()

p4 <- ggplot(df, aes(x = c1, y = c2, z = gp)) +
  geom_density_2d_filled() +
  theme_light()

plot_grid(p3, p4)
