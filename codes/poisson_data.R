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
library(cowplot)
library(mvtnorm)
library(microbenchmark)
source("./nngp-notes/codes/gp_nngp_ll.R")

generate_data <- function(n, seed = NULL) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n))
  
  # To do: make generic to data without ordering
  ord <- order(coords[,1])
  coords <- coords[ord,]
  
  d <- dist(coords) %>%
    as.matrix()
  w <- rmvnorm(1, sigma = sigma^2 * exp(-phi^2 * d)) %>%
    c()
  eps <- rnorm(n, 0, tau)
  y <- w + eps
  l <- list(y = y, w = w, coords = coords)
  return(l)
}

confint_par <- function(l) {
  se <- solve(l$hessian) %>%
    diag() %>%
    sqrt()
  ci <- tibble(
    lb = l$par - 1.96 * se,
    ub = l$par + 1.96 * se
  )
  return(ci)
}

n <- 5
m <- 3

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma <- 1.7
phi <- 3.2

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma^2 * exp(-phi^2 * d)) %>%
  c()
# C <- sigma^2 * exp(-phi^2 * d)
# w <- c(t(chol(C)) %*% rnorm(n))
y <- rpois(n, exp(w))

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
    x = expression(sigma), y = "neg log-likelihood", colour = "",
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
    x = expression(Phi), y = expression(sigma), colour = "neg log-likelihood",
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
    x = expression(Phi), y = expression(sigma), colour = "neg log-likelihood",
    title = "NNGP") +
  scale_colour_viridis_b(direction = -1) +
  theme_light()

legend <- get_legend(p1)

plot_grid(
  p1 + theme(legend.position = "none"), 
  p2 + theme(legend.position = "none"),
  legend, ncol = 3, rel_widths = c(1, 1, .4))

# -----------------------
# Benchmarking of processing time (this can take a while)
# -----------------------
size_grid <- expand.grid(
  n = c(10, 50, 100, 1000),
  m = c(3, 10, 50)
) %>%
  filter(n > m) %>%
  arrange(n)
init_parms <- rep(1, 3)

bench <- list()
for (i in 1:nrow(size_grid)) {
  paste0(i, "/", nrow(size_grid), " | n: ", size_grid$n[i], " m: ", size_grid$m[i]) %>%
    print()
  df <- generate_data(size_grid$n[i])
  b <- microbenchmark(
    times = 10, 
    unit = "milliseconds",
    "gp" = {
      gp <- optim(
        init_parms, gp_ll, coords = df$coords, y = df$y, method = "L-BFGS-B",
        lower = rep(10e-5, 3))
    },
    "nngp" = {
      nngp <- optim(
        init_parms, nngp_ll, coords = df$coords, y = df$y, m = size_grid$m[i], 
        method = "L-BFGS-B", lower = rep(10e-5, 3))
    }
  )
  bench[[i]] <- b 
}

# Saving data
# saveRDS(bench, "./nngp-notes/data/bench.rds")
# saveRDS(size_grid, "./nngp-notes/data/size_grid.rds")

# Loading data
bench <- readRDS("./nngp-notes/data/bench.rds")
size_grid <- readRDS("./nngp-notes/data/size_grid.rds")

p <- list()
for (i in 1:nrow(size_grid)) {
  title <- paste0("n=", size_grid$n[i], " m=", size_grid$m[i])
  p[[i]] <- autoplot(bench[[i]]) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}
plot_grid(plotlist = p)

# -----------------------
# Comparing estimate of spatial effect
# -----------------------
size_grid2 <- size_grid %>%
  filter(n == max(size_grid$n))

df <- generate_data(size_grid2$n[1], 523)
list_name <- paste0("nngp", size_grid2$m)

estim <- list()
gp <- optim(
  init_parms, gp_ll, coords = df$coords, y = df$y, method = "L-BFGS-B",
  lower = rep(10e-5, 3), hessian = T)
estim[["gp"]] <- gp

for (i in 1:nrow(size_grid2)) {
  paste0(i, "/", nrow(size_grid2), " | n: ", size_grid2$n[i], " m: ", size_grid2$m[i]) %>%
    print()
  nngp <- optim(
    init_parms, nngp_ll, coords = df$coords, y = df$y, m = size_grid$m[i], 
    method = "L-BFGS-B", lower = rep(10e-5, 3), hessian = T)
  estim[[list_name[i]]] <- nngp
}

preds <- list()
for (i in 1:length(estim)) {
  preds[[i]] <- apply(
    df$coords, MARGIN = 1, pred_cond_mvn, x = df$coords, y = df$y, 
    p = estim[[i]]$par) %>%
    t()
}

df_plot <- tibble(
  c1 = df$coords[,1],
  c2 = df$coords[,2],
  w = df$w,
  gp = preds[[1]][,1],
  nngp1 = preds[[2]][,1],
  nngp2 = preds[[3]][,1],
  nngp3 = preds[[4]][,1]
)

p_w <- ggplot(df_plot, aes(x = c1, y = c2, z = w)) +
  geom_density_2d_filled() +
  labs(title = "True W", x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_gp <- ggplot(df_plot, aes(x = c1, y = c2, z = gp)) +
  geom_density_2d_filled() +
  labs(title = "GP", x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_nngp1 <- ggplot(df_plot, aes(x = c1, y = c2, z = nngp1)) +
  geom_density_2d_filled() +
  labs(title = "NNGP 3", x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_nngp2 <- ggplot(df_plot, aes(x = c1, y = c2, z = nngp2)) +
  geom_density_2d_filled() +
  labs(title = "NNGP 10", x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_nngp3 <- ggplot(df_plot, aes(x = c1, y = c2, z = nngp3)) +
  geom_density_2d_filled() +
  labs(title = "NNGP 50", x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

plot_grid(
  p_w + theme(legend.position = "none"), 
  p_gp + theme(legend.position = "none"), 
  p_nngp1 + theme(legend.position = "none"), 
  p_nngp2 + theme(legend.position = "none"), 
  p_nngp3 + theme(legend.position = "none"))

# ggsave(
#   "./nngp-notes/images/spatial_estimate.png", plot = last_plot(), 
#   width = 10, height = 8)

# -----------------------
# Comparing parameter estimates
# -----------------------
ci <- lapply(estim, confint_par) %>%
  do.call(bind_rows, .) %>%
  mutate(
    name = rep(c("sigma", "phi", "tau"), length(preds)),
    type = rep(c("gp", "nngp3", "nngp10", "nngp50"), 3)
  )

ci <- map(estim, ~ .x[1]) %>%
  do.call(bind_rows, .) %>%
  bind_cols(ci, .)

p_sigma <- filter(ci, name == "sigma") %>%
  ggplot(., aes(y = par, x = type, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = sigma, linetype = 2) +
  geom_pointrange() +
  labs(title = expression(sigma), x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_phi <- filter(ci, name == "phi") %>%
  ggplot(., aes(y = par, x = type, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = phi, linetype = 2) +
  geom_pointrange() +
  labs(title = expression(Phi), x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

p_tau <- filter(ci, name == "tau") %>%
  ggplot(., aes(y = par, x = type, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = tau, linetype = 2) +
  geom_pointrange() +
  labs(title = expression(tau), x = "", y = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

plot_grid(p_sigma, p_phi, p_tau, nrow = 1)

# ggsave(
#   "./nngp-notes/images/par_estimate.png", plot = last_plot(),
#   width = 12, height = 6)
