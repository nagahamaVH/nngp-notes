library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
source("./nngp-notes/codes/poisson_functs.R")

# -------------------------------------------------
# Gaussian approximation of posterior:
# Checking predicted vs true random effects
# -------------------------------------------------
n <- 300

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
y <- rpois(n, exp(w))

hist(y)

table(y == 0) %>%
  prop.table()

sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
w_init <- rep(0, n)

gauss <- estim_gaussian_proxy(w_init, y, sigma_w)

gauss_df <- tibble(
  id = 1:n,
  w = w,
  w_hat = gauss$mu,
  resid = w - w_hat,
  std_resid = resid / sd(resid)
)

ggplot(gauss_df, aes(x = w, y = w_hat)) +
  geom_abline(slope = 1, col = "red") +
  geom_point()

ggplot(gauss_df, aes(x = id, y = std_resid)) +
  geom_hline(yintercept = 0, col = "red") +
  geom_point()

# -------------------------------------------------
# Gaussian approximation of posterior:
# Checking if the algorithm is finding the mode of posterior
# -------------------------------------------------
n <- 2

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
y <- rpois(n, exp(w))

sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, sigma_w)

df_plot <- NULL
w_grid <- list()
h <- 0.05
for (i in 1:n) {
  se <- fit_proxy$sigma[i, i] %>%
    sqrt()
  w_grid[[i]] <- seq(fit_proxy$mu[i] - 3 * se, fit_proxy$mu[i] + 3 * se, by = h)
}

df_plot <- expand.grid(w_grid) %>%
  as_tibble()
names(df_plot) <- paste0("w", 1:n)

unnorm_exact <- apply(df_plot, 1, posterior, y = y, sigma = sigma_w)
unnorm_proxy <- apply(df_plot, 1, dmvnorm, fit_proxy$mu, fit_proxy$sigma) 

df_plot <- df_plot %>%
  mutate(
    exact = unnorm_exact / sum(unnorm_exact * h^2),
    proxy = unnorm_proxy / sum(unnorm_proxy * h^2)
  )

df_plot %>%
  pivot_longer(cols = c(exact, proxy), names_to = "type", values_to = "density") %>%
  ggplot(., aes(x = w1, y = w2)) +
  geom_contour(aes(z = density)) +
  facet_wrap(vars(type)) +
  geom_point(aes(x = w[1], y = w[2]), pch = 3)

# -------------------------------------------------
# Gaussian approximation of posterior:
# Checking the marginal distribution
# -------------------------------------------------
n <- 300

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
y <- rpois(n, exp(w))

sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, sigma_w)

h <- 0.05
n_sample <- 6
set.seed(523)
sampled <- sample(1:n, n_sample, replace = F)
plots <- list()
for (i in 1:n_sample) {
  id <- sampled[i]
  se <- fit_proxy$sigma[id, id] %>%
    sqrt()
  w_grid <- seq(fit_proxy$mu[id] - 3 * se, fit_proxy$mu[id] + 3 * se, by = h)
  
  df_plot <- rep(w, each = length(w_grid)) %>%
    matrix(ncol = n) %>%
    as_tibble()
  df_plot[,id] <- w_grid
  names(df_plot) <- paste0("w", 1:n)
  
  unnorm_exact <- apply(df_plot, 1, posterior, y = y, sigma = sigma_w)
  unnorm_proxy <- apply(df_plot, 1, dmvnorm, fit_proxy$mu, fit_proxy$sigma) 
  
  df_plot <- df_plot %>%
    mutate(
      exact = unnorm_exact / sum(unnorm_exact * h^2),
      proxy = unnorm_proxy / sum(unnorm_proxy * h^2)
    )
  
  # p <- ggplot(df, aes(x = b, y = exact)) +
  #   geom_line() +
  #   geom_line(aes(y = proxy), linetype = 2) +
  #   geom_point(aes_string(x = b0_grid[i], y = 0)) +
  #   labs(y = "Density") +
  #   xlim(c(-2, 3)) +
  #   theme_light()
  
  p <- ggplot(df_plot, aes_string(x = paste0("w", id))) +
    geom_line(aes(y = exact)) +
    geom_line(aes(y = proxy), linetype = 2, col = "red") +
    geom_point(aes(x = w[id], y = 0))
  plots[[i]] <- p
}

plot_grid(plotlist = plots, nrow = 2)

# -------------------------------------------------
# Checking Laplace proxy
# -------------------------------------------------
n <- 300

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 1
phi_true <- 2.3

# Generate data
d <- dist(coords) %>%
  as.matrix()
w <- rmvnorm(1, sigma = sigma_true^2 * exp(-phi_true^2 * d)) %>%
  c()
y <- rpois(n, exp(w))

hist(y)
table(y == 0) %>%
  prop.table()

parms_grid <- tibble(
  log_sigma = seq(log(0.6), log(1.5), length.out = 10),
  log_phi = seq(log(1.7), log(4), length.out = 10)
) %>%
  expand.grid()

w_init <- rep(0, n)
sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
fit_proxy <- estim_gaussian_proxy(w_init, y, sigma_w)

ll <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_hat = fit_proxy$mu, 
  sigma_hat = fit_proxy$sigma, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(ll)

ggplot(parms_grid, aes(x = log_sigma, y = log_phi)) +
  geom_contour_filled(aes(z = ll)) +
  geom_point(data = parms_grid[which.max(parms_grid$ll),], size = 5) +
  geom_point(
    data = tibble(sigma_true, phi_true), 
    aes(x = log(sigma_true), y = log(phi_true)),
    shape = 3, size = 5
  )

parms_grid %>%
  filter(ll == max(ll)) %>%
  mutate(
    sigma = exp(log_sigma),
    phi = exp(log_phi),
    sigma_true,
    phi_true
  )
