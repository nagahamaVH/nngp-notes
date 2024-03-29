library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
library(msm)
source("./nngp-notes/codes/stat_utils.R")

n <- 100

set.seed(126)
# coords <- cbind(runif(n), runif(n))
coords <- cbind(rep(1:10, 10), rep(1:10, each = 10))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 3
phi_true <- .5
tau_true <- 1

# Generate data
d <- dist(coords) %>%
  as.matrix()
sigma_w <- sigma_true^2 * exp(-phi_true^2 * d)
w <- rmvnorm(1, sigma = sigma_w) %>%
  c()
y <- rnorm(n, w, tau_true)

# Decay
plot(d[,1], sigma_w[,1])

# -----------------------------------------------------------------------
gp_ll = function(p, coords, y) {
  sigma <- p[1] %>%
    as.numeric() %>%
    exp()
  phi <- p[2] %>%
    as.numeric() %>%
    exp()
  tau <- p[3] %>%
    as.numeric() %>%
    exp()
  d <- dist(coords) %>%
    as.matrix()
  Sigma <- sigma^2 * exp(-phi^2 * d) + diag(tau^2, nrow(coords))
  mu <- rep(0, nrow(coords))
  ll <- dmvnorm(y, mu, Sigma, log = TRUE)
  return(ll)
}

posterior <- function(y, w, sd_y, sigma_w, mu = rep(0, length(w)), log = F) {
  ll <- sum(dnorm(y, w, sd_y, log = T)) + dmvnorm(w, mu, sigma_w, log = T)
  if (!log) {
    ll <- exp(ll)
  }
  return(ll)
}

gradient <- function(w, y, tau) {
  1 / tau^2 * (y - w)
}

hessian <- function(tau) {
  -1 / tau^2
}

estim_gaussian_proxy <- function(
    w0, y, tau, sigma_w, control = list(it = 1000, tol = 1e-4)) {
  for (i in 1:control$it) {
    if (any(hessian(tau) > 0)) {
      print(hessian(tau))
      stop("hessian error")
    }

    hess_m <- diag(hessian(tau), length(w0))
    sigma <- chol2inv(chol((chol2inv(chol(sigma_w)) - hess_m)))
    mu <- as.vector(sigma %*% (gradient(w0, y, tau) - hess_m %*% w0))

    if (max(abs(w0 - mu)) < control$tol) break
    if (i == control$it) stop("Max iteraction number reached")
    w0 <- mu
  }

  parms <- list(
    mu = mu, sigma = sigma)
  return(parms)
}

# Laplace approximation
# f(theta | y) \approx f(y | w_hat, theta) * f(w_hat | theta) / MVN(w_hat, sigma)
laplace <- function(y, coords, p_init, w_init, mu_w = rep(0, length(y))) {
  sigma <- exp(p_init[1])
  phi <- exp(p_init[2])
  tau <- exp(p_init[3])
  d <- dist(coords) %>%
    as.matrix()
  sigma_w <- sigma^2 * exp(-phi^2 * d)
  
  # Gaussian proxy
  gauss <- estim_gaussian_proxy(w_init, y, tau, sigma_w)
  w_hat <- gauss$mu
  sigma_hat <- gauss$sigma
  
  denom <- -(length(y) / 2) * log(2 * pi) - 0.5 * determinant(sigma_hat)$modulus
  ll <- posterior(y, w_hat, tau, sigma_w, mu_w, log = T) - denom
  
  if (ll == -Inf) {
    ll <- -sqrt(.Machine$double.xmax)
  }
  
  return(ll)
}

# -------------------------------------------------
# Gaussian approximation of posterior:
# Checking the marginal distribution
# -------------------------------------------------
# Assuming true parameters
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, tau_true, sigma_w)

h <- 0.01
n_sample <- 6
set.seed(52)
sampled <- sample(1:n, n_sample, replace = F)
plots <- list()
for (i in 1:n_sample) {
  id <- sampled[i]

  se <- fit_proxy$sigma[id, id] %>%
    sqrt()
  w_grid <- seq(fit_proxy$mu[id] - 3 * se, fit_proxy$mu[id] + 3 * se, by = h)
  
  df_plot <- rep(w, each = length(w_grid)) %>%
    matrix(ncol = n) %>%
    as_tibble(.name_repair = make.names)
  df_plot[,id] <- w_grid
  names(df_plot) <- paste0("w", 1:n)
  
  unnorm_exact <- apply(df_plot, 1, posterior, sd_y = tau_true, y = y, sigma = sigma_w)
  unnorm_proxy <- apply(df_plot, 1, dmvnorm, fit_proxy$mu, fit_proxy$sigma) 
  
  df_plot <- df_plot %>%
    mutate(
      exact = unnorm_exact / sum(unnorm_exact * h^2),
      proxy = unnorm_proxy / sum(unnorm_proxy * h^2)
    )

  p <- ggplot(df_plot, aes_string(x = paste0("w", id))) +
    geom_line(aes(y = exact)) +
    geom_line(aes(y = proxy), linetype = 2, col = "red") +
    geom_point(data = tibble(w = w[id]), aes(x = w, y = 0))
  plots[[i]] <- p
}

plot_grid(plotlist = plots, nrow = 2)

plot(w, fit_proxy$mu); abline(a = 0, b = 1)

# -------------------------------------------------
# Gaussian approximation of posterior:
# Comparing parameters of posterior of latent field and its Gaussian proxy
# -------------------------------------------------
sigma <- .1
phi <- phi_true
tau <- tau_true
w_init <- rep(0, n)
sigma_w <- sigma^2 * exp(-phi^2 * d)
fit_proxy <- estim_gaussian_proxy(w_init, y, tau, sigma_w)

# Exact pi(w | y, theta)
sigma_y <- diag(tau^2, n)
prec_post <- solve(sigma_w) + solve(sigma_y)
sigma_post <- solve(prec_post)
mu_post <- (sigma_post %*% solve(sigma_y) %*% y) %>%
  c()

plot(fit_proxy$mu, mu_post); abline(0, 1)

id <- lower.tri(fit_proxy$sigma, diag = T)
lm(fit_proxy$sigma[id] ~ sigma_post[id])

# -------------------------------------------------
# Laplace approximation
# -------------------------------------------------
# -----------------------------------------------------------
# Sigma
# -----------------------------------------------------------
parms_grid <- tibble(
  log_sigma = seq(.01, 50, by = 1),
  phi = phi_true,
  tau = tau_true
) %>%
  mutate_all(log)

# Assuming true parameters
w_init <- rep(0, n)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)
la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_init = w_init, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  filter(ll == max(ll)) %>%
  mutate(sigma = exp(log_sigma))

parms_grid %>%
  filter(la == max(la)) %>%
  mutate(sigma = exp(log_sigma))

ggplot(parms_grid, aes(x = log_sigma)) +
  geom_line(aes(y = ll)) +
  geom_line(aes(y = la), linetype = 2, color = "red")

# -----------------------------------------------------------
# Phi
# -----------------------------------------------------------
parms_grid <- tibble(
  sigma = sigma_true,
  log_phi = seq(.2, 3, by = 0.01),
  tau = tau_true
) %>%
  mutate_all(log)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)
la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_init = w_init, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  filter(ll == max(ll)) %>%
  mutate(phi = exp(log_phi))

parms_grid %>%
  filter(la == max(la)) %>%
  mutate(phi = exp(log_phi))

ggplot(parms_grid, aes(x = log_phi)) +
  geom_line(aes(y = ll)) +
  geom_line(aes(y = la), linetype = 2, color = "red")

# -----------------------------------------------------------
# Tau
# -----------------------------------------------------------
parms_grid <- tibble(
  sigma = sigma_true,
  phi = phi_true,
  log_tau = seq(.5, 3, by = 0.01)
) %>%
  mutate_all(log)

w_init <- rep(0, n)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)
la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_init = w_init, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  filter(ll == max(ll)) %>%
  mutate(tau = exp(log_tau))

parms_grid %>%
  filter(la == max(la)) %>%
  mutate(tau = exp(log_tau))

ggplot(parms_grid, aes(x = log_tau)) +
  geom_line(aes(y = ll)) +
  geom_line(aes(y = la), linetype = 2, color = "red")

# Bivariate plot
parms_grid <- tibble(
  log_sigma = seq(.1 * sigma_true, 1.5 * sigma_true, length.out = 10),
  log_phi = seq(.1 * phi_true, 1.5 * phi_true, length.out = 10),
  tau = tau_true
) %>%
  mutate_all(log) %>%
  expand.grid()

la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_init = w_init,
  mu_w = rep(0, n))

ll <- apply(
  parms_grid, 1, gp_ll, y = y, coords = coords)

parms_grid <- parms_grid %>%
  mutate(ll, la)

p1 <- ggplot(parms_grid, aes(x = log_sigma, y = log_phi)) +
  geom_contour_filled(aes(z = ll)) +
  geom_point(data = parms_grid[which.max(parms_grid$ll),], size = 5) +
  geom_point(
    data = tibble(sigma_true, phi_true), 
    aes(x = log(sigma_true), y = log(phi_true)),
    shape = 3, size = 5
  ) +
  theme(legend.position = "none")

p2 <- ggplot(parms_grid, aes(x = log_sigma, y = log_phi)) +
  geom_contour_filled(aes(z = la)) +
  geom_point(data = parms_grid[which.max(parms_grid$ll),], size = 5) +
  geom_point(
    data = tibble(sigma_true, phi_true), 
    aes(x = log(sigma_true), y = log(phi_true)),
    shape = 3, size = 5
  ) +
  theme(title = "laplace approximation")

plot_grid(p1, p2, nrow = 1)

# -----------------------------------------------------------
# Checking Laplace approximation estimation
# -----------------------------------------------------------
# Likelihood
gp <- optim(
  log(c(1, 1, 1)), gp_ll, coords = coords, y = y, method = "Nelder-Mead",
  control = list(fnscale = -1), hessian = T)

sigma_hat <- deltamethod(
  list(~exp(x1), ~exp(x2), ~exp(x3)), gp$par, -solve(gp$hessian), ses = F)

confint2(exp(gp$par), sigma_hat) %>%
  mutate(
    par = c("sigma", "phi", "tau"),
    true = c(sigma_true, phi_true, tau_true)
  )

# Laplace approximation
p_init <- rep(1, 3) %>% 
  log()
w_init <- rep(0, n)
mu_w <- 0
fitted_parms <- optim(
  p_init, laplace, y = y, coords = coords, w_init = w_init, mu_w = rep(mu_w, n), 
  method = "Nelder-Mead", hessian = T, control = list(fnscale = -1))

sigma_hat <- deltamethod(
  list(~exp(x1), ~exp(x2), ~exp(x3)), fitted_parms$par, 
  solve(-fitted_parms$hessian), ses = F)

ci_parms <- confint2(exp(fitted_parms$par), sigma_hat) %>%
  mutate(
    par = c("sigma", "phi", "tau"),
    true = c(sigma_true, phi_true, tau_true)
  )

ggplot(data = ci_parms, aes(x = "")) +
  facet_wrap(~par, scales = "free") +
  geom_linerange(aes(ymin = lb, ymax = ub)) +
  geom_point(aes(y = mu)) +
  geom_hline(aes(yintercept = true), lty = 2, col = "red") +
  labs(y = "Estimate", x = "")
