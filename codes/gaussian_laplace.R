library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
library(msm)
source("./nngp-notes/codes/stat_utils.R")

n <- 100

set.seed(126)
coords <- cbind(runif(n), runif(n))

# To do: make generic to data without ordering
ord <- order(coords[,1])
coords <- coords[ord,]

# Parameters
sigma_true <- 3
phi_true <- 2
tau_true <- .1

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
  Sigma <- sigma^2 * exp(-phi^2 * d) + tau^2 * diag(nrow(coords))
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
    sigma <- solve(solve(sigma_w) - hess_m)
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
laplace <- function(y, coords, p_init, w_hat, sigma_hat, mu_w = rep(0, length(y))) {
  sigma <- exp(p_init[1])
  phi <- exp(p_init[2])
  tau <- exp(p_init[3])
  d <- dist(coords) %>%
    as.matrix()
  sigma_w <- sigma^2 * exp(-phi^2 * d)

  denom <- -(length(y) / 2) * log(2 * pi) - 0.5 * log(det(sigma_hat) + sqrt(.Machine$double.xmin))
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
    as_tibble()
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
# Checking Laplace approximation
# Ploting the likelihood function with respect of parameters
# given the true parameters and mode estimation of latent field 
# -------------------------------------------------
# -----------------------------------------------------------
# Sigma
# -----------------------------------------------------------
parms_grid <- tibble(
  log_sigma = seq(1, 8, by = 0.01),
  phi = phi_true,
  tau = tau_true
) %>%
  mutate_all(log)

# Assuming true parameters
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, tau_true, sigma_w)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)

la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_hat = fit_proxy$mu, 
  sigma_hat = fit_proxy$sigma, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  pivot_longer(cols = c(ll, la), names_to = "type", values_to = "ll") %>%
  ggplot(., aes(x = log_sigma, y = ll)) +
  geom_line() +
  facet_wrap(vars(type), scales = "free")

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
  log_phi = seq(.5, 4, by = 0.01),
  tau = tau_true
) %>%
  mutate_all(log)

# Assuming true parameters
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, tau_true, sigma_w)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)

la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_hat = fit_proxy$mu, 
  sigma_hat = fit_proxy$sigma, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  pivot_longer(cols = c(ll, la), names_to = "type", values_to = "ll") %>%
  ggplot(., aes(x = log_phi, y = ll)) +
  geom_line() +
  facet_wrap(vars(type), scales = "free")

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
  log_tau = seq(.005, 1, by = 0.01)
) %>%
  mutate_all(log)

# Assuming true parameters
w_init <- rep(0, n)
fit_proxy <- estim_gaussian_proxy(w_init, y, tau_true, sigma_w)

ll <- apply(parms_grid, 1, gp_ll, coords = coords, y = y)

la <- apply(
  parms_grid, 1, laplace, y = y, coords = coords, w_hat = fit_proxy$mu,
  sigma_hat = fit_proxy$sigma, mu_w = rep(0, n))

parms_grid <- parms_grid %>%
  mutate(la, ll)

parms_grid %>%
  pivot_longer(cols = c(ll, la), names_to = "type", values_to = "ll") %>%
  ggplot(., aes(x = log_tau, y = ll)) +
  geom_line() +
  facet_wrap(vars(type), scales = "free")

parms_grid %>%
  filter(ll == max(ll)) %>%
  mutate(tau = exp(log_tau))

parms_grid %>%
  filter(la == max(la)) %>%
  mutate(tau = exp(log_tau))

ggplot(parms_grid, aes(x = log_tau)) +
  geom_line(aes(y = ll)) +
  geom_line(aes(y = la), linetype = 2, color = "red")

# Check if works
gp <- optim(
  log(c(1, 1, 1)), gp_ll, coords = coords, y = y, method = "BFGS",
  control = list(fnscale = -1), hessian = T)

sigma_hat <- deltamethod(
  list(~exp(x1), ~exp(x2), ~exp(x3)), gp$par, -solve(gp$hessian), ses = F)

confint2(exp(gp$par), sigma_hat) %>%
  mutate(
    par = c("sigma", "phi", "tau"),
    true = c(sigma_true, phi_true, tau_true)
  )

fit <- function() {
  n_it <- 100
  p_init <- c(log(1), log(1), log(1))
  w_init <- rep(0, n)
  mu_w <- 0
  for (i in 1:n_it) {
    print(paste0("it: ", i))
    sigma <- exp(p_init[1])
    phi <- exp(p_init[2])
    tau <- exp(p_init[3])
    sigma_w <- sigma^2 * exp(-phi^2 * d)
    
    # Gaussian proxy
    gauss <- estim_gaussian_proxy(w_init, y, tau, sigma_w)
    w_hat <- gauss$mu
    sigma_hat <- gauss$sigma

    # Laplace
    fitted_parms <- optim(
      p_init, laplace, y = y, coords = coords, w_hat = w_hat,
      sigma_hat = sigma_hat, mu_w = rep(mu_w, n), method = "BFGS", hessian = T,
      control = list(fnscale = -1))
    if (fitted_parms$convergence != 0) stop("Convergence error Laplace")
    
    e <- abs(p_init - fitted_parms$par)
    if (max(e) < 1e-4) break
    if (i == n_it) stop("Max iteraction number reached")
    p_init <- fitted_parms$par
    w_init <- w_hat
    print(paste0("parms: ", exp(fitted_parms$par)))
  }
  return(fitted_parms)
}

fitted_model <- fit()

exp(fitted_model$par)

sigma_hat <- deltamethod(
  list(~exp(x1), ~exp(x2), ~exp(x3)), fitted_model$par, -solve(fitted_model$hessian), ses = F)

ci_parms <- confint2(exp(fitted_model$par), sigma_hat) %>%
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
