# https://mbjoseph.github.io/posts/2018-12-27-gaussian-predictive-process-models-in-stan/
# https://mc-stan.org/users/documentation/case-studies/nngp.html
# https://mc-stan.org/docs/2_23/stan-users-guide/posterior-predictive-simulation-in-stan.html

library(rstan)
library(bayesplot)
source("./codes/NNMatrix.R")

# ------------------- Setup ---------------------------------------------------
options(mc.cores = parallel::detectCores())

model_name <- "poisson_nngp_zhang"

data_board <- pins::board_folder("./data", versioned = T)
model_board <- pins::board_folder("models", versioned = T)
# -----------------------------------------------------------------------------

n <- 1000
max_s <- 1.5

set.seed(163)
coords <- cbind(runif(n, max = max_s), runif(n, max = max_s))

ord <- order(coords[,1])
coords <- coords[ord,]

sigma <- .9
l <- .4
beta <- c(3, .5)

D <- dist(coords) |>
  as.matrix()
C <- sigma^2 * exp(-D / (2 * l^2))
x <- mvtnorm::rmvnorm(1, sigma = C, method = "chol", checkSymmetry = F) |>
  c()
Z <- cbind(1, rnorm(n))
eta <- exp(tcrossprod(Z, t(beta)) + x)
y <- rpois(n, eta)

m <- 5

nn <- NNMatrix(coords, m)

stan_data <- list(
  N = n,
  Y = y,
  X = Z,
  P = dim(Z)[2],
  NN_ind = nn$NN_ind,
  NN_dist = nn$NN_dist,
  NN_distM = nn$NN_distM,
  M = m)

# stan_data <- list(
#   n = n,
#   y = y,
#   Z = Z,
#   p = dim(Z)[2],
#   nn = nn$NN_ind,
#   d_pairs = nn$NN_dist,
#   d_nn_pairs = nn$NN_distM,
#   m = m)

hist(y)

# ------------------------ Stan parameters ------------------------------------
n_chain <- 4
n_it <- 1500
model_file <- "./codes/poisson_nngp_zhang.stan"
# -----------------------------------------------------------------------------

t_init <- proc.time()
stan_fit <- stan(
  file = model_file,
  data = stan_data,
  chains = n_chain,
  iter = n_it,
  seed = 171
)
t_total <- proc.time() - t_init

# Save data
data_meta <- list(
  cov_kernel = "se",
  n = n,
  max_s = max_s,
  cov_par = list(
    sigma = sigma,
    l = l
  ),
  beta = beta
)

pins::pin_write(
  data_board, stan_data, model_name, type = "rds", metadata = data_meta)

# Save model
data_version <- data_board |>
  pins::pin_versions(model_name) |>
  tail(1) |>
  dplyr::pull(version)

model_meta <- list(
  model = model_name,
  file = model_file,
  data = data_version,
  n_it = n_it,
  n_chain = n_chain,
  time = t_total,
  model_code = readLines(model_file) |>
    paste(collapse = "\n")
)

pins::pin_write(
  model_board, stan_fit, model_name, type = "rds", metadata = model_meta)
