# Comparing the joint likelihood of GP and NNGP
# p(y_1) * p(y_2 | y_1) * p(y_3 | y_2, y_1) * p(y_n | y_{n-1}, ..., y_1)

library(dplyr)
source("./nngp-notes/codes/nn_matrix.R")

n <- 5
m <- 3
coords <- cbind(runif(n), runif(n))
d <- dist(coords) %>%
  as.matrix()

# Exponential covariance function
sigma_sq <- 2
phi <- 3
Sigma <- sigma_sq * exp(-phi * d)

plot(seq(-3, 5, by = 0.1), dnorm(seq(-3, 5, by = 0.1), 1.2, 0.5), "l")

# ------------------------------------------
# Full GP
# -----------------------------------------
gp_density <- function(x, sigma, seed = NULL){
  y <- NULL
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  y[1] <- dnorm(x[1], sd = sqrt(sigma[1, 1]))
  for (i in 2:nrow(sigma)) {
    ynew_mean <- sigma[i, 1:(i - 1)] %*% solve(sigma[1:(i - 1), 1:(i - 1)], y)
    ynew_sigma <- sigma[i, i] - sigma[i, 1:(i - 1)] %*% 
      solve(sigma[1:(i - 1), 1:(i - 1)], sigma[i, 1:(i - 1)])
    y[i] <- dnorm(x[i], ynew_mean, sqrt(ynew_sigma))
  }
  return(y)
}

gp_density(rep(0, 5), sigma)

# ------------------------------------------
# NNGP
# -----------------------------------------
# order coords
# knn considering previous subsets
# zero the cov for other vertices

# Creating the matrix by scratch
# m_nearest_neighboor <- function(coords, m){
#   ord <- order(coords[,1])
#   coords_ord <- coords[ord,]
#   d <- dist(coords) %>%
#     as.matrix()
# }
# 
# m_nearest_neighboor <- function(idx, d, m){
#   if (idx < m) {
#     l <- idx
#   } else{
#     l <- m
#   }
#   diag(d) <- Inf
# }

nngp_density <- function(x, sigma, seed = NULL){
  y <- NULL
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  y[1] <- dnorm(x[1], sd = sqrt(sigma[1, 1]))
  for (i in 2:nrow(sigma)) {
    ynew_mean <- sigma[i, 1:(i - 1)] %*% solve(sigma[1:(i - 1), 1:(i - 1)], y)
    ynew_sigma <- sigma[i, i] - sigma[i, 1:(i - 1)] %*% 
      solve(sigma[1:(i - 1), 1:(i - 1)], sigma[i, 1:(i - 1)])
    y[i] <- dnorm(x[i], ynew_mean, sqrt(ynew_sigma))
  }
  return(y)
}
