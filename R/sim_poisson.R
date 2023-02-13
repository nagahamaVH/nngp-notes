board <- pins::board_folder("./data", versioned = T)

# sample_size <- c(100, 500, 1000, 10000)
sample_size <- c(10, 20)

max_s <- 1
sigma <- .9
l <- .4
beta <- c(3, .5)

set.seed(163)

for (i in 1:length(sample_size)) {
  n <- sample_size[i]
  coords <- cbind(runif(n, max = max_s), runif(n, max = max_s))
  
  ord <- order(coords[,1])
  coords <- coords[ord,]
  D <- dist(coords) |>
    as.matrix()
  C <- sigma^2 * exp(-D / (2 * l^2))
  x <- mvtnorm::rmvnorm(1, sigma = C, method = "chol", checkSymmetry = F) |>
    c()
  Z <- cbind(1, rnorm(n))
  eta <- exp(tcrossprod(Z, t(beta)) + x)
  y <- rpois(n, eta)
  
  data <- list(
    N = n,
    Y = y,
    X = Z,
    P = dim(Z)[2],
    coords = coords
  )
  
  meta <- list(
    model = "y ~ Poisson(exp(b0 + b1 * z + x))",
    cov_function = "SE",
    n = n,
    max_s = max_s,
    par = list(
      cov_par = list(
        sigma = sigma,
        l = l
      ),
      beta = beta 
    ),
    code = readLines("R/sim_poisson.R") |>
      paste(collapse = "\n")
  )

  pins::pin_write(
    board, data, paste("poisson-gp", n, sep = "_"), type = "rds", metadata = meta)
}
