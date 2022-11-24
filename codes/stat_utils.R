# Common utility functions for statistics

confint_par <- function(par, hess, alpha = 0.05) {
  q <- abs(qnorm(alpha / 2))
  se <- sqrt(diag(solve(hess)))

  ci <- data.frame(
    lb = par - q * se,
    ub = par + q * se
  )
  return(ci)
}
