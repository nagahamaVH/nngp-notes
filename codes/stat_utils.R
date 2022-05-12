# Common utility functions for statistics

library(dplyr)

confint_par <- function(est, heiss, alpha = 0.05) {
  se <- solve(heiss) %>%
    diag() %>%
    sqrt()
  q <- qnorm(alpha / 2) %>% 
    abs()
  ci <- tibble(
    lb = est - q * se,
    ub = est + q * se
  )
  return(ci)
}