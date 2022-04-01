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

# write as linear combination
# Cholesky(as(Sigma, "dsCMatrix"), LDL = TRUE)
# C^{-1} = (I - A)^T D^{-1} (I - A)

A <- matrix(rep(0, n^2), nrow = n)
D <- diag(n)
D[1, 1] <- Sigma[1, 1]
C <- Sigma
for(i in 1:(n-1)){
  A[i+1,1:i] = solve(C[1:i,1:i], C[1:i,i+1])
  D[i+1,i+1] = C[i+1,i+1] - sum(C[i+1,1:i] * A[i+1,1:i])
}
