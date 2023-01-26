functions{
  // Reconstruct distance matrix for point i and its neighbours
  // matrix get_D(vector d, vector d_m) {
  //   int m = dims(d)[1];
  //   int n = m + 1;
  //   matrix[n, n] D = rep_matrix(0, n, n);
  //   
  //   D[2:n, 1] = d;
  //   D[1, 2:n] = to_row_vector(D[1, 2:n]);
  //   for (i in 2:m) {
  //     for (j in (i + 1):n) {
  //       D[i, j] = d_m[j - i];
  //       D[j, i] = D[i, j];
  //     }
  //   }
  //   return D;
  // }
  
  real nngp_lpdf(
    vector x, matrix d_pairs, matrix d_nn_pairs, int[,] nn, int n, int m, real sigma, real l) {
    int m_i;
    int d_size = m + 1;
    matrix[n, n] C;
    // matrix[n, n] L;
    matrix[n, n] I = diag_matrix(rep_vector(1, n));
    matrix[d_size, d_size] C_i;
    matrix[n, n] A = rep_matrix(0, n, n);
    vector[n] dd;
    dd[1] = sigma^2;
    
    // Update A and d
    for (i in 1:(n - 1)) {
      matrix[d_size, d_size] d = rep_matrix(0, d_size, d_size);
      int h = 1;
      
      // Neighbor size for i-th location
      if (i < m) {
        m_i = i;
      } else {
        m_i = m;
      }
  
      // Reconstruct distance matrix (d) for i-th location and its neighbors
      d[1, 2:d_size] = d_pairs[i,];
      d[2:d_size, 1] = to_vector(d_pairs[i,]);
      for (j in 2:m) {
        for (k in (j + 1):d_size) {
          d[k, j] = d_nn_pairs[i, h];
          d[j, k] = d[k, j];
          h = h + 1;
        }
      }
  
      // Covariance function
      C_i = sigma^2 * exp(-d / (2 * l^2));
      // L = cholesky_decompose(C_i);

      A[i + 1, nn[i, 1:m_i]] = C_i[1, 2:(m_i + 1)] * inverse(C_i[2:(m_i + 1), 2:(m_i + 1)]);
      dd[i + 1] = C_i[1, 1] - C_i[1, 2:(m_i + 1)] * A[i + 1, nn[i, 1:m_i]]';
      // Compute U based on A, dd and apply Cholesky
    }
    
    return -0.5 * (sum(log(dd)) + (x' * (I - A)') * diag_matrix(1 ./ dd) * ((I - A) * x));
  }
}

data{
  int<lower = 1> n;
  int<lower = 1> m;
  int<lower = 1> p;
  int<lower = 0> y[n];
  matrix[n, p] Z;
  int nn[n - 1, m];
  matrix[n - 1, m] d_pairs;
  matrix[n - 1, m * (m - 1) / 2] d_nn_pairs;
}

parameters{
  // Betas
  vector[p] beta;
  // Latent field
  vector[n] x;
  // GP standard deviation parameters
  real<lower=0> sigma;
  // GP length-scale parameters
  real<lower=0> l;
}

model{
  matrix[n, n] C;
  vector[n] mu;
  
  x ~ nngp_lpdf(d_pairs, d_nn_pairs, nn, n, m, sigma, l);
  
  // Priors
  l ~ inv_gamma(3, 1);
  l ~ inv_gamma(3, 1);
  beta ~ normal(0, 100);
  
  // Data likelihood
  mu = Z * beta + x;
  y ~ poisson_log(mu);
}
