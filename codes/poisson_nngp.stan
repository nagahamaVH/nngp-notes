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
  
  matrix nngp_chol(
    matrix d_pairs, matrix d_nn_pairs, int[,] nn, int m, real sigma, real l) {
    int n = dims(d_pairs)[1] + 1;
    int m_i;
    int d_size = m + 1;
    matrix[n, n] C;
    matrix[n, n] I = diag_matrix(rep_vector(1, n));
    matrix[d_size, d_size] d = rep_matrix(0, d_size, d_size);
    matrix[d_size, d_size] C_i;
    // Initialize A = [0] and dd = diag(D)
    matrix[n, n] A = rep_matrix(0, n, n);
    vector[n] dd;
    dd[1] = sigma^2;

    // Update A and d
    for (i in 1:(n - 1)) {
      // Reconstruct distance matrix (d) for point i and its neighbours
      d[2:d_size, 1] = to_vector(d_pairs[i,]);
      d[1, 2:d_size] = d_pairs[i,];
      for (j in 2:m) {
        for (k in (j + 1):d_size) {
          d[k, j] = d_nn_pairs[i, k - j];
          d[j, k] = d[k, j];
        }
      }

      // Avoid selecting 0s in nn array
      if (i < m) {
        m_i = i;
      } else {
        m_i = m;
      }
      
      C_i = sigma^2 * exp(-d / (2 * l^2));
      A[i + 1, nn[i, 1:m_i]] = to_row_vector(
        C_i[1, 2:(m_i + 1)] * inverse(C_i[2:(m_i + 1), 2:(m_i + 1)]));
      dd[i + 1] = C_i[1, 1] - A[i + 1, nn[i, 1:m_i]] * C_i[1, 2:(m_i + 1)];
    }

    A = inverse(I - A);
    C = A * diag_matrix(dd) * A';
    return C;
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
  matrix[n - 1, m] d_nn_pairs;
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
  
  C = nngp_chol(d_pairs, d_nn_pairs, nn, m, sigma, l);
  x ~ multi_normal(rep_vector(0, n), C);
  
  // Priors
  l ~ inv_gamma(3, 0.5);
  sigma ~ normal(0, 3);
  beta ~ normal(0, 1);
  
  // Data likelihood
  mu = Z * beta + x;
  y ~ poisson_log(mu);
}
