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
  
  matrix nngp_chol(matrix D, matrix D_m, matrix nn, int m, real sigma, real l) {
    int n = dims(nn)[1] + 1;
    matrix[m + 1, m + 1] D_i;
    matrix[m + 1, m + 1] C_i;
    matrix[n, n] C;
    matrix[n, n] I;

    // Initialize A = [0] and dd = diag(D)
    matrix[n, n] A = rep_matrix(0, n, n);
    vector[n] dd;
    dd[1] = sigma^2;

    // Update A and dd
    for (i in 1:(n - 1)) {
      // Reconstruct  distance matrix for point i and its neighbours
      matrix[n, n] DD = rep_matrix(0, n, n);

      DD[2:n, 1] = to_vector(D[i,]);
      DD[1, 2:n] = DD[1, 2:n];
      for (j in 2:m) {
        for (k in (j + 1):n) {
          DD[j, k] = D_m[i, k - j];
          DD[k, j] = DD[j, k];
        }
      }
      C_i = sigma^2 * exp(-DD / (2 * l^2));
      A[i + 1, 1:m] = to_row_vector(C_i[1, -1] * inverse(C_i[2:m, 2:m]));
      dd[i + 1] = C_i[1, 1] - A[i + 1, nn[i,]] * C_i[1, -1];
    }

    I = diag_matrix(rep_vector(1, n));
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
  matrix[n - 1, m] nn;
  matrix[n - 1, m] D;
  matrix[n - 1, m] D_m;
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
  
  C = nngp_chol(D, D_m, nn, m, sigma, l);
  x ~ multi_normal(rep_vector(0, n), C);
  
  // Priors
  l ~ inv_gamma(3, 0.5);
  sigma ~ normal(0, 3);
  beta ~ normal(0, 1);
  
  // Data likelihood
  mu = Z * beta + x;
  y ~ poisson_log(mu);
}
