functions{
  matrix se_kernel(matrix D, real sigma, real l) {
    int n = dims(D)[1];
    matrix[n, n] C;
    
    for (i in 1:(n - 1)) {
      C[i, i] = sigma^2;
      for (j in (i + 1):n) {
        C[i, j] = sigma^2 * exp(-D[i, j] / (2 * l^2));
        C[j, i] = C[i, j];
      }
    }
    C[n, n] = sigma^2;
    return C;
  }
}

data{
  int<lower = 1> n;
  int<lower = 1> p;
  int<lower = 0> y[n];
  matrix[n, n] D;
  matrix[n, p] Z;
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
  
  C = se_kernel(D, sigma, l);
  x ~ multi_normal(rep_vector(0, n), C);
  
  // Priors
  l ~ inv_gamma(3, 0.5);
  sigma ~ normal(0, 3);
  beta ~ normal(0, 1);

  // Data likelihood
  mu = Z * beta + x;
  y ~ poisson_log(mu);
}
