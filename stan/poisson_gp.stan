/* Latent GP model
https://mc-stan.org/users/documentation/case-studies/nngp.html
https://github.com/mbjoseph/gpp-speed-test/blob/master/stan/full_pois.stan
*/

data {
    int<lower = 1> N;
    int<lower = 1> P;
    int<lower = 0> Y[N];
    matrix[N, P] X;
    row_vector[2] coords[N];
}

parameters{
    vector[P] beta;
    real<lower = 0> sigma;
    real<lower = 0> l;
    vector[N] z;
}

transformed parameters {
  matrix[N, N] K;
  vector[N] w;
  K = cov_exp_quad(coords, sigma, l);
  for (i in 1:N) {
    K[i, i] += 1e-12;
  }
  w = cholesky_decompose(K) * z;
}

model{
  z ~ normal(0, 1);
  l ~ inv_gamma(2, 1);
  sigma ~ inv_gamma(2, 1);
  beta ~ normal(0, 1);
  Y ~ poisson_log(X * beta + w);
}
