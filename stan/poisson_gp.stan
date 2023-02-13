/* Latent GP model
https://mc-stan.org/users/documentation/case-studies/nngp.html
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
    vector[N - 1] w_raw;
}

transformed parameters {
  // Hard sum-to-zero constrain
  vector[N] w = append_row(w_raw, -sum(w_raw));
}

model{
  matrix[N, N] K = cov_exp_quad(coords[1:N], sigma, l);

  w ~ multi_normal(rep_vector(0, N), K);
  l ~ inv_gamma(2, 1);
  sigma ~ inv_gamma(2, 1);
  beta ~ normal(0, 1);
  Y ~ poisson_log(X * beta + w);
}
