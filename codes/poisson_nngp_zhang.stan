/* Latent NNGP model*/
  
functions{
  real nngp_w_lpdf(vector w, real sigmasq, real l, matrix NN_dist,
                   matrix NN_distM, int[,] NN_ind, int N, int M){
    
    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;
    
    for (i in 2:N) {
      
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNdistM;
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNCholL;
      vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
      vector[ i < (M + 1)? (i - 1) : M] v;
      row_vector[i < (M + 1)? (i - 1) : M] v2;
      
      dim = (i < (M + 1))? (i - 1) : M;
      
      if(dim == 1){iNNdistM[1, 1] = 1;}
      else{
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            iNNdistM[j, k] = exp(-NN_distM[(i - 1), h] / (2 * l^2));
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for(j in 1:dim){
          iNNdistM[j, j] = 1;
        }
      }
      
      iNNCholL = cholesky_decompose(iNNdistM);
      iNNcorr = to_vector(exp(-NN_dist[(i - 1), 1:dim] / (2 * l^2)));
      
      v = mdivide_left_tri_low(iNNCholL, iNNcorr);
      
      V[i] = 1 - dot_self(v);
      
      v2 = mdivide_right_tri_low(v', iNNCholL);

            I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];

        }
        V[1] = 1;
        return - 0.5 * ( 1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                        sum(log(V)) + N * log(sigmasq));
    }
}

data {
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> P;
    int<lower = 0> Y[N];
    matrix[N, P] X;
    int NN_ind[N - 1, M];
    matrix[N - 1, M] NN_dist;
    matrix[N - 1, M] NN_distM;
}

parameters{
    vector[P] beta;
    real<lower = 0> sigma;
    real<lower = 0> l;
    vector[N] w;
}

transformed parameters {
    real sigmasq = square(sigma);
}

model{
  vector[N] mu;

  l ~ inv_gamma(3, 0.5);
  sigma ~ normal(0, 3);
  beta ~ normal(0, 1);
  w ~ nngp_w(sigmasq, l, NN_dist, NN_distM, NN_ind, N, M);
  mu = X * beta + w;
  Y ~ poisson_log(mu);
}

