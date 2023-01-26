/* Latent NNGP model */

functions{
  /* The term inside of exponental of nngp can be written as 
  -1/2 * U^T * D^{-1} * U,
  where U = (I - A) * X and A is a lower triangular matrix and has at max M 
  non-zero entries. Moreover, U can be vectorised and A can be a matrix
  storing non-zero elements. */
  real nngp_w_lpdf(vector w, real sigmasq, real lsq, matrix NN_dist, 
      matrix NN_distM, int[,] NN_ind, int N, int M){
    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;
    
    V[1] = 1;
    
    // For each location i compute the u_i = (I - A) * W
    for (i in 2:N) {
      matrix[i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1) : M] iNNdistM;
      matrix[i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1) : M] iNNCholL;
      vector[i < (M + 1) ? (i - 1) : M] iNNcorr;
      vector[i < (M + 1) ? (i - 1) : M] v;
      row_vector[i < (M + 1) ? (i - 1) : M] v2;
      dim = (i < (M + 1)) ? (i - 1) : M;
      
      // Scaled covariance matrix of neighbors of i-th location - C(c_i, c_i)
      if(dim == 1){
        iNNdistM[1, 1] = 1;
      }
      else{
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            iNNdistM[j, k] = exp(-NN_distM[(i - 1), h] / (2 * lsq));
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for(j in 1:dim){
          iNNdistM[j, j] = 1;
        }
      }
      
      // C(c_i, c_i) = L * L^T
      iNNCholL = cholesky_decompose(iNNdistM);
      
      // Scaled covariance vector between i-th location and its neighbors - C(s_i, c_i)
      iNNcorr = to_vector(exp(-NN_dist[(i - 1), 1:dim] / (2 * lsq)));
      
      // Stan: inverse(tri(A)) * b
      v = mdivide_left_tri_low(iNNCholL, iNNcorr);
      
      // Diagonal elements of D, i.e, d_i
      V[i] = 1 - dot_self(v);
      
      // Compute a_i. Stan: b * inverse(tri(A))
      v2 = mdivide_right_tri_low(v', iNNCholL);

      // u_i = w_i - a_i * w_{c_i}
      I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
    }
    
    // I_Aw ./ V is the element-wise division, i.e, u_i / d_i
    return -0.5 * (1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) + sum(log(V)) + 
      N * log(sigmasq));
  }
}

data {
    int<lower = 1> N;
    int<lower = 1> M;
    int<lower = 1> P;
    int<lower = 0> Y[N];
    matrix[N, P] X;
    int NN_ind[N - 1, M];
    matrix[N - 1, M] NN_dist;
    matrix[N - 1, (M * (M - 1) / 2)] NN_distM;
}

parameters{
    vector[P] beta;
    real<lower = 0> sigma;
    real<lower = 0> l;
    vector[N] w;
}

transformed parameters {
    real sigmasq = square(sigma);
    real lsq = square(l);
}

model{
  l ~ inv_gamma(3, 1);
  sigma ~ inv_gamma(3, 1);
  beta ~ normal(0, 100);
  w ~ nngp_w(sigmasq, lsq, NN_dist, NN_distM, NN_ind, N, M);
  Y ~ poisson_log(X * beta + w);
}

