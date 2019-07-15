data {
  int P;
  int<lower=0> N;
  matrix[N, P] X;
  matrix[1, N] y;
}
parameters {
  vector<lower=0>[N] w;
  real<lower=0> alpha; 
  real<lower=0> rho;
}
transformed parameters{
 matrix[N, N] Sigma = cov_exp_quad(X, alpha, rho);
}
model {
  target += gamma_lpdf(alpha | 2, .01);
  target += gamma_lpdf(rho | 2, 1);
  target += multi_gp_lpdf(matrix y | Sigma, w);
}
generated quantities{
  vector[N] ySim;
  {
  vector[N] invw;
  for (i in 1:N){
    invw[i] = pow(w[i],-1);
  }
  ySim = multi_normal_rng(rep(0, N), invw, Sigma);
  }
}