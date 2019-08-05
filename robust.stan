data {
  int<lower=1> K;     // number of mixture components
  int<lower=1> N;         // number of data points
  int<lower=1> P;         // number of predictors
  real<lower=0, upper=1> pi;
  real y[N];              // observations
  matrix[N,P] X;
}
parameters {
  real Intercept;
  vector[P] beta;
  // scales of mixture components
  real<lower=0> tau_1;    
  real<lower=0> tau_2; 
}
transformed parameters{
  real sigma_1 = pow(tau_1, -0.50);
  real sigma_2 = pow(tau_2, -0.50);
}
model {
  target += student_t_lpdf(Intercept | 1, 0, 1000);
  target += normal_lpdf(beta | 0, 1);
  target += gamma_lpdf(tau_1 | 0.01, 0.01);
  target += gamma_lpdf(tau_2 | 0.01, 0.01);
  for (n in 1:N)
      target += log_mix(pi, normal_lpdf(y[n] | Intercept + X*beta, sigma_1), double_exponential_lpdf(y[n] | Intercept + X*beta, sigma_2));
}
