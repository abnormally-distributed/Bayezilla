functions {

  real rho_quantile(real y, real quantile) {
     if (y < 0) {
       return y * (quantile - 1);
     } else {
       return y * quantile;
     }
   }
  real asym_double_exp_lpdf(vector y, vector mu, real scale, real quantile) { 
    vector[num_elements(y)] l;
    for (i in 1:num_elements(y)){
      l[i] = log(quantile * (1 - quantile)) - log(scale) - rho_quantile((y[i] - mu[i]) / scale, quantile); 
    }
   return(sum(l));
   }
   
   real asym_double_multi_lpdf(vector y, vector mu10, vector mu25, vector mu50, vector mu75, vector mu90, real scale10, real scale25, real scale50, real scale75, real scale90) { 
    vector[num_elements(y)] l10;
    vector[num_elements(y)] l25;
    vector[num_elements(y)] l50;
    vector[num_elements(y)] l75;
    vector[num_elements(y)] l90;
    vector[num_elements(y)] l;
    
    for (i in 1:num_elements(y)){
      l10[i] = log(0.10 * (1 - 0.10)) - log(scale10) - rho_quantile((y[i] - mu10[i]) / scale10, 0.10); 
      l25[i] = log(0.25 * (1 - 0.25)) - log(scale25) - rho_quantile((y[i] - mu25[i]) / scale25, 0.25);  
      l50[i] = log(0.50 * (1 - 0.50)) - log(scale50) - rho_quantile((y[i] - mu50[i]) / scale50, 0.50);  
      l75[i] = log(0.75 * (1 - 0.75)) - log(scale75) - rho_quantile((y[i] - mu75[i]) / scale75, 0.75);  
      l90[i] = log(0.90 * (1 - 0.90)) - log(scale90) - rho_quantile((y[i] - mu90[i]) / scale90, 0.90);  
      l[i] = l10[i] + l25[i] + l50[i] + l75[i] + l90[i];
    }
    return(sum(l));
    }
   
   
  real asym_double_multi_single_lpdf(real y, real mu10, real mu25, real mu50, real mu75, real mu90, real scale10, real scale25, real scale50, real scale75, real scale90) { 
    real l10;
    real l25;
    real l50;
    real l75;
    real l90;
    real l;
      l10 = log(0.10 * (1 - 0.10)) - log(scale10) - rho_quantile((y - mu10) / scale10, 0.10); 
      l25 = log(0.25 * (1 - 0.25)) - log(scale25) - rho_quantile((y - mu25) / scale25, 0.25);  
      l50 = log(0.50 * (1 - 0.50)) - log(scale50) - rho_quantile((y - mu50) / scale50, 0.50);  
      l75 = log(0.75 * (1 - 0.75)) - log(scale75) - rho_quantile((y - mu75) / scale75, 0.75);  
      l90 = log(0.90 * (1 - 0.90)) - log(scale90) - rho_quantile((y - mu90) / scale90, 0.90);  
      l = l10 + l25 + l50 + l75 + l90;
    return(l);
    }
   
  real asym_double_exp_single_lpdf(real y, real mu, real scale, real quantile) { 
     return log(quantile * (1 - quantile)) - log(scale) - rho_quantile((y - mu) / scale, quantile); 
   }

  real asym_double_exp_rng(real mu, real scale, real quantile){
        real q;
        real p = beta_rng(1, 1);
        if (p < quantile){
          q =  mu + ((scale * log(p/quantile))/(1 - quantile));
        } 
        if (p > quantile){
          q = mu - ((scale * log((1 - p)/(1 - quantile)))/quantile);
        }
        return(q);
   }
}
data {
  int N; //the number of observations
  int P; //the number of columns in the model matrix
  vector[N] y; //the response
  matrix[N,P] X; //the model matrix
  real<lower=0> s;
  real<lower=1> df;
}
parameters {
  real Intercept10;
  real Intercept25;
  real Intercept50;
  real Intercept75;
  real Intercept90;
  vector[P] beta10;
  vector[P] beta25;
  vector[P] beta50;
  vector[P] beta75;
  vector[P] beta90;
  real<lower=0> tau10;
  real<lower=0> tau25;
  real<lower=0> tau50;
  real<lower=0> tau75;
  real<lower=0> tau90;
  real<lower=0> eta;
}
transformed parameters{
  vector[N] mu10;
  vector[N] mu25;
  vector[N] mu50;
  vector[N] mu75;
  vector[N] mu90;
  real scale10 = sqrt(pow(tau10, -1));
  real scale25 = sqrt(pow(tau25, -1));
  real scale50 = sqrt(pow(tau50, -1));
  real scale75 = sqrt(pow(tau75, -1));
  real scale90 = sqrt(pow(tau90, -1));
    for (i in 1:N){
    mu10[i] = X[i] * beta10 + Intercept10;
    mu25[i] = X[i] * beta25 + Intercept25;
    mu50[i] = X[i] * beta50 + Intercept50;
    mu75[i] = X[i] * beta75 + Intercept75;
    mu90[i] = X[i] * beta90 + Intercept90;
  }
}
model{
  real sqrtinveta = sqrt(pow(eta, -1));
  target += gamma_lpdf(eta | df * 0.50, pow(s, 2) * (df * 0.50));
  target += gamma_lpdf(tau10 | 0.01, 0.01);
  target += gamma_lpdf(tau25 | 0.01, 0.01);
  target += gamma_lpdf(tau50 | 0.01, 0.01);
  target += gamma_lpdf(tau75 | 0.01, 0.01);
  target += gamma_lpdf(tau90 | 0.01, 0.01);
  target += normal_lpdf(Intercept10 | 0, 1e5);
  target += normal_lpdf(Intercept25 | 0, 1e5);
  target += normal_lpdf(Intercept50 | 0, 1e5);
  target += normal_lpdf(Intercept75 | 0, 1e5);
  target += normal_lpdf(Intercept90 | 0, 1e5);
  target += normal_lpdf(beta10 | 0, sqrtinveta);
  target += normal_lpdf(beta25 | 0, sqrtinveta);
  target += normal_lpdf(beta50 | 0, sqrtinveta);
  target += normal_lpdf(beta75 | 0, sqrtinveta);
  target += normal_lpdf(beta90 | 0, sqrtinveta);
  target += asym_double_multi_lpdf(y | mu10, mu25, mu50, mu75, mu90, scale10, scale25, scale50, scale75, scale90);
}
generated quantities{
  vector[P] sdbeta;
  real Deviance;
  vector[N] ySim10;
  vector[N] ySim25;
  vector[N] ySim50;
  vector[N] ySim75;
  vector[N] ySim90;
  vector[N] log_lik;
  {
  matrix[P, 5] b;
  
  for (p in 1:P){
    b[p, 1] = beta10[p];
    b[p, 2] = beta25[p];
    b[p, 3] = beta50[p];
    b[p, 4] = beta75[p];
    b[p, 5] = beta90[p];
  }
  for (p in 1:P){
    sdbeta[p] = sd(b[p, ]);
  }
  }
  Deviance = -2 * asym_double_multi_lpdf(y | mu10, mu25, mu50, mu75, mu90, scale10, scale25, scale50, scale75, scale90);
  for(i in 1:N){
    ySim10[i] = asym_double_exp_rng(mu10[i], scale10, 0.10);
    ySim25[i] = asym_double_exp_rng(mu25[i], scale25, 0.25);
    ySim50[i] = asym_double_exp_rng(mu50[i], scale50, 0.50);
    ySim75[i] = asym_double_exp_rng(mu75[i], scale75, 0.75);
    ySim90[i] = asym_double_exp_rng(mu90[i], scale90, 0.90);
  }
  for(i in 1:N){
    log_lik[i] = asym_double_multi_single_lpdf(y[i] | mu10[i], mu25[i], mu50[i], mu75[i], mu90[i], scale10, scale25, scale50, scale75, scale90);
  }
}
