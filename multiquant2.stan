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
  
  real asym_double_multi_lpdf(vector y, vector mu15, vector mu50, vector mu85, real scale15, real scale50, real scale85) { 
    vector[num_elements(y)] l15;
    vector[num_elements(y)] l50;
    vector[num_elements(y)] l85;
    vector[num_elements(y)] l;
    
    for (i in 1:num_elements(y)){
      l15[i] = log(0.15 * (1 - 0.15)) - log(scale15) - rho_quantile((y[i] - mu15[i]) / scale15, 0.15); 
      l50[i] = log(0.50 * (1 - 0.50)) - log(scale50) - rho_quantile((y[i] - mu50[i]) / scale50, 0.50);  
      l85[i] = log(0.85 * (1 - 0.85)) - log(scale85) - rho_quantile((y[i] - mu85[i]) / scale85, 0.85);  
      l[i] = l15[i] + l50[i] + l85[i];
    }
    return(sum(l));
  }
  
  
  real asym_double_multi_single_lpdf(real y, real mu15, real mu50, real mu85, real scale15, real scale50, real scale85) { 
    real l15;
    real l50;
    real l85;
    real l;
    l15 = log(0.15 * (1 - 0.15)) - log(scale15) - rho_quantile((y - mu15) / scale15, 0.15); 
    l50 = log(0.50 * (1 - 0.50)) - log(scale50) - rho_quantile((y - mu50) / scale50, 0.50);  
    l85 = log(0.85 * (1 - 0.85)) - log(scale85) - rho_quantile((y - mu85) / scale85, 0.85);  
    l = l15 + l50 + l85;
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
  real Intercept15;
  real Intercept50;
  real Intercept85;
  vector[P] beta15;
  vector[P] beta50;
  vector[P] beta85;
  real<lower=0> tau15;
  real<lower=0> tau50;
  real<lower=0> tau85;
  real<lower=0> eta;
}
transformed parameters{
  vector[N] mu15;
  vector[N] mu50;
  vector[N] mu85;
  real scale15 = sqrt(pow(tau15, -1));
  real scale50 = sqrt(pow(tau50, -1));
  real scale85 = sqrt(pow(tau85, -1));
  for (i in 1:N){
    mu15[i] = X[i] * beta15 + Intercept15;
    mu50[i] = X[i] * beta50 + Intercept50;
    mu85[i] = X[i] * beta85 + Intercept85;
  }
}
model{
  real sqrtinveta = sqrt(pow(eta, -1));
  target += gamma_lpdf(eta | df * 0.50, pow(s, 2) * (df * 0.50));
  target += gamma_lpdf(tau15 | 0.01, 0.01);
  target += gamma_lpdf(tau50 | 0.01, 0.01);
  target += gamma_lpdf(tau85 | 0.01, 0.01);
  target += normal_lpdf(Intercept15 | 0, 1e5);
  target += normal_lpdf(Intercept50 | 0, 1e5);
  target += normal_lpdf(Intercept85 | 0, 1e5);
  target += normal_lpdf(beta15 | 0, 1);
  target += normal_lpdf(beta50 | 0, sqrtinveta);
  target += normal_lpdf(beta85 | 0, 1);
  target += asym_double_multi_lpdf(y | mu15, mu50, mu85, scale15, scale50, scale85);
}
generated quantities{
  vector[P] sdbeta;
  real Deviance;
  vector[N] ySim15;
  vector[N] ySim50;
  vector[N] ySim85;
  vector[N] log_lik;
  {
    matrix[P, 3] b;
    
    for (p in 1:P){
      b[p, 1] = beta15[p];
      b[p, 2] = beta50[p];
      b[p, 3] = beta85[p];
    }
    for (p in 1:P){
      sdbeta[p] = sd(b[p, ]);
    }
  }
  Deviance = -2 * asym_double_multi_lpdf(y | mu15, mu50, mu85, scale15, scale50, scale85);
  for(i in 1:N){
    ySim15[i] = asym_double_exp_rng(mu15[i], scale15, 0.15);
    ySim50[i] = asym_double_exp_rng(mu50[i], scale50, 0.50);
    ySim85[i] = asym_double_exp_rng(mu85[i], scale85, 0.85);
  }
  for(i in 1:N){
    log_lik[i] = asym_double_multi_single_lpdf(y[i] | mu15[i], mu50[i], mu85[i], scale15, scale50, scale85);
  }
}
