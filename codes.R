one_sample_poisson <- "model {
  y ~ dpois(lambda * t)
  lambda ~ dgamma(0.001, 0.001)
  ySim~ dpois(lambda * t)
}"

two_sample_poisson<- "model {
  for(group_i in 1:2) {
    y[group_i] ~ dpois(lambda[group_i] * t[group_i])
    lambda[group_i] ~ dgamma(0.001, 0.001)
    ySim[group_i] ~ dpois(lambda[group_i] * t[group_i])
  }
  rate_diff <- lambda[1] - lambda[2]
  rate_ratio <- lambda[1] / lambda[2]
}"

two_sample_binom_model <- "model {
  for(i in 1:length(x)) {
    y[i] ~ dbinom(phi[i], n[i])
    phi[i] ~ dbeta(a, b)
    ySim[i] ~ dbinom(phi[i], n[i])
  }
  phi_diff <- phi[1] - phi[2]
}"

binom_model <- "
data{
  compval <- a / (a+b)
}
model {
  y ~ dbinom(phi, n)
  phi ~ dbeta(a, b)
  ySim ~ dbinom(phi, n)
  phi_diff <- phi - compval
}"
