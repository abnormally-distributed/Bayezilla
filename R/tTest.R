#' Bayesian 't-tests' 
#'
#' @description 
#' 
#' 
#' INDEPENDENT SAMPLES T-TEST: 
#' 
#' The independent samples t-test places cauchy priors on the means of both groups with a median
#' set at the overall mean of the measurement (null hypothesis = same means). The posterior distribution of the group
#' differences is taken as the difference of each posterior group mean. The likelihood function is a student t distribution
#' with a gamma(2, .01) prior on the normality parameter and a scaled gamma prior on the precisions of each group, implying
#' half-cauchy priors on the standard deviations of each group. Note the use of a student-t likelihood is wholly unrelated to 
#' this being analagous to the frequentist t-test. In the t-test, the student-t distribution is the asymptotic sampling distribution of the 
#' t-statistic, but here it is used to model the heaviness of the tails in the data. When there are no outliers the resulting estimations and inferences
#' remain Gaussian. This is an adaptive prior and not one that could prove to be a disadvantage. For this reason the Gaussian likelihood is not offered here.
#'
#' The Bayes Factor is estimated as the exponentiated difference in log density of the effect evaluated at zero and the log density evaluated
#' at each point within the posterior distribution (the Savage-Dickey ratio method). Hence, this is an estimator of the log-marginal likelihood function of the group differences.
#' Use the median Bayes Factor if you use the median as your point estimate, and the mean if using the posterior mean as your point estimate. 
#' 
#' Several effect size measures are also provided. A quantity variously known as Cohen's d or Hedge's g (1981)* is given as the primary effect size measure.
#' Due to the ambiguity in which term is correct, it is simply labeled as "effectSize" in the output. The specific formula used is given below: 
#' 
#' \deqn{\sqrt{\frac{\left[\sigma_{1}^{2} \times\left(n_{1}-1\right)\right]+\left[\sigma_{2}^{2} \times\left(n_{2}-1\right)\right]}{n_{1}+n_{2}-2}}}
#' 
#' This formula for the effect size is robust to differences in the variances across groups.
#' 
#' Cohen's U3 (1977) is defined as a measure of non-overlap, which quantifies the percentage of data points
#' in group A that are smaller than the median of group B. It is given by the probability density below the observed effect size via the normal cumulative
#' distribution function; phi(effectSize). 
#' 
#' A similar measure of effect size known as the probability of superiority, or alternatively the common language effect size (CLES), is given by the effect size
#' divided by the square root of two entered into the normal cumulative distribution function; phi(effectSize / sqrt(2)). It is very similar in size to Cohen's U3 and 
#' carries a similar meaning, but is more strightforward. The CLES gives the probability that an observation sampled at random from group A will have a higher value 
#' than an observation sampled at random from group B. The common language effect size was developed to facilitate an easy way to explain the importance of a statistical 
#' result to the layman who may not have an intuition for Cohen's d or Hedge's g, which are defined as changes in standard deviation. It is labeled "CL" in the posterior output. 
#'
#' Finally, the output contains posterior predictive simulations for each participant (column-wise). The row-wise posterior predictive distribution gives simulated samples of size N 
#' that allow you to visualize what future data sets of the same size might look like.   
#'         
#' *The labels Cohen's d and Hedge's g have both been applied to various measures of effect size. See the citation below for a discussion of this.
#'         
#'         Difference between Cohen's d and Hedges' g for effect size metrics, URL (version: 2018-09-28): https://stats.stackexchange.com/q/338043
#'
#'
#' REPEATED MEASURES T-TEST:
#' 
#' Much of the information is the same as above, except Cohen's U3 is not provided as it lacks an intuitive definition for this case. 
#' 
#' The effect size is simply the standardized mean difference (z-score) 
#' 
#' ONE SAMPLE T-TEST:
#' 
#' Same as above, but the CLES is also dropped because it is not applicable to one group. The effect size is calculated as the mean of y - the comparison value  / sd(y)
#'  
#' @param formula If using an independent samples t-test, supply a formula of the form y ~ group 
#' @param data the data frame containing the outcome variable and group labels
#' @param y if using the one sample t-test, the vector of observations. If using repeated measures, the vector of group differences.
#' @param model one of "is" (independent samples t-test) , "rm" (repeated measures t-test), or "os" (one sample t-test)
#' @param compval the hypothesized null value for a one sample t-test
#' @param iter 
#' @param warmup 
#' @param adapt 
#' @param chains 
#' @param thin 
#' @param method 
#' @param cl 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 

tTest = function(formula, data, y = NULL, model = "is", iter=10000, warmup=1000, adapt=5000, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){
  
  RNGlist = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")
  
  if (chains > 4){
    chains = 4
  }
  
  if (model == "is"){
  
  two_sample_t <- 
    
  "model {

  nu ~ dgamma(2, .01) T(1, )

  for (g in 1:2){
    tau[g] ~ dscaled.gamma(ysd, 1)
    mu[g]  ~ dt(ymean , .5 * 1/sqrt(ysd), 1)
  }    
    
  for (i in 1:N){
       y[i] ~ dt(mu[group[i]], tau[group[i]], nu)
       ySim[i] ~ dt(mu[group[i]], tau[group[i]], nu)
  }

  sigma[1] <- sqrt(1 / tau[1])
  sigma[2] <- sqrt(1 / tau[2])
  mu_diff <- mu[1] - mu[2]
  sigma_diff <- sigma[1] - sigma[2] 
  effectSize <- (mu_diff) / sqrt( ( (pow(sigma[1],2)*(N1-1)) + (pow(sigma[2],2)*(N2-1)) ) / (N1+N2-2) )
  U3 <- phi(effectSize)
  CL <- phi(effectSize / sqrt(2))
  prior_effectSize <- logdensity.norm(0, 0, 1)
  posterior_effectSize <- logdensity.norm( 0, effectSize, 1)
  BF <- exp(prior_effectSize - posterior_effectSize)
}
"

group = as.numeric(as.factor(model.matrix(formula, data)[,2]))
N = length(y)
N1 = as.vector(table(group))[1]
N2 = as.vector(table(group))[2]
y = model.frame(formula, data)[,1]
write_lines(two_sample_t, "two_sample_t.txt")
jagsdata = list("N" = N, "N1" = N1, "N2" = N2, "group" = group, "y" = y, "ysd" = mad(y), "ymean" = median(y))
monitor = c("mu_diff", "sigma_diff", "mu", "sigma", "nu", "effectSize", "U3", "CL", "BF", "ySim")
inits = lapply(1:chains, function(z) list("mu" = c(median(y), median(y)), "ySim"  = y, "tau" = c(1/var(y),1/var(y)), "nu" = 3,
                                          .RNG.name=RNGlist[z], .RNG.seed= sample(1:10000, 1)))
out = run.jags(model = "two_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
if (!is.null(cl)){
  parallel::stopCluster(cl = cl)
}
file.remove("two_sample_t.txt")
return(out)
  
}

if (model == "rm") {
  
  repeated_measures_t <- "
  
  model {
    nu ~ dgamma(2, .01) T(1, )
    mu ~ dt(0, 1, 1)
    tau ~ dscaled.gamma(ysd, 1)
    
    for (i in 1:N) {
      y[i] ~ dt(mu, tau , nu)
      ySim[i] ~ dt(mu, tau, nu)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - 0) / sigma
    CL <- phi(effectSize / sqrt(2))
    prior_effectSize <- logdensity.norm(0, 0, 1)
    posterior_effectSize <- logdensity.norm(0, effectSize, 1)
    BF <- exp(prior_effectSize - posterior_effectSize)
}"
 
  if (is.null(y)) {
    cat(crayon::red("Please provide a vector of paired differences "))
  }
  
  write_lines(repeated_measures_t, "repeated_measures_t.txt")
  jagsdata = list("N" = length(y), "y" = y, "ysd" = mad(y))
  monitor = c("mu", "sigma", "nu", "effectSize", "CL", "BF", "ySim")
  inits = lapply(1:chains, function(z) list("mu" = median(y), "ySim"  = y, "tau" = 1/var(y), "nu" = 3,
                                            .RNG.name=RNGlist[z], .RNG.seed= sample(1:10000, 1)))
  out = run.jags(model = "repeated_measures_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("repeated_measures_t.txt")
  return(out)
  
}

if (model == "os") {
  
  one_sample_t <- "
  
  model {
    nu ~ dgamma(2, .01) T(1, )
    mu ~ dt(0, 1, 1)
    tau ~ dscaled.gamma(ysd, 1)
    
    for (i in 1:N) {
      y[i] ~ dt(mu, tau , nu)
      ySim[i] ~ dt(mu, tau, nu)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - compval) / sigma
    prior_effectSize <- logdensity.norm(compval, compval, 1)
    posterior_effectSize <- logdensity.norm(compval, effectSize, 1)
    BF <- exp(prior_effectSize - posterior_effectSize)
}"
  
  if (is.null(y)) {
    cat(crayon::red("Please provide a vector of observations"))
  }
  if (is.null(compval)) {
    cat(crayon::red("Please provide a comparison value against which to test the mean of y"))
  }
  
  write_lines(one_sample_t, "one_sample_t.txt")
  jagsdata = list("N" = length(y), "y" = y, "ysd" = mad(y))
  monitor = c("mu", "sigma", "nu", "effectSize", "BF", "ySim")
  
  inits = lapply(1:chains, function(z) list("mu" = median(y), "ySim"  = y, "tau" = 1/var(y), "nu" = 3,
                                            .RNG.name=RNGlist[z], .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "one_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, 
                 burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, ...)
  
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  
  file.remove("one_sample_t.txt")
  return(out)

  }
}

