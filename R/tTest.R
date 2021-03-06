#' Bayesian 't-tests' 
#'
#' @description 
#' INDEPENDENT SAMPLES T-TEST: 
#' 
#' The difference in means is modeled as mu_1 - mu_2. The difference in the standard deviations of both groups
#' is modeled as sqrt(1/tau_1 - 1/tau_2).
#' 
#' The log-Bayes Factor is estimated as the  difference in log density of the effect evaluated at zero and the log density evaluated
#' at each point within the posterior distribution (the Savage-Dickey ratio method). It is presented on the log-density scale for
#' numerical stability. Note that the log-Bayes Factor is here presented as the log-evidence in favor of the alternative hypothesis.
#' Hence, a larger log-BF is evidence in favor of a non-zero difference. Use the median Bayes Factor if you use the median as your point estimate, and the mean if using the posterior mean as your point estimate. 
#' 
#' Several effect size measures are also provided. A quantity variously known as Cohen's d or Hedge's g (1981)* is given as the primary effect size measure.
#' Due to the ambiguity in which term is correct, it is simply labeled as "effSize" in the output. The specific formula used is given below: 
#' 
#' \deqn{\sqrt{\frac{\left[\sigma_{1}^{2} \times\left(n_{1}-1\right)\right]+\left[\sigma_{2}^{2} \times\left(n_{2}-1\right)\right]}{n_{1}+n_{2}-2}}}
#' 
#' This formula for the effect size is robust to differences in the variances across groups.
#' 
#' Cohen's U3 (1977) is defined as a measure of non-overlap, which quantifies the percentage of data points
#' in group A that are smaller than the median of group B. It is given by the probability density below the observed effect size via the normal cumulative
#' distribution function; phi(effSize). 
#' 
#' A similar measure of effect size known as the probability of superiority, or alternatively the common language effect size (CLES), is given by the effect size
#' divided by the square root of two entered into the normal cumulative distribution function; phi(effSize / sqrt(2)). It is very similar in size to Cohen's U3 and 
#' carries a similar meaning, but is more strightforward. The CLES gives the probability that an observation sampled at random from group A will have a higher value 
#' than an observation sampled at random from group B. The common language effect size was developed to facilitate an easy way to explain the importance of a statistical 
#' result to the layman who may not have an intuition for Cohen's d or Hedge's g, which are defined as changes in standard deviation. It is labeled "CL" in the posterior output. 
#'
#' Finally, the output contains posterior predictive simulations for each participant (column-wise). The row-wise posterior predictive distribution gives simulated samples of size N 
#' that allow you to visualize what future data sets of the same size might look like.   
#'         
#' *The labels Cohen's d and Hedge's g have both been applied to various measures of effect size. See the citation below for a discussion of this.
#'        
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
#'        
#' @param formula If using an independent samples t-test, supply a formula of the form y ~ group 
#' @param data the data frame containing the outcome variable and group labels for independent samples t-test. 
#' If using the one sample t-test, the vector of observations. 
#' If using repeated measures, the vector of group differences.
#' @param model one of "is" (independent samples t-test), "rm" (repeated measures t-test), or "os" (one sample t-test)
#' @param median change to TRUE to compare the medians instead of means
#' @param compval the hypothesized null value for a one sample t-test
#' @param iter the number of iterations. defaults to 10000.
#' @param warmup number of burnin samples. defaults to 2500.
#' @param adapt number of adaptation steps. defaults to 2500.
#' @param chains number of chains. defaults to 4.
#' @param thin the thinning interval. defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... other arguments to run.jags
#'
#' @return
#' a runjags object
#' @export
#'
#' @examples
#' tTest(len ~ supp, ToothGrowth)
#' 
tTest = function(formula = NULL, data, compval = 0, model = "is", median = FALSE, iter=10000, warmup=2500, adapt=2500, chains=4, thin=3, method = "parallel", cl = makeCluster(2), ...){
 
  if (median == FALSE){

    if (model == "is"){
      
      two_sample_t <- 
        
        "model {

    for (g in 1:2){
      tau[g] ~ dgamma(1e-100, 1e-100)
      mu[g]  ~ dnorm(0, 1e-100)
      sigma[g] <- sqrt(1 / tau[g])
  }    
    
    for (i in 1:N){
      y[i] ~ dnorm(mu[group[i]], tau[group[i]])
      ySim[i] ~ dnorm(mu[group[i]], tau[group[i]])
  }

  muDiff <- mu[1] - mu[2]
  sigmaDiff <- sigma[1] - sigma[2]
  effSize <- (muDiff) / sqrt( ( (pow(sigma[1],2)*(N1-1)) + (pow(sigma[2],2)*(N2-1)) ) / (N1+N2-2) )
  U3 <- phi(effSize)
  CL <- phi(effSize / sqrt(2))
  prior_effSize <- logdensity.norm(0, 0, 1)
  posterior_effSize <- logdensity.norm(0, effSize, 1)
  logBF <- prior_effSize - posterior_effSize
}
"
    y = model.frame(formula, data)[,1]
    group = as.numeric(as.factor(model.matrix(formula, data)[,2]))
    N = length(y)
    N1 = as.vector(table(group))[1]
    N2 = as.vector(table(group))[2]
    write_lines(two_sample_t, "two_sample_t.txt")
    jagsdata = list("N" = N, "N1" = N1, "N2" = N2, "group" = group, "y" = y)
    monitor = c("muDiff", "sigmaDiff", "mu", "sigma",  "effSize", "U3", "CL", "logBF", "ySim")
    inits = lapply(1:chains, function(z) list("mu" = c(mean(y), mean(y)), "ySim"  = y, "tau" = c(1/var(y),1/var(y)),
                                              .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:1000, 1)))
    out = run.jags(model = "two_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    if (!is.null(cl)){
      parallel::stopCluster(cl = cl)
    }
    file.remove("two_sample_t.txt")
    return(out)
    
  }
  
  if (model == "rm") {
    
    repeated_measures_t <- "
  
  model {
    mu ~ dnorm(0, 1e-100)
    tau ~ dgamma(1e-100, 1e-100)
    
    for (i in 1:N) {
      y[i] ~ dnorm(mu, tau)
      ySim[i] ~ dnorm(mu, tau)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - 0) / sigma
    CL <- phi(effSize / sqrt(2))
    prior_effSize <- logdensity.norm(0, 0, 1)
    posterior_effSize <- logdensity.norm(0, effSize, 1)
    logBF <- prior_effSize - posterior_effSize
  }"
    
    if (is.null(data)) {
      cat(crayon::red("Please provide a vector of paired differences "))
    }
    
    y = data
    write_lines(repeated_measures_t, "repeated_measures_t.txt")
    jagsdata = list("N" = length(y), "y" = y)
    monitor = c("mu", "sigma",  "effSize", "CL", "logBF", "ySim")
    inits = lapply(1:chains, function(z) list("mu" = mean(y), "ySim"  = y, "tau" = 1/var(y),
                                              .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
    out = run.jags(model = "repeated_measures_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    if (!is.null(cl)){
      parallel::stopCluster(cl = cl)
    }
    file.remove("repeated_measures_t.txt")
    return(out)
    
  }
  
  if (model == "os") {
    
    one_sample_t <- "
  
  model {
    mu ~ dnorm(0, 1e-100)
    tau ~ dgamma(1e-100, 1e-100)
    
    for (i in 1:N) {
      y[i] ~ dnorm(mu, tau)
      ySim[i] ~ dnorm(mu, tau)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - compval) / sigma
    prior_effSize <- logdensity.norm(compval, compval, 1)
    posterior_effSize <- logdensity.norm(compval, effSize, 1)
    logBF <- prior_effSize - posterior_effSize
}"
    
    if (is.null(data)) {
      cat(crayon::red("Please provide a vector of observations"))
    }
    if (is.null(compval)) {
      cat(crayon::red("Please provide a comparison value against which to test the mean of y"))
    }
    
    y = data
    write_lines(one_sample_t, "one_sample_t.txt")
    jagsdata = list("N" = length(y), "y" = y, compval = compval)
    monitor = c("mu", "sigma",  "effSize", "logBF", "ySim")
    
    inits = lapply(1:chains, function(z) list("mu" = mean(y), "ySim"  = y, "tau" = 1/var(y),
                                              .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
    
    out = run.jags(model = "one_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, 
                   burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    
    if (!is.null(cl)){
      parallel::stopCluster(cl = cl)
    }
    
    file.remove("one_sample_t.txt")
    return(out)
    
    }
  
  }

  if (median == TRUE){
    
    if (model == "is"){
      
      two_sample_t <- 
        
        "model {

    for (g in 1:2){
      tau[g] ~ dscaled.gamma(1e100, 1)
      mu[g]  ~ dnorm(0, 1e-100)
      sigma[g] <- sqrt(1 / tau[g])
    }    
    
    for (i in 1:N){
      y[i] ~ ddexp(mu[group[i]], tau[group[i]])
      ySim[i] ~ ddexp(mu[group[i]], tau[group[i]])
    }

  muDiff <- mu[1] - mu[2]
  sigmaDiff <- sigma[1] - sigma[2]
  effSize <- (muDiff) / sqrt( ( (pow(sigma[1],2)*(N1-1)) + (pow(sigma[2],2)*(N2-1)) ) / (N1+N2-2) )
  U3 <- phi(effSize)
  CL <- phi(effSize / sqrt(2))
  prior_effSize <- logdensity.norm(0, 0, 1)
  posterior_effSize <- logdensity.norm(0, effSize, 1)
  logBF <- prior_effSize - posterior_effSize
}
"
    y = model.frame(formula, data)[,1]
    group = as.numeric(as.factor(model.matrix(formula, data)[,2]))
    N = length(y)
    N1 = as.vector(table(group))[1]
    N2 = as.vector(table(group))[2]
    write_lines(two_sample_t, "two_sample_t.txt")
    jagsdata = list("N" = N, "N1" = N1, "N2" = N2, "group" = group, "y" = y)
    monitor = c("muDiff", "sigmaDiff", "mu", "sigma",  "effSize", "U3", "CL", "logBF", "ySim")
    inits = lapply(1:chains, function(z) list("mu" = c(mean(y), mean(y)), "ySim"  = y, "tau" = c(1/var(y),1/var(y)),
                                              .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:1000, 1)))
    out = run.jags(model = "two_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
    if (!is.null(cl)){
      parallel::stopCluster(cl = cl)
    }
    file.remove("two_sample_t.txt")
    return(out)
    
    }

if (model == "rm") {
  
  repeated_measures_t <- "
  
  model {
    mu ~ dnorm(0, 1e-100)
    tau ~ dscaled.gamma(1e100, 1)
    
    for (i in 1:N) {
      y[i] ~ ddexp(mu, tau)
      ySim[i] ~ ddexp(mu, tau)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - 0) / sigma
    CL <- phi(effSize / sqrt(2))
    prior_effSize <- logdensity.norm(0, 0, 1)
    posterior_effSize <- logdensity.norm(0, effSize, 1)
    logBF <- prior_effSize - posterior_effSize
}"
  
  if (is.null(data)) {
    cat(crayon::red("Please provide a vector of paired differences "))
  }
  
  y = data
  write_lines(repeated_measures_t, "repeated_measures_t.txt")
  jagsdata = list("N" = length(y), "y" = y)
  monitor = c("mu", "sigma",  "effSize", "CL", "logBF", "ySim")
  inits = lapply(1:chains, function(z) list("mu" = mean(y), "ySim"  = y, "tau" = 1/var(y),
                                            .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  out = run.jags(model = "repeated_measures_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("repeated_measures_t.txt")
  return(out)
  
}

if (model == "os") {
  
  one_sample_t <- "
  
  model {
    mu ~ dnorm(0, 1e-100)
    tau ~ dscaled.gamma(1e100, 1)
    
    for (i in 1:N) {
      y[i] ~ ddexp(mu, tau)
      ySim[i] ~ ddexp(mu, tau)
    }
    
    sigma <- sqrt(1/tau)
    effSize <- (mu - compval) / sigma
    prior_effSize <- logdensity.norm(compval, compval, 1)
    posterior_effSize <- logdensity.norm(compval, effSize, 1)
    logBF <- prior_effSize - posterior_effSize
    
}"
  
  if (is.null(data)) {
    cat(crayon::red("Please provide a vector of observations"))
  }
  if (is.null(compval)) {
    cat(crayon::red("Please provide a comparison value against which to test the mean of y"))
  }
  
  y = data
  write_lines(one_sample_t, "one_sample_t.txt")
  jagsdata = list("N" = length(y), "y" = y, compval = compval)
  monitor = c("mu", "sigma",  "effSize", "logBF", "ySim")
  
  inits = lapply(1:chains, function(z) list("mu" = mean(y), "ySim"  = y, "tau" = 1/var(y),
                                            .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
  
  out = run.jags(model = "one_sample_t.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, 
                 burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  
  file.remove("one_sample_t.txt")
  return(out)
  
}

  }


}

