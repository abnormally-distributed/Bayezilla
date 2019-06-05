#' One or Two Sample Binomial Tests
#'
#' @param z a numeric value giving the number of succesful trials. If two sample, a vector such as c(1, 3)
#' @param N a numeric value giving the total number of trials. If two sample, a vector such as c(5, 5)
#' @param shape the prior value for the shape parameters of the beta prior (defaults to c(0.50, 0.50))
#' @param one.sample Defaults to TRUE
#' @param iter the number of iterations. defaults to 10000.
#' @param warmup number of burnin samples. defaults to 2500.
#' @param adapt number of adaptation steps. defaults to 2500.
#' @param chains number of chains. defaults to 4.
#' @param thin the thinning interval. defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... other arguments to run.jags
#'
#' @return a runjags object
#' @export
#' @examples binomTest()
#' 
binomTest<- function(z, N, one.sample = TRUE, shape = c(0.50, 0.50), iter=10000,
                          warmup=2500, adapt=2500, chains=4, thin=3, 
                          method = "parallel", cl = makeCluster(2), summarise = FALSE, ...){
  if (isTRUE(one.sample)) {
    
    jags_binomial <- "model {
                        phi ~ dbeta(a, b)
                        ySim ~ dbinom(phi, N)
                        y ~ dbinom(phi, N)
                    }"
    
    jagsdata = list("N" = N, "y" = z, "a" = shape[1], "b" = shape[2])
    write_lines(jags_binomial, "jags_binomial.txt")
    monitor = c("phi", "ySim")
    inits = lapply(1:chains, function(z) list("phi" = z/N, .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
    
  }
  else if (isFALSE(one.sample)) {
    
    jags_binomial <- "model {
                            for(i in 1:2) {
                                y[i] ~ dbinom(phi[i], N[i])
                                phi[i] ~ dbeta(a, b)
                                ySim[i] ~ dbinom(phi[i], N[i])
                              }
                                phiDiff <- phi[1] - phi[2]
                          }"
             jagsdata = list("N" = N, "y" = z, "a" = shape[1], "b" = shape[2])
             write_lines(jags_binomial, "jags_binomial.txt")
             monitor = c("phi", "phiDiff", "ySim")
             inits = lapply(1:chains, function(z) list("phi" = c(.5,.5), .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
             
  }
  
  out = run.jags(model = "jags_binomial.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_binomial.txt")
  return(out)

}


#' Poisson Test
#'
#' @param x number of events. A vector of length one or two.
#' @param t	 time base for event count. A vector of length one or two.
#' @param r	 comparison value (only applies to one sample) on scale of lambda. Ie, on the scale of lambda_raw/t. Defaults to 1 and should typically stay 1 unless you have a specific hypothesis.
#' @param shra the shape and rate parameters of the gamma prior. Defaults to c(1e-6, 1e-6) which is non-informative and yields
#' a posterior mean typically equal to the maximum likelihood estimate.
#' @param one.sample Defaults to TRUE
#' @param iter the number of iterations. defaults to 10000.
#' @param warmup number of burnin samples. defaults to 2500.
#' @param adapt number of adaptation steps. defaults to 2500.
#' @param chains number of chains. defaults to 4.
#' @param thin the thinning interval. defaults to 3.
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel" or. Otherwise, "rjags" (single core run).
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores.
#' @param ... other arguments to run.jags
#' 
#' @return an rjags object
#' 
#' @export
#' @examples poisTest()
#' 
poisTest <- function(x, t = c(1, 1), r = 1, shra = c(1e-6, 1e-6), one.sample = TRUE, iter=10000, warmup=2500, adapt=2500, chains=4, thin=3, method = "parallel", cl = makeCluster(2), summarise = FALSE, ...){
  
  if (isTRUE(one.sample)) {
    
    jags_poisson <- "model {
        lambda ~ dgamma(sh, ra)
        y ~ dpois(lambda * t)
        ySim ~ dpois(lambda * t)
        rateDiff <-  lambda - compval
        rateRatio <- lambda / compval
      }"
    
        jagsdata = list("y" = x[1], "t" = t[1], "compval" = r, "sh" = shra[1], "ra" = shra[2])
        write_lines(jags_poisson, "jags_poisson.txt")
        monitor = c("lambda", "rateDiff", "rateRatio", "ySim")
        inits = lapply(1:chains, function(z) list("lambda" = 1, .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
    
    
  }
  else if (isFALSE(one.sample)) {
    
    two_sample_poisson<- "model {
        for(group in 1:2) {
            y[group] ~ dpois(lambda[group] * t[group])
            lambda[group] ~ dgamma(sh, ra)
            ySim[group] ~ dpois(lambda[group] * t[group])
        }
            rateDiff <-  lambda[1] - lambda[2]
            rateRatio <- lambda[1] / lambda[2]
        }"
    
        jagsdata = list("y" = x, "t" = t, "sh" = shra[1], "ra" = shra[2])
        write_lines(jags_poisson, "jags_poisson.txt")
        monitor = c("lambda", "rateDiff", "rateRatio", "ySim")
        inits = lapply(1:chains, function(z) list("lambda" = c(1,1), .RNG.name = "lecuyer::RngStream", .RNG.seed = sample(1:10000, 1)))
  }
  out = run.jags(model = "jags_poisson.txt", modules = "glm", monitor = monitor, data = jagsdata, inits = inits, burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
  if (!is.null(cl)){
    parallel::stopCluster(cl = cl)
  }
  file.remove("jags_poisson.txt")
  return(out)
}