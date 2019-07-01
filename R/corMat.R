#' Correlation Matrix Estimation (also returns partial correlations)
#' 
#' Just a simple multivariate norma-wishart conjugate model that returns the standardized inverse precision matrix (correlation matrix)
#' and standardized precision matrix (partial correlations). 
#' 
#' @param x a data frame or matrix 
#' @param df degrees of freedom for wishart prior. Defaults to ncol(X) + 1
#' @param iter the number of iterations. defaults to 4000.
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
#' corMat(iris$Sepal.Width, iris$Petal.Length)
#' 
corMat = function(x, df = "default", iter = 4000, warmup=2500, adapt=2500, chains=4, thin=1, 
                   method = "parallel", cl = makeCluster(2), ...){

x = scale(x)
x = as.matrix(x)
  
jags_cor_test = "

model {
  for ( i in 1:N ) {
    y[i,1:P] ~ dmnorm(mu[1:P] , Tau[1:P,1:P]) 
  }
  for ( varIdx in 1:P ) {
    mu[varIdx] ~ dnorm( 0 , 1e-200) 
  }
  
  # Estimate precision matrix
  Tau ~ dwish(priorScale[1:P, 1:P], df)
  
  # Convert precision matrix to correlation:
  Sigma <- inverse(Tau)
  
  for ( varIdx in 1:P ) { 
    sqrtSigma[varIdx] <- sqrt(Sigma[varIdx,varIdx]) 
    }
  for ( varIdx1 in 1:P ) { 
    for ( varIdx2 in 1:P ) {
    Rho[varIdx1,varIdx2] <- ( Sigma[varIdx1,varIdx2] 
                               / (sqrtSigma[varIdx1]*sqrtSigma[varIdx2]) )
    }
  }
  
  # Convert precision matrix to partial correlation:
  for ( varIdx in 1:P ) { 
    sqrtinvTau[varIdx] <- sqrt(Tau[varIdx,varIdx]) 
  }
  
  for ( varIdx1 in 1:P ) { 
    for ( varIdx2 in 1:P ) {
    parRho[varIdx1,varIdx2] <-  -1 *  ( Tau[varIdx1,varIdx2] 
                               / (sqrtinvTau[varIdx1]*sqrtinvTau[varIdx2]))
    }
  }
}" 
write_lines(jags_cor_test, "jags_cor_test.txt")
if (df == "default"){
  df = ncol(x) + 1
} else {
  df = df
}

jagsdata = list(
  "y" = x,
  "N"=  nrow(x),
  "P" = ncol(x),
  # For wishart (dwish) prior on inverse covariance matrix:
  "priorScale" = diag(rep(1, ncol(x))),  # Rmat = diag(apply(y,2,var))
  "df" = df
)
monitor = c("Rho", "parRho")
inits = lapply(1:chains, function(z) list("Tau" = pseudoinverse(cov(x)), "mu" = colMeans(x), .RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
out = run.jags(model = "jags_cor_test.txt", modules = "glm", monitor = monitor, data = jagsdata, n.chains = chains,
               burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
if (!is.null(cl)){
  parallel::stopCluster(cl = cl)
}
file.remove("jags_cor_test.txt")
return(out)
}