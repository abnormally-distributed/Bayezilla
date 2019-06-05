#' Bivariate Correlation Test
#' 
#'    
#' 
#' @param x a data frame or matrix containing two variables to be correlated, OR a single vector to be paired with y.
#' @param y if x is a vector, put the second variable to be correlated here.  
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
#' corTest(iris$Sepal.Width, iris$Petal.Length)
#' 
corTest = function(x, y = NULL, iter=10000, warmup=2500, adapt=2500, chains=4, thin=3, 
                   method = "parallel", cl = makeCluster(2), ...){
  

if (is.data.frame(x) || is.matrix(x)){
  y = x
} else if (!is.null(y)) {
  y = cbind(x, y)
}

jags_cor_test = "
model {
  for ( i in 1:N ) {
    y[i,1:2] ~ dmnorm(mu[1:2] , Tau[1:2,1:2]) 
  }
  for ( varIdx in 1:2 ) {
    mu[varIdx] ~ dnorm( 0 , 1/2^2 ) 
  }
  
  # Estimate Precision Matrix
  Tau ~ dwish(priorScale[1:2, 1:2], 3)
  
  # Convert invCovMat to sd and correlation:
  Sigma <- inverse( Tau )
  
  for ( varIdx in 1:2 ) { 
    sqrtSigma[varIdx] <- sqrt(Sigma[varIdx,varIdx]) 
    }
  for ( varIdx1 in 1:2 ) { 
    for ( varIdx2 in 1:2 ) {
    Rho[varIdx1,varIdx2] <- ( Sigma[varIdx1,varIdx2] 
                               / (sqrtSigma[varIdx1]*sqrtSigma[varIdx2]) )
    }
  }
  rho <- Rho[2,1]
}
" 
write_lines(jags_cor_test, "jags_cor_test.txt")

jagsdata = list(
  "y" = as.matrix(scale(y)) ,
  "N"=  nrow(y) ,
  # For wishart (dwish) prior on inverse covariance matrix:
  priorScale = diag(x=1,nrow=ncol(y))  # Rmat = diag(apply(y,2,var))
)

monitor = c("rho")
inits = lapply(1:chains, function(z) list(.RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:100000, 1)))
out = run.jags(model = "jags_cor_test.txt", modules = "glm", monitor = monitor, data = jagsdata, n.chains = chains,
               burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
if (!is.null(cl)){
  parallel::stopCluster(cl = cl)
}
file.remove("jags_cor_test.txt")
return(out)
}
