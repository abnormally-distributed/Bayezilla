#' Estimate a Correlation Matrix 
#' 
#' 
#' @param x a data frame or matrix containing ONLY numeric variables to be correlated 
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
corMat = function(x, iter=10000, warmup=2500, adapt=2500, chains=4, thin=3, 
                   method = "parallel", cl = makeCluster(2), ...){
  
  if (!is.matrix(x)){
    y = as.matrix(x)
  } 
  
  jags_cor_matrix = "
model {
  for ( i in 1:N ) {
    y[i,1:P] ~ dmnorm(mu[1:P] , Tau[1:P,1:P]) 
  }
  
  for ( varIdx in 1:P ) {
    mu[varIdx] ~ dnorm( 0 , 1/2^2 ) 
  }
  
  # Estimate Precision Matrix
  Tau ~ dwish(priorScale[1:P, 1:P], df)
  
  # Convert Precision Matrix to sd and correlation:
  Sigma <- inverse( Tau )
  
  for ( varIdx in 1:P ) { 
    sqrtSigma[varIdx] <- sqrt(Sigma[varIdx,varIdx]) 
  }
    
  for ( varIdx1 in 1:P ) { 
    for ( varIdx2 in 1:P ) {
    Rho[varIdx1,varIdx2] <- ( Sigma[varIdx1,varIdx2] 
                               / (sqrtSigma[varIdx1]*sqrtSigma[varIdx2]) )
    }
  }
}
" 
write_lines(jags_cor_matrix, "jags_cor_matrix.txt")

jagsdata = list(
  "y"  =  as.matrix(y),
  "N"  =  nrow(y),
  "P"  =  ncol(y),
  "df" =  ncol(y)^2  /  4,
  # For wishart (dwish) prior on inverse covariance matrix:
  priorScale = diag(diag(pseudoinverse(as.matrix(cov(y)))))
)

monitor = c("Rho")
inits = lapply(1:chains, function(z) list(.RNG.name="lecuyer::RngStream", .RNG.seed= sample(1:10000, 1)))
out = run.jags(model = "jags_cor_matrix.txt", modules = "glm", monitor = monitor, data = jagsdata, n.chains = chains,
               burnin = warmup, sample = iter, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE, ...)
if (!is.null(cl)){
  parallel::stopCluster(cl = cl)
}
file.remove("jags_cor_matrix.txt")
return(out)
}
