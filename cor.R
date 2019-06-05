y = cbind(data$bmi, data$hdl)

jagsdata = list(
  "y" = as.matrix(scale(y)) ,
  "N"=  nrow(y) ,
  # For wishart (dwish) prior on inverse covariance matrix:
  priorScale = diag(x=1,nrow=ncol(y))  # Rmat = diag(apply(y,2,var))
)

jags_cor_test = "
model {
  for ( i in 1:N ) {
    y[i,1:2] ~ dmnorm(mu[1:2] , Tau[1:2,1:2]) 
  }
  for ( varIdx in 1:2 ) {
    mu[varIdx] ~ dnorm( 0 , 1/2^2 ) 
  }
  
  # Estimate Precision Matrix
  Tau ~ dwish(priorScale[1:2,1:2], 3)
  
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
" # close quote for modelString

monitor = c("rho")
write_lines(jags_cor_test, "jags_cor_test.txt")
out = run.jags("jags_cor_test.txt", monitor = monitor, data = jagsdata, n.chains = 4)
post_summary(out)
