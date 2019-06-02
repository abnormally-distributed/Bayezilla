#' Generate Random Correlation Matrices
#'
#' @param p The number of variabes to simulate
#' @param eta The regularization parameter. Larger values result in smaller correlations. Only applies to method = "lkj"
#' @param df The degrees of freedom. Only applies to method = "wishart"
#' @param One of "lkj" (the default) or "wishart"
#'
#' @return
#' a correlation matrix
#' @export
#'
#' @examples
#' rcorr(10, 2)
#' 
rcorr <-function(p, eta=1, df = NULL, method = "lkj") { 
  
  if (method == "lkj"){
  
  p<-as.integer(p)
  
  if(p<=0 || !is.integer(p))
    
  { 
    stop(cat(crayon::red("The dimension 'p' must be positive.\n")))
  }
  if(p==1){
    stop(cat(crayon::cyan("The dimension 'p' must be greater than 1.\n")))
  }
  if(eta<=0){
    stop(cat(crayon::blue("'eta' must be positive.\n")))
  }
  
  if(p==2){ 
    rho<-2*rbeta(1,eta,eta)-1
    rho<-matrix(c(1,rho,rho,1),2,2)
    return(rho) 
  }
  
  rho<-matrix(0,p,p)
  beta<-eta+(p-2)/2
  
  # step 1
  r12<-2*rbeta(1,beta,beta)-1
  rho<-matrix(c(1,r12,r12,1),2,2)
  
  # iterative steps
  for(m in 2:(p-1)){ 
    beta<-beta-0.5
    y<-rbeta(1,m/2,beta)
    z<-rnorm(m,0,1)
    znorm<-sqrt(sum(z^2))
    # random on surface of unit sphere
    z<-z/znorm
    w=sqrt(y)*z
    rhalf<-chol(rho)
    qq<-w%*%rhalf
    rho<-cbind(rho,t(qq))
    rho<-rbind(rho,c(qq,1))
  }
  
  # return rho
  return(rho)
}

  if (method == "wishart"){
      
      if (df < p){
        cat(crayon::magenta(crayon::bold("Please select a degrees of freedom equal to or greater than the number of dimensions.")))
      }
    
     cov2cor(rWishart(1, df = df, Sigma = diag(rep(1,p)))[,,1])
    
  }
  
}