N <- 60
P <- 12
betas = c(.89, 0, 0, 0, 1.49, 0, -3.51, 0, 3.08, 
          -2.18113776, -3.035, 0)
mu <- 0
set.seed(418)
X = MASS::mvrnorm(N, mu = rep(0, P), Sigma = rcorr(P), empirical = TRUE)
X = as.matrix(scale(X))
y = scale(mu + as.vector(X %*% matrix(betas, ncol = 1)))
betas = Bayezilla::lmSolve(y ~ ., cbind.data.frame(y = y, X))[-1]
noise = runif(N*P, -.132, .132)
y = scale(mu + as.vector(X %*% matrix(betas, ncol = 1)) + rnorm(N, 0, .5))
X = scale(X + matrix(noise, N, P))
demo.dat3 = cbind.data.frame(y, X)
colnames(demo.dat3) = c("y", paste0("X", 1:P))
betas = round(c(0, betas), 2)
names(betas) = c("Intercept", paste0("V", 1:P))
attr(demo.dat3, "true.betas") = betas

N <- 60
P <- 12
betas = c(.89, -.56, .2644, .682351, 1.493321231, 
          -.114213215, -3.51, -.0921231231, 3.08, 
          -2.18113776, -2.935, -.35361351465465)
mu <- 0
set.seed(418)
X = MASS::mvrnorm(N, mu = rep(0, P), Sigma = rcorr(P), empirical = TRUE)
X = as.matrix(scale(X))
y = scale(mu + as.vector(X %*% matrix(betas, ncol = 1)))
betas = Bayezilla::lmSolve(y ~ ., cbind.data.frame(y = y, X))[-1]
noise = runif(N*P, -.132, .132)
y = scale(mu + as.vector(X %*% matrix(betas, ncol = 1)) + rnorm(N, 0, .5))
X = scale(X + matrix(noise, N, P))
demo.dat4 = cbind.data.frame(y, X)
demo.dat4$X2 <- NULL
colnames(demo.dat4) = c("y", paste0("X", 1:P-1))
betas = round(c(0, betas[-2]), 2)
names(betas) = c("Intercept", paste0("V", 1:P-1))
attr(demo.dat4, "true.betas") = betas
