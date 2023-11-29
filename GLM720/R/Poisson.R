Poisson <- function(y, X, beta){
  mu <- exp(X%*%beta)
  f <- dpois(y, mu)
  u <- t(X)%*%(y - mu)
  W <- diag(as.vector(mu))
  H <- -t(X)%*%W%*%X
  output_list <- list(f = f, u = u, H = H)
  return(output_list)
}

