BinomialNegativa <- function(y, X, beta, k){
  n <- length(y)
  mu <- exp(X%*%beta)
  f <- gamma(y+k)/(gamma(k)*gamma(y+1))*(mu/(mu+k))^y*(k/(mu+k))^k
  f_g <- diag(n)
  diag(f_g) <- mu
  W <- diag(n)
  diag(W) <- mu/(1+mu*k)
  u <- t(X)%*%W%*%solve(f_g)%*%(y-mu)
  H <- -t(X)%*%W%*%X
  U_k <- sum(digamma(k+y)-digamma(k)-(y+k)/(k+mu)+log(k/(k+mu))+1)
  L_kk <- sum(trigamma(k+y)+(y-2*mu-k)/(k+mu)^2) +n/k*(1-k*trigamma(k))
  output_list <- list(f = f, mu = mu, u = u, H = H, U_k = U_k, L_kk = L_kk)
  return(output_list)
}
