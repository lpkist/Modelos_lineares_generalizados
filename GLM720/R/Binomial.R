Binomial <- function(y, X, beta, link = "logit"){
  eta <- X%*%beta
  if(link == "logit"){
    pi <- exp(eta)/(1 + exp(eta))
    u <- t(X)%*%(y - pi)
    W <- diag(as.vector(pi*(1-pi)))
  } else if(link == "probit"){
    pi <- pnorm(eta)
    u <- t(X)%*%diag(as.vector(dnorm(eta)))%*%diag(as.vector(1/(pi*(1-pi))))%*%(y-pi)
    W <- diag(as.vector(dnorm(eta)^2/(pi*(1-pi))))
  } else if(link == "cloglog"){
    pi <- 1-exp(-exp(eta))
    u <- t(X)%*%diag(as.vector(exp(eta)/((1+exp(eta))^2)))%*%diag(as.vector(1/(pi*(1-pi))))%*%(y-pi)
    W <- diag(as.vector((exp(eta)/((1+exp(eta))^2))^2/(pi*(1-pi))))
  } else stop("Função de ligação incorreta. link deve estar em c('logit', 'probit', 'cloglog')")

  f <- dbinom(y, size = 1, prob = pi)


  H <- -t(X)%*%W%*%X
  output_list <- list(f = f, u = u, H = H)
  return(output_list)
}
