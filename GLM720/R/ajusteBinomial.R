ajusteBinomial <-function(y, X, init, eps=1e-6, maxiter=50, gamma = 0.95, link = "logit"){

  n <- dim(X)[1]
  p <- dim(X)[2]

  outNewton <- newton(function(input) Binomial(y = y, X = X, beta = input, link = link),
                      init, eps, maxiter)
  est <- outNewton$params
  out <- outNewton$out
  iter <- outNewton$iter

  preditor <- X%*%est
  if(link == "logit"){
    v.ajust <- exp(preditor)/(1 + exp(preditor))
    u <- t(X)%*%(y - v.ajust)
    W <- diag(as.vector(v.ajust*(1-v.ajust)))
  } else if(link == "probit"){
    v.ajust <- pnorm(preditor)
    u <- t(X)%*%diag(as.vector(dnorm(preditor)))%*%diag(as.vector(1/(v.ajust*(1-v.ajust))))%*%(y-v.ajust)
    W <- diag(as.vector(dnorm(preditor)^2/(v.ajust*(1-v.ajust))))
  } else if(link == "cloglog"){
    v.ajust <- 1-exp(-exp(preditor))
    u <- t(X)%*%diag(as.vector(exp(-preditor)/((1+exp(-preditor))^2)))%*%diag(as.vector(1/v.ajust*(1-v.ajust)))%*%(y-v.ajust)
    W <- diag(as.vector((exp(-preditor)/((1+exp(-preditor))^2))^2/(v.ajust*(1-v.ajust))))
  }
  I <- t(X)%*%W%*%X
  se <- sqrt(diag(solve(I)))

  di <- -2*dbinom(y, size = 1, prob = v.ajust, log = T)
  desvio <- sum(di)
  AIC <- desvio + 2*p
  BIC <- desvio + 2*p*log(n)

  rP <- (y-v.ajust)/sqrt(v.ajust*(1-v.ajust))
  rD <- ifelse((y - v.ajust) > 0, sqrt(di), -sqrt(di))
  H <- diag(sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W))
  rPs <- rP/sqrt((1 - H))
  rDs <- rD/sqrt((1 - H))

  z.value <- est/se
  p.value <- 2*pnorm(abs(z.value), lower.tail = FALSE)

  z.alpha <- -qnorm((1-gamma)/2)
  LI <- est - z.alpha*se
  LS <- est + z.alpha*se

  saida <- list(out = out, est = est, iter = iter, se = se, preditor = preditor,
                v.ajust = v.ajust, z.value = z.value, p.value = p.value, LI = LI,
                LS = LS, AIC = AIC, BIC = BIC, desvio = desvio, gl.res = (n-p),
                comp.desvio = rD, residuos.pearson = rP,
                comp.desvio.st = rDs, residuos.pearson.st = rPs, Hhat = H)

  return(saida)

}
