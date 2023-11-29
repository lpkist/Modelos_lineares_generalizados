ajusteBinomialNegativa <-function(y, X, init, phi_init, eps=1e-6, maxiter=50, gamma = 0.95, f, dados){

  n <- dim(X)[1]
  p <- dim(X)[2]


  outNewton <- newtonBN(function(input, phi_init) BinomialNegativa(y = y, X = X, beta = input, k = phi_init),
                        init, phi_init, eps, maxiter)


  est <- outNewton$params
  k <- outNewton$k
  out <- outNewton$out
  iter <- outNewton$iter
  out_k <- outNewton$out_k

  preditor <- X%*%est
  v.ajust <- exp(preditor)
  W <- diag(n)
  diag(W) <- v.ajust/(1+v.ajust/k)
  I <- t(X)%*%W%*%X
  se <- sqrt(diag(solve(I)))

  di <- 2*(ifelse(y==0,0,y*log(y/v.ajust)-(y+k)*log((1+y/k)/(1+v.ajust/k))))
  desvio <- sum(di)
  AIC <- -2*sum(dnbinom(y, mu = v.ajust, size = k, log = T))+2*p
  BIC <- -2*sum(dnbinom(y, mu = v.ajust, size = k, log = T))+2*p*log(n)

  rP <- (y-v.ajust)/sqrt(v.ajust+v.ajust^2/k)
  rD <- ifelse((y - v.ajust) > 0, sqrt(di), -sqrt(di))
  H <- diag(sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W))
  rPs <- rP/sqrt((1 - H))
  rDs <- rD/sqrt((1 - H))

  z.value <- est/se
  p.value <- 2*pnorm(abs(z.value), lower.tail = FALSE)

  z.alpha <- -qnorm((1-gamma)/2)
  LI <- est - z.alpha*se
  LS <- est + z.alpha*se

  saida <- list(out = out, est = est, k=k, iter = iter, se = se, preditor = preditor,
                v.ajust = v.ajust, z.value = z.value, p.value = p.value, LI = LI,
                LS = LS, AIC = AIC, BIC = BIC, desvio = desvio, gl.res = (n-p),
                comp.desvio = rD, residuos.pearson = rP,
                comp.desvio.st = rDs, residuos.pearson.st = rPs, Hhat = H)

  return(saida)

}
