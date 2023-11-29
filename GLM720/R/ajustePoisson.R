ajustePoisson <-function(y, X, init, eps=1e-6, maxiter=50, gamma = 0.95){

  n <- dim(X)[1]
  p <- dim(X)[2]

  outNewton <- newton(function(input) Poisson(y = y, X = X, beta = input),
                      init, eps, maxiter)
  est <- outNewton$params
  out <- outNewton$out
  iter <- outNewton$iter

  preditor <- X%*%est
  v.ajust <- exp(preditor)
  W <- diag(n)
  diag(W) <- v.ajust
  I <- t(X)%*%W%*%X
  se <- sqrt(diag(solve(I)))

  di <- 2*(y*(ifelse(y == 0, 0, log(y)) - log(v.ajust)) - y + v.ajust)
  desvio <- sum(di)
  AIC <- -2*(sum(dpois(y, v.ajust, log = TRUE))) + 2*p
  BIC <- -2*(sum(dpois(y, v.ajust, log = TRUE))) + 2*p*log(n)

  rP <- (y - v.ajust)/sqrt(v.ajust)
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
