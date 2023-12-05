#' @title Ajusta modelos com resposta BinomialNegativa(\eqn{\mu},\eqn{\phi})
#' @description
#' Esta função auxiliar ajusta modelos de regressão com resposta
#' Binomial Negativa. É possível (até o momento) utilizar
#'  apenas a função de ligação `log`.
#' @param y `vector` com a contagem de eventos (números
#' inteiros);
#' @param X `matrix` com as variáveis preditoras;
#' @param init `matrix` com \eqn{p} linhas e uma coluna com os
#' chutes iniciais dos valores dos \eqn{p} parâmetros a serem
#'  estimados;
#' @param phi_init `numeric` com o chute inicial do parâmetro
#' de dispersão \eqn{\phi} a ser estimado;
#' @param eps Valor de \eqn{\varepsilon}, a distância máxima
#'  entre duas iterações consecutivas do método de Newton para
#'  se dizer que houve convergência do método;
#' @param maxiter Número máximo de iterações a serem realizadas
#'  no método de Newton para obtenção das estimativas dos
#'  parâmetros;
#' @param gamma Variável do tipo 'numeric' entre 0 e 1, indicando a confiança desejada no intervalo a ser construído para os parâmetros.
#' @param link Função de ligação a ser utilizada no modelo. As
#'  funções de ligação implementadas (e a forma de utilizá-las)
#'  até o momento são:
#'  \itemize{
#'  \item{"logit": }{Ligação logística, dada por
#'   \eqn{g(\mu) = \log\left(\frac{\mu}{1-\mu} \right);}
#'  \item{"probit": }{Ligação probito, a inversa da f.d.a. da
#'   normal padrão, dada por \eqn{g(\mu) = \phi^{-1}(\mu)},
#'   onde \eqn{\phi(.)} denota o valor da f.d.a. de uma normal
#'   padrão avaliada em (.);}
#'   \item {"cloglog": } Ligação log-log complementar, dada por
#'   \eqn{g(1-\mu) = \log(-\log(1-\mu))};}
#'  }
#' @return Retorna uma lista com as informações do ajuste. A
#' descrição de cada item está apresentada abaixo:
#' \itemize{
#' \item out: `matrix` com as estimativas dos parâmetros geradas
#'  a cada iteração do método de Newton;
#' \item est: `matrix` com as estimativas (finais) dos parâmetros
#'   do modelo;
#' \item iter: Número de iterações necessárias até a convergência;
#' \item se: `vector` com os erros-padrão dos estimadores dos
#' parâmetros do modelo;
#' \item preditor: `matrix` com o valor de \eqn{\eta} (preditor
#'  linear) de cada observação;
#' \item v.ajust: `matrix` com os valores das médias ajustadas
#' para cada observação;
#' \item z.value: `matrix` com os valores das estatísticas-Z
#' do teste individual para cada \eqn{\beta_i = 0};
#' \item p.value: `matrix` com os p-valores dos testes individuais para cada \eqn{\beta_i = 0  \times  \beta_i \neq 0};
#' \item LI: `matrix` com os limites inferiores dos intervalos de confiança de \eqn{\gamma}% para os parâmetros;
#' \item LS: `matrix` com os limites superiores dos intervalos de confiança de \eqn{\gamma}% para os parâmetros;
#' \item AIC: AIC do modelo ajustado;
#' \item BIC: BIC do modelo ajustado;
#' \item desvio: Desvio do modelo ajustado;
#' \item gl.res: Número de graus de liberdade dos resíduos;
#' \item comp.desvio: Resíduos componentes do desvio do modelo
#'  ajustado;
#' \item residuos.pearson: Resíduos de Pearson do desvio do modelo
#'  ajustado;
#' \item comp.desvio.st: Resíduos componentes do desvio padronizados
#'  do modelo ajustado;
#' \item residuos.pearson.st: Resíduos de Pearson padronizados
#'  do modelo ajustado;
#'  \item Hhat: Matriz hessiana ajustada.
#' }
#'

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
