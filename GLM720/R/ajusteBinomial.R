#' @title Ajusta modelos com resposta Binomial(1,.)
#' @description
#' Esta função auxiliar ajusta modelos de regressão com resposta
#' Binomial a partir de dados desagrupados (ou seja, assume-se)
#' que \eqn{n_i =1} para todo \eqn{i=1,2,...,n}. É possível (até
#' o momento) utilizar as funções de ligação `logit`, `probit` e
#' `cloglog`, que são, respectivamente, as inversas das f.d.a.s
#' das distribuições logística, normal e valor extremo padrão.
#' @param y `vector` com o número de sucessos (0 ou 1);
#' @param X `matrix` com as variáveis preditoras;
#' @param init `matrix` com \eqn{p} linhas e uma coluna com os
#' chutes iniciais dos valores dos \eqn{p} parâmetros a serem
#'  estimados;
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


ajusteBinomial <- function(y, X, init, eps=1e-6, maxiter=50, gamma = 0.95, link = "logit"){

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
