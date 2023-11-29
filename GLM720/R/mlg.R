#' @title Ajusta modelos lineares generalizados (MLG)
#' @description
#' A função realiza o ajuste de vários modelos lineares generalizados
#' manualmente, retornando um intervalo de confiança para os
#' parâmetros estimados, além dos objetos retornados pela função
#' `glm`.
#' @param f Fórmula do ajuste a ser realizado;
#' @param dados Conjunto de dados do tipo `data.frame` com as colunas especificadas na fórmula (parâmetro `f`);
#' @param init `matrix` com uma coluna e `p` linhas, onde `p` é o
#'  número de parâmetros incluídos no modelo. Esses valores serão
#'   utilizados como chutes iniciais para o ajuste via método de
#'   Newton-Raphson;
#' @param phi_init objeto do tipo `numeric` com o valor inicial
#' para o parâmetro \eqn{\phi} quando se está estimando um modelo
#' binomial negativo;
#' @param link Função de ligação a ser utilizada no modelo. As funções de ligação implementadas até o momento são:
#' \itemize{
#'  \item{Binomial}{"logit" (ligação logpistica), "probit" (ligação
#'   probito, a inversa da f.d.a. da normal padrão) e "cloglog" (
#'   ligação log-log complementar);}
#'  \item{Binomial negativa}{"log"}
#'  \item{Poisson}{"log"}
#' }
#' @param distr Distribuição assumida para a variável resposta. Até
#' o momento, estão implementadas as distribuições "Binomial",
#' "Binomial negativa" e "Poisson";
#' @param eps Valor de \eqn{\varepsilon}, a distância máxima entre
#' duas iterações consecutivas do método de Newton para se dizer
#' que houve convergência do método;
#' @param maxiter Número máximo de iterações a serem realizadas no
#' método de Newton para obtenção das estimativas dos parâmetros;
#' @param show.result Variável do tipo `logical` indicando se é
#' de interesse imprimir os valores dos parâmetros e outras
#' informações sumárias;
#' @param gamma Variável do tipo `numeric` entre 0 e 1, indicando
#' a confiança desejada no intervalo a ser construído para os
#' parâmetros.
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
#' \item mv.ajust: `matrix` com os valores das médias ajustadas
#' para cada observação;
#' \item z.values: `matrix` com os valores das estatísticas-Z
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
#' @importFrom MASS rnegbin
#' @export
mlg <- function(f, dados = NULL, init, phi_init, link = "logit", distr, eps=1e-6, maxiter=100, show.result = TRUE, gamma = 0.95){

  y <- as.matrix(model.frame(f, data = dados)[,1])
  X <- model.matrix(f, data = dados)
  n <- dim(X)[1]
  p <- dim(X)[2]

  # Validações
  if(p != nrow(init)) stop("Valores inválidos para os parâmetros iniciais")
  if(maxiter <= 0 | maxiter%%1 != 0) stop("maxiter precisa ser um inteiro positivo")
  if(eps <=0 | eps > 1) stop("O erro precisa pertencer ao  intervalo (0,1]")

  # Outputs
  if(distr == "Poisson") out = ajustePoisson(y, X, init, eps, maxiter, gamma = gamma)
  if(distr == "Binomial") out = ajusteBinomial(y, X, init, eps, maxiter, gamma = gamma, link = link)
  if(distr == "Binomial Negativa") out = ajusteBinomialNegativa(y = y, X = X, init = init, phi_init = phi_init, eps = eps, maxiter = maxiter, gamma = gamma, f=f, dados = dados)

  # Criterios
  desvio <- round(out$desvio,3)
  desvioinf <- matrix(c(desvio, out$gl.res), 1, 2)
  AIC    <- round(out$AIC,3)
  BIC    <- round(out$BIC,3)
  iter   <- out$iter

  Est <- cbind(round(out$est, 5), round(out$se, 5), round(out$z.value, 2), out$p.value)
  colx <- ncol(as.matrix(X))
  namesx <- colnames(X)
  LI <- round(out$LI, 3); LS <- round(out$LS, 3)
  IC <- paste0('[', LI[1], '; ', LS[1], ']')
  if(colx > 1){
    for(i in 2:colx) IC <- rbind(IC, paste0('[', LI[i], '; ', LS[i], ']')) }

  Estimativas <- data.frame(Est, IC = IC)
  colnames(Estimativas) <- c('Estimativas', 'SE', 'valor z', 'p-valor', 'ICs')
  Estimativas["Signif."] <-
    ifelse(Estimativas["p-valor"] < 1e-3, "***",
           ifelse(Estimativas["p-valor"] < 1e-2, "**",
                  ifelse(Estimativas["p-valor"] < .05, "*",
                         ifelse(Estimativas["p-valor"] < .1, ".",""))))
  criteria <- as.matrix(c(AIC, BIC))
  dimnames(criteria) <- list(c('AIC:', 'BIC:'), c(' '))
  dimnames(desvioinf) = list(c('Desvio:  '), c(' ', 'g.l'))
  dtn <- deparse(substitute(dados))

  ligFunc <- ifelse(distr == 'Poisson', 'log', 'logito')

  if(show.result == TRUE){
    cat('MLG para o modelo ', paste0(distr, ' com fução de ligação ', ligFunc),  '\n')
    cat('Modelo: ', paste0(deparse(f), ', dados = ', dtn), '\n')
    cat('---\n')
    print(Estimativas)
    cat("Códigos de significância:\n 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
    cat("\n")
    cat('---\n')
    print(desvioinf)
    print(criteria)
    cat('---\n')
    cat('Número de iterações: ', iter)
  }
  if(distr == "Binomial Negativa"){
    cat(c("\nEstimativa de k: ", round(out$k, 3)))
    cat(c("\nEstimativa de phi: ", round(1/out$k, 3)))
  }
  class(out) <- "mlgPoisson"
  return(invisible(out))

}

