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

