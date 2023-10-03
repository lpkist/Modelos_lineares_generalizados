lm_Lucas <- function(formula, data, alpha = .05){
  # Informações básicas
  Y <- as.matrix(model.frame(formula, data = data)[,1])
  X <- model.matrix(formula, data = data)
  n <- nrow(X)
  p <- ncol(X)
  
  # Definição das matrizes de interesse
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  J <- matrix(1, nrow = n, ncol = n)
  A <- H - 1/n*J
  B <- diag(n) - H
  
  # Estimação
  betas <- solve(t(X)%*%X)%*%t(X)%*%Y
  sigma2 <- 1/(n-p)*t((Y-H%*%Y))%*%(Y-H%*%Y)
  sigma <- sqrt(sigma2)
  cov_betas <- sigma2[[1]]*solve(t(X)%*%X)
  var_betas <- diag(cov_betas)  
  
  # Testes t
  estat_t <- betas/sqrt(var_betas)
  p_valor <- pt(abs(estat_t), n-p, lower.tail = F)
  
  # ICS
  li <- qt(alpha/2, n-p)
  IC_inf <- betas + li*var_betas
  IC_sup <- betas - li*var_betas
  
  # Resíduos
  res <- Y - H%*%Y

  # ANOVA
  SQM <- t(Y)%*%A%*%Y
  SQR <- t(Y)%*%B%*%Y
  SQT <- SQM+SQR
  
  # R2
  R2 <- round(SQM/SQT,4)
  R2_adj <- round(1 - SQR/SQT*(n-1)/(n-p),4)
  
  cat("Resíduos:\n")
  print(summary(c(res)))
  cat("\n")
  cat("Coeficientes")
  
  estimativas <- cbind(round(betas, 5), round(sqrt(var_betas), 5),
                       round(estat_t, 3), p_valor,
                       paste0("[",round(IC_inf, 3), ";", round(IC_sup, 3), "]"),
        ifelse(p_valor < 1e-3, "***",
               ifelse(p_valor < 1e-2, "**",
                      ifelse(p_valor < .05, "*",
                             ifelse(p_valor < .1, ".","")))))

  colnames(estimativas) <- c("Estimativa", "SD", "Estat. t", "p-valor",
                             paste0("ICs de ", 100*(1-alpha),"%"), "")
  print(estimativas)
  cat("---")
  print("Códigos de significância: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  cat("\n")
  cat(c("Desvio-padrão dos resíduos: ", round(sigma, 4), "com", n-(p+1), "graus de liberdade\n"))
  cat(c("R² múltiplo:", R2, "R² ajustado:", R2_adj, "\n"))
  cat("Estatística F:", round((SQM/(p-1))/(SQR/(n-p)), 2), "com",
      p-1, "e", n-p, "graus de liberdade, p-valor",
      pf((SQM/(p-1))/(SQR/(n-p)), df1 = p-1, df2 = n-p, lower.tail = F))  

  output <- list(X = X, Y = Y, betas = betas, sigma2 = sigma2, 
                 coefs = c(betas, sigma2), res = res, R2 = R2,
                 R2_adj = R2_adj, n = n, p = p, cov = cov_betas,
                 H = H, A = A, B = B, estat_t = estat_t,
                 p_valor = p_valor, LIC = IC_inf, LSC = IC_sup,
                 ajustados = H%*%Y, SQM = SQM, SQR = SQR, SQT = SQT)
  return(output)
}
library(wooldridge)
lm_Lucas(educ ~ motheduc + fatheduc + abil + I(abil^2), htv, .05)
summary(lm(educ ~ motheduc + fatheduc + abil + I(abil^2), htv))
modelo <- lm_Lucas(educ ~ motheduc + fatheduc + abil + I(abil^2), htv, .05)
C <- matrix(c(matrix(0,nrow=4), diag(4)),nrow = 4)


teste_cbeta <- function(modelo, C){
  betas <- modelo$betas
  sigma2 <- modelo$sigma2
  X <- modelo$X
  q <- nrow(C)
  p <- ncol(C)
  n <- nrow(X)
  
  cbeta <- C%*%betas
  Q_estrela <- 1/sigma2*t(cbeta)%*%solve(C%*%solve(t(X)%*%X)%*%t(C))%*%cbeta
  p_valor <- pf(Q_estrela, q, n-p, lower.tail = F)
  
  cat(c("Teste C beta\n"))
  cat("Matriz C:\n")
  print(C)
  cat("Betas:\n")
  print(betas)
  cat("\nEstatística do teste:", Q_estrela)
  cat("\np-valor do teste:", p_valor)
  output <- list(q = q, n = n, p = p, betas = betas, sigma2 = sigma2,
                 Q_star = Q_estrela, p_valor = p_valor)
  return(output)
}
teste_cbeta(modelo, C)
