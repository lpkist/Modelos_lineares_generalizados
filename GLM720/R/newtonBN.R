newtonBN <- function(func, init, k, eps=1e-16, maxiter=50){

  params <- init
  k <- k
  out <- matrix(NA, nrow = maxiter + 1, ncol = length(t(params)))
  out[1, ] <- t(params)
  out_k <- matrix(NA, nrow = maxiter + 1, ncol = 1)
  out_k[1, ] <- k
  i <- 1
  continue <- T

  while(continue){
    i <- i + 1
    funcOut <- func(params, k)
    params <- params - solve(funcOut$H)%*%funcOut$u
    if(sum(is.na(params)) > 0) stop("NA nas estimativas")
    out[i, ] <- t(params)
    k <- k - funcOut$U_k/funcOut$L_kk
    if(is.na(k)) stop("NA nas estimativas")
    out_k[i, ] <- k
    continue <- (sqrt(t(funcOut$u)%*%funcOut$u) > eps) && (sqrt(t(funcOut$U_k)%*%funcOut$U_k) > eps)      &&(i <= maxiter)
  }

  if (i > maxiter) warning("Máximo número de iterações atingido")

  out <- out[!is.na(out[, 1]), ]
  out_k <- out_k[!is.na(out_k[, 1]), ]
  output_list <- list(params = params, k = k,
                      out = out, iter = i, out_k = out_k)
  return(output_list)

}
