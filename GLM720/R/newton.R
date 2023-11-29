newton <- function(func, init, eps=1e-16, maxiter=50){

  params <- init
  out <- matrix(NA, nrow = maxiter + 1, ncol = length(t(params)))
  out[1, ] <- t(params)
  i <- 1
  continue <- T

  while(continue){
    i <- i + 1
    funcOut <- func(params)
    params <- params - solve(funcOut$H)%*%funcOut$u
    if(sum(is.na(params)) > 0) stop("NA nas estimativas")
    out[i, ] <- t(params)
    continue <- (sqrt(t(funcOut$u)%*%funcOut$u) > eps) && (i <= maxiter)
  }

  if (i > maxiter) warning("Máximo número de iterações atingido")

  out <- out[!is.na(out[, 1]), ]
  output_list <- list(params = params, out = out, iter = i)
  return(output_list)

}
