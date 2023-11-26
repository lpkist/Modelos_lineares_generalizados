data {
  int<lower=1> n;     // tamanho da amostra
  int<lower=0> y[n];  // observação i
  int<lower=0> m[n];  // m da binomial
  vector[n] x;   // matriz de covariáveis
}

parameters {
  real beta0;  //parâmetros
  real beta1;
}

transformed parameters{
  vector[n] pest;
  for (i in 1:n)
  {
    pest[i] = exp(beta0 + beta1*x[i])/(1+exp(beta0+beta1*x[i])); //  probabilidade estimada
  }
}

model{
  beta0 ~ normal(0,10);// média, desvio-padrão
  beta1 ~ normal(0,10);// média, desvio-padrão
  for (i in 1:n)
  {
    y[i] ~ binomial(m[i],inv_logit(beta0 + x[i]*beta1));
  }
}

generated quantities {
  real<lower=0> yrep[n];
  real log_lik[n];
  real prep[n];
  for(i in 1:n)
  {
    yrep[i] = binomial_rng(m[i],pest[i]);
    prep[i] = yrep[i]/m[i];
    log_lik[i] = binomial_lpmf(y[i]|m[i],pest[i]);
  }
}

