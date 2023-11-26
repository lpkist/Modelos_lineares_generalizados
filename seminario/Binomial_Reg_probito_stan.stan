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
    pest[i] = Phi(beta0 + beta1*x[i]); 
  }
}

model{
  beta0 ~ normal(0,10); // média, desvio-padrão
  beta1 ~ normal(0,10); // média, desvio-padrão
  for (i in 1:n)
  {
    y[i] ~ binomial(m[i],pest[i]);
  }
}

generated quantities {
  int<lower=0> yrep[n];
  real log_lik[n];
  for(i in 1:n)
  {
    yrep[i] = binomial_rng(m[i],pest[i]);
    log_lik[i] = binomial_lpmf(y[i]|m[i],pest[i]);
  }
}
