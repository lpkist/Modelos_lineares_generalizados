library(coda)
library(bayesplot)
library(loo)
library(rstan)
library(TeachingDemos)
library(moments)
library(tidyverse)
library(patchwork)

set.seed(1234223)
# Dados
# coloca tudo no mesmo diretório e seta lá
dados <- read_table("meninas.txt")
x <- dados$idade_m # covariável
m <- dados$n_entrevistadas # m da binomial(m,p) para cada observação
y <- dados$n_apresentou # número de sucessos
xp <- x-mean(x) # tem que centralizar a covariável pra convergir
data <- list(y=y,x=xp,m=m,n=length(y)) # tem que deixar numa lista pra rodar

# caminho para o arquivo com o modelo 
folder <- getwd()

# Parâmetros MCMC

pars1 <- c("beta0","beta1","pest") # cadeias de interesse
pars3 <- c("beta0","beta1") # para rodar os gráficos
pars2 <- c("yrep","log_lik") # outras coisas de interesse
pars <- c(pars1,pars2)


############## logito
# Parâmetros MCMC
nchains <- 1 # são 3 pra ver se converge
niter <- 20500 # número de iterações total
nburnin <- 500 # número de observações iniciais jogadas fora
nthin <- 10 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
# caminho do .stan
model<- paste0(folder, "\\Binomial_Reg_stan.stan")

simullBSgs_inf<-stan(file = model, data = data, pars=pars,
                     chains = nchains, iter = niter,
                     warmup = nburnin, thin = nthin)
simullBSgs_inf
summary(simullBSgs_inf)

# Extraindo a cadeia (caso se queira usar outros
# pacotes e para fazer alguas contas)
mchain <- (rstan::extract(simullBSgs_inf,pars=pars,
                          permuted=FALSE,inc_warmup=FALSE))
result1 <- data.frame(mcmc(mchain[,1,]))
str(result1)


meu_tema <- theme_bw() + theme(text = element_text(size = 17))
# Estimativas básicas
mcmc_intervals(simullBSgs_inf,pars=pars3) + meu_tema # IC p/ os betas


#pdf(paste(file.save,sep="","HistInf",nome_file,".pdf"))
# plot_hist <- mcmc_hist(simullBSgs_inf,pars=pars3, facet_args = list(labeller = label_bquote(cols = list(beta_0, beta_1 ))))+meu_tema # hist das posterioris
plot_hist <- data.frame("beta0" =
                          simullBSgs_inf@sim$samples[[1]]$beta0[51:2050],
           "beta1" = simullBSgs_inf@sim$samples[[1]]$beta1[51:2050]) %>%
  pivot_longer(cols = 1:2, names_to = "beta", values_to = "aa") %>% 
  mutate(expr = ifelse(beta == "beta0", "beta[0]", "beta[1]")) %>% 
  ggplot(aes(x = aa))+
  geom_histogram(color = "blue", fill = "lightblue")+
  facet_wrap(~expr, scales = "free",
             labeller = label_parsed)+meu_tema+
  theme(axis.title = element_blank())
#dev.off()
#pdf(paste(file.save,sep="","DensInf",nome_file,".pdf"))
#plot_dens <- mcmc_dens(simullBSgs_inf,pars=pars3)+meu_tema # dens das posterioris
#dev.off()

plot_dens <- data.frame("beta0" =
                          simullBSgs_inf@sim$samples[[1]]$beta0[51:2050],
  "beta1" = simullBSgs_inf@sim$samples[[1]]$beta1[51:2050]) %>%
  pivot_longer(cols = 1:2, names_to = "beta", values_to = "aa") %>% 
  mutate(expr = ifelse(beta == "beta0", "beta[0]", "beta[1]")) %>% 
  ggplot(aes(x = aa))+
  geom_density(color = "blue", fill = "lightblue")+
  facet_wrap(~expr, scales = "free",
             labeller = label_parsed)+meu_tema+
  theme(axis.title = element_blank())


pdf("post_logito.pdf", width = 7, height = 4)
plot_hist / plot_dens
dev.off()

res<-summary(simullBSgs_inf)$summary

# Critérios de Informação e ajuste do modelo
# Checagem preditiva a posteriori
yrep<-(rstan::extract(simullBSgs_inf,
                      pars=c("yrep"))$yrep[order(sample(500)),])
y<-data$y

densmedia <- ppc_stat(y,yrep,stat="mean")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Média")

densvar <- ppc_stat(y,yrep,stat="var")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Variância")

densca <- ppc_stat(y,yrep,stat=function(y) skewness(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +labs(title = "Coef. de assimetria")

denscurt <-ppc_stat(y,yrep,stat=function(y) kurtosis(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+ labs(title = "Curtose")

pdf("momentos_logito.pdf", width = 10, height = 4)
densmedia + densvar + densca + denscurt
dev.off()

########Usando a proporção
m<-data$m
p<- y/m
prep<- yrep/t(matrix(m,length(y),500))

densprop<-ppc_dens_overlay(p,prep)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

cdfprop<-ppc_ecdf_overlay(p,prep)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

intervalprop<- ppc_intervals(p, prep)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ppc_ribbon(p, prep)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

scatterprop <- ppc_scatter_avg(p,prep)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("prop_logito.pdf", width = 7, height = 4)
densprop + cdfprop + intervalprop + scatterprop
dev.off()

# LOO e WAIC (comparação de modelos)
m_log_lik <-extract_log_lik(simullBSgs_inf,
                            parameter_name="log_lik",
                            merge_chains = FALSE)
loo(m_log_lik)
waic(m_log_lik)

library(xtable)
# Estimativas
result_est<-cbind(res[1:length(pars3),c(1,6,3,4,8)],
                  matrix(t(apply(mchain[,1,1:length(pars3)],2,emp.hpd,0.95)),length(pars3),2))
colnames(result_est)[c(6,7)] <-c("HPD(2.5%)","HPD(97.5%)")
xtable(result_est)
#











############## cloglog
# Parâmetros MCMC
nchains <- 1 # são 3 pra ver se converge
niter <- 210000 # número de iterações total
nburnin <- 10000 # número de observações iniciais jogadas fora
nthin <- 100 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
# caminho do .stan
model_cll <- paste0(folder, "\\Binomial_Reg_cloglog_stan.stan")

simullBSgs_inf_cll <-stan(file = model_cll, data = data, pars=pars,
                     chains = nchains, iter = niter,
                     warmup = nburnin, thin = nthin)
simullBSgs_inf_cll
summary(simullBSgs_inf_cll)

# Extraindo a cadeia (caso se queira usar outros
# pacotes e para fazer alguas contas)
mchain_cll <- (rstan::extract(simullBSgs_inf_cll,pars=pars,
                          permuted=FALSE,inc_warmup=FALSE))
result1_cll <- data.frame(mcmc(mchain_cll[,1,]))
str(result1_cll)


meu_tema <- theme_bw() + theme(text = element_text(size = 17))
# Estimativas básicas
mcmc_intervals(simullBSgs_inf_cll,pars=pars3) + meu_tema # IC p/ os betas


#pdf(paste(file.save,sep="","HistInf",nome_file,".pdf"))
plot_hist_cll <- mcmc_hist(simullBSgs_inf_cll,pars=pars3)+meu_tema # hist das posterioris
#dev.off()
#pdf(paste(file.save,sep="","DensInf",nome_file,".pdf"))
plot_dens_cll <- mcmc_dens(simullBSgs_inf_cll,pars=pars3)+meu_tema # dens das posterioris
#dev.off()

pdf("post__cll.pdf", width = 7, height = 4)
plot_hist_cll / plot_dens_cll
dev.off()

res_cll<-summary(simullBSgs_inf_cll)$summary

# Critérios de Informação e ajuste do modelo
# Checagem preditiva a posteriori
yrep_cll<-(rstan::extract(simullBSgs_inf_cll,
                      pars=c("yrep"))$yrep[order(sample(500)),])
y<-data$y

densmedia_cll <- ppc_stat(y,yrep_cll,stat="mean")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Média")

densvar_cll <- ppc_stat(y,yrep_cll,stat="var")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Variância")

densca_cll <- ppc_stat(y,yrep_cll,stat=function(y) skewness(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +labs(title = "Coef. de assimetria")

denscurt_cll <-ppc_stat(y,yrep_cll,stat=function(y) kurtosis(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+ labs(title = "Curtose")

pdf("momentos__cll.pdf", width = 10, height = 4)
densmedia_cll + densvar_cll + densca_cll + denscurt_cll
dev.off()

########Usando a proporção
m<-data$m
p<- y/m
prep_cll<- yrep_cll/t(matrix(m,length(y),500))

densprop_cll<-ppc_dens_overlay(p,prep_cll)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

cdfprop_cll<-ppc_ecdf_overlay(p,prep_cll)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

intervalprop_cll<- ppc_intervals(p, prep_cll)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ppc_ribbon(p, prep_cll)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

scatterprop_cll <- ppc_scatter_avg(p,prep_cll)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("prop_cll.pdf", width = 7, height = 4)
densprop_cll + cdfprop_cll + intervalprop_cll + scatterprop_cll
dev.off()

# LOO e WAIC (comparação de modelos)
m_log_lik_cll <-extract_log_lik(simullBSgs_inf_cll,
                            parameter_name="log_lik",
                            merge_chains = FALSE)
loo(m_log_lik_cll)
waic(m_log_lik_cll)

library(xtable)
# Estimativas
result_est_cll<-cbind(res_cll[1:length(pars3),c(1,6,3,4,8)],
                  matrix(t(apply(mchain_cll[,1,1:length(pars3)],2,emp.hpd,0.95)),length(pars3),2))
colnames(result_est_cll)[c(6,7)] <-c("HPD(2.5%)","HPD(97.5%)")
xtable(result_est_cll)
#






############## probito
# Parâmetros MCMC
nchains <- 1 # são 3 pra ver se converge
niter <- 20500 # número de iterações total
nburnin <- 500 # número de observações iniciais jogadas fora
nthin <- 10 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
# caminho do .stan
model_prob <- paste0(folder, "\\Binomial_Reg_probito_stan.stan")

simullBSgs_inf_prob <-stan(file = model_prob, data = data, pars=pars,
                          chains = nchains, iter = niter,
                          warmup = nburnin, thin = nthin)
simullBSgs_inf_prob
summary(simullBSgs_inf_prob)

# Extraindo a cadeia (caso se queira usar outros
# pacotes e para fazer alguas contas)
mchain_prob <- (rstan::extract(simullBSgs_inf_prob,pars=pars,
                              permuted=FALSE,inc_warmup=FALSE))
result1_prob <- data.frame(mcmc(mchain_prob[,1,]))
str(result1_prob)


meu_tema <- theme_bw() + theme(text = element_text(size = 17))
# Estimativas básicas
mcmc_intervals(simullBSgs_inf_prob,pars=pars3) + meu_tema # IC p/ os betas


#pdf(paste(file.save,sep="","HistInf",nome_file,".pdf"))
plot_hist_prob <- mcmc_hist(simullBSgs_inf_prob,pars=pars3)+meu_tema # hist das posterioris
#dev.off()
#pdf(paste(file.save,sep="","DensInf",nome_file,".pdf"))
plot_dens_prob <- mcmc_dens(simullBSgs_inf_prob,pars=pars3)+meu_tema # dens das posterioris
#dev.off()

pdf("post__prob.pdf", width = 7, height = 4)
plot_hist_prob / plot_dens_prob
dev.off()

res_prob<-summary(simullBSgs_inf_prob)$summary

# Critérios de Informação e ajuste do modelo
# Checagem preditiva a posteriori
yrep_prob<-(rstan::extract(simullBSgs_inf_prob,
                          pars=c("yrep"))$yrep[order(sample(500)),])
y<-data$y

densmedia_prob <- ppc_stat(y,yrep_prob,stat="mean")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Média")

densvar_prob <- ppc_stat(y,yrep_prob,stat="var")+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+labs(title = "Variância")

densca_prob <- ppc_stat(y,yrep_prob,stat=function(y) skewness(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +labs(title = "Coef. de assimetria")

denscurt_prob <-ppc_stat(y,yrep_prob,stat=function(y) kurtosis(y))+meu_tema+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+ labs(title = "Curtose")

pdf("momentos__prob.pdf", width = 10, height = 4)
densmedia_prob + densvar_prob + densca_prob + denscurt_prob
dev.off()

########Usando a proporção
m<-data$m
p<- y/m
prep_prob<- yrep_prob/t(matrix(m,length(y),500))

densprop_prob<-ppc_dens_overlay(p,prep_prob)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

cdfprop_prob<-ppc_ecdf_overlay(p,prep_prob)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

intervalprop_prob<- ppc_intervals(p, prep_prob)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ppc_ribbon(p, prep_prob)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

scatterprop_prob <- ppc_scatter_avg(p,prep_prob)+meu_tema#+
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("prop_prob.pdf", width = 7, height = 4)
densprop_prob + cdfprop_prob + intervalprop_prob + scatterprop_prob
dev.off()

# LOO e WAIC (comparação de modelos)
m_log_lik_prob <-extract_log_lik(simullBSgs_inf_prob,
                                parameter_name="log_lik",
                                merge_chains = FALSE)
loo(m_log_lik_prob)
waic(m_log_lik_prob)

library(xtable)
# Estimativas
result_est_prob<-cbind(res_prob[1:length(pars3),c(1,6,3,4,8)],
                      matrix(t(apply(mchain_prob[,1,1:length(pars3)],2,emp.hpd,0.95)),length(pars3),2))
colnames(result_est_prob)[c(6,7)] <-c("HPD(2.5%)","HPD(97.5%)")
xtable(result_est_prob)
#

result_est %>% xtable()
result_est_cll %>% xtable()
result_est_prob %>% xtable()


loo(m_log_lik)$estimates %>% xtable()
loo(m_log_lik_prob)$estimates %>% xtable()
loo(m_log_lik_cll)$estimates %>% xtable()

waic(m_log_lik)$estimates %>% xtable()
waic(m_log_lik_prob)$estimates %>% xtable()
waic(m_log_lik_cll)$estimates %>% xtable()

