library(coda)
library(bayesplot)
library(loo)
library(rstan)
# library(TeachingDemos)
library(moments)
library(tidyverse)

set.seed(40028922)
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
nchains <- 3 # são 3 pra ver se converge
niter <- 20500 # número de iterações total
nburnin <- 500 # número de observações iniciais jogadas fora
nthin <- 10 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
pars1 <- c("beta0","beta1","pest") # cadeias de interesse
pars3 <- c("beta0","beta1") # para rodar os gráficos
pars2 <- c("yrep","log_lik") # outras coisas de interesse
pars <- c(pars1,pars2)


# Execução do algoritmo MCMC

# caminho do .stan 
model<- file.path(folder, "Binomial_Reg_stan.stan")

simullBSgs_conv<-stan(file = model, data = data, pars=pars,
                 chains = nchains, iter = niter,
                 warmup = nburnin, thin = nthin)
simullBSgs_conv
summary(simullBSgs_conv)

# Extraindo as cadeias
mchain <- (rstan::extract(simullBSgs_conv,pars=pars3,
                          permuted=FALSE))
result1 <- data.frame(mcmc(mchain[,1,]))
result2 <- data.frame(mcmc(mchain[,2,]))
result3 <- data.frame(mcmc(mchain[,3,]))

# para usar, se for o caso, outros softwares MCMC
am_parametros <- data.frame("Cadeia1,beta0" = result1$beta0,
                      "Cadeia2,beta0" = result2$beta0,
                      "Cadeia3,beta0" = result3$beta0,
                      "Cadeia1,beta1" = result1$beta1,
                      "Cadeia2,beta1" = result2$beta1,
                      "Cadeia3,beta1" = result3$beta1, 
                      "idx" = 1:((niter - nburnin)/nthin))
am_parametros <- am_parametros %>% pivot_longer(cols = 1:6, names_to = "Cadeia", values_to = "Valor") %>% 
  separate("Cadeia", into = c("Cadeia", "Parâmetro"), sep = "\\.")

# Estudo da convergência
#pdf(paste(file.save,sep="","TracePlot",nome_file,".#pdf"))
mcmc_trace(simullBSgs_conv,pars=pars3)
#dev.off()
# Outra forma de ver esse gráfico

pdf("TracePlot.pdf", width = 7, height = 4)
am_parametros %>% ggplot(aes(x = idx, y = Valor, color = Cadeia))+
  geom_line()+
  facet_wrap(~Parâmetro, scales = "free")+
  theme_bw()
dev.off()

pdf("acf.pdf", width = 7, height = 4)
mcmc_acf(simullBSgs_conv,pars=pars3)
dev.off()

#pdf("violin.pdf")

mcmc_violin(simullBSgs_conv,pars=pars3)
#dev.off()

rhatMCMC <- rhat(simullBSgs_conv,pars=pars1)
rhatMCMC_aux <- rhatMCMC %>% data.frame()
colnames(rhatMCMC_aux) <- c("r")

pdf("r.pdf", width = 7, height = 4)
ggplot(rhatMCMC_aux)+
  geom_histogram(aes(x = r), color = 1, fill = "darkblue", bins = 10)+
  theme_bw()+
  geom_vline(xintercept = 1, color = "red", size = 1.5)+
  labs(x = "R", y = "Contagem")
dev.off()


pdf("dens.pdf", width = 7, height = 4)
mcmc_dens_overlay(simullBSgs_conv,pars=pars3)
dev.off()
#pdf(paste(file.save,sep="","Hist",nome_file,".#pdf"))
mcmc_hist_by_chain(simullBSgs_conv,pars=pars3)
#dev.off()




################## NÃO FAÇAM ISSO !!!
######### modelo cloglog

nchains <- 3 # são 3 pra ver se converge
niter <- 2010 # número de iterações total
nburnin <- 10 # número de observações iniciais jogadas fora
nthin <- 1 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
model_cloglog<- file.path(folder, "Binomial_Reg_cloglog_stan.stan")

simullBSgs_conv_err <-stan(file = model_cloglog, data = data,
                           pars=pars,
                           chains = nchains, iter = niter,
                           warmup = nburnin, thin = nthin)
simullBSgs_conv_err
summary(simullBSgs_conv_err)

# Extraindo as cadeias
mchain_err <- (rstan::extract(simullBSgs_conv_err,pars=pars3,
                              permuted=FALSE))
result1_err <- data.frame(mcmc(mchain_err[,1,]))
result2_err <- data.frame(mcmc(mchain_err[,2,]))
result3_err <- data.frame(mcmc(mchain_err[,3,]))

# para usar, se for o caso, outros softwares MCMC
am_parametros_err <- data.frame("Cadeia1,beta0" = result1_err$beta0,
                                "Cadeia2,beta0" = result2_err$beta0,
                                "Cadeia3,beta0" = result3_err$beta0,
                                "Cadeia1,beta1" = result1_err$beta1,
                                "Cadeia2,beta1" = result2_err$beta1,
                                "Cadeia3,beta1" = result3_err$beta1, 
                                "idx" = 1:((niter - nburnin)/nthin))
am_parametros_err <- am_parametros_err %>% pivot_longer(cols = 1:6, names_to = "Cadeia", values_to = "Valor") %>% 
  separate("Cadeia", into = c("Cadeia", "Parâmetro"), sep = "\\.")

# Estudo da convergência
#pdf(paste(file.save,sep="","TracePlot",nome_file,".#pdf"))
mcmc_trace(simullBSgs_conv_err,pars=pars3)
#dev.off()
# Outra forma de ver esse gráfico

pdf("TracePlot_err.pdf", width = 7, height = 4)
am_parametros_err %>% ggplot(aes(x = idx, y = Valor, color = Cadeia))+
  geom_line()+
  facet_wrap(~Parâmetro, scales = "free")+
  theme_bw()
dev.off()

pdf("acf_err.pdf", width = 7, height = 4)
mcmc_acf(simullBSgs_conv_err,pars=pars3)
dev.off()

#pdf(paste(file.save,sep="","Violin",nome_file,".#pdf"))
mcmc_violin(simullBSgs_conv_err,pars=pars3)
#dev.off()

rhatMCMC <- rhat(simullBSgs_conv_err,pars=pars1)
rhatMCMC_aux <- rhatMCMC %>% data.frame()
colnames(rhatMCMC_aux) <- c("r")

pdf("r_err.pdf", width = 7, height = 4)
ggplot(rhatMCMC_aux)+
  geom_histogram(aes(x = r), color = 1, fill = "darkblue", bins = 10)+
  theme_bw()+
  geom_vline(xintercept = 1, color = "red", size = 1.5)
dev.off()

pdf("dens_err.pdf", width = 7, height = 4)
mcmc_dens_overlay(simullBSgs_conv_err,pars=pars3)
dev.off()
#pdf(paste(file.save,sep="","Hist",nome_file,".#pdf"))
mcmc_hist_by_chain(simullBSgs_conv_err,pars=pars3)
#dev.off()





###########################################
######### modelo cloglog

nchains <- 3 # são 3 pra ver se converge
niter <- 210000 # número de iterações total
nburnin <- 10000 # número de observações iniciais jogadas fora
nthin <- 100 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
model_cloglog<- file.path(folder, "Binomial_Reg_cloglog_stan.stan")

simullBSgs_conv_cll <-stan(file = model_cloglog, data = data,
                           pars=pars,
                      chains = nchains, iter = niter,
                      warmup = nburnin, thin = nthin)
simullBSgs_conv_cll
summary(simullBSgs_conv_cll)

# Extraindo as cadeias
mchain_cll <- (rstan::extract(simullBSgs_conv_cll,pars=pars3,
                          permuted=FALSE))
result1_cll <- data.frame(mcmc(mchain_cll[,1,]))
result2_cll <- data.frame(mcmc(mchain_cll[,2,]))
result3_cll <- data.frame(mcmc(mchain_cll[,3,]))

# para usar, se for o caso, outros softwares MCMC
am_parametros_cll <- data.frame("Cadeia1,beta0" = result1_cll$beta0,
                            "Cadeia2,beta0" = result2_cll$beta0,
                            "Cadeia3,beta0" = result3_cll$beta0,
                            "Cadeia1,beta1" = result1_cll$beta1,
                            "Cadeia2,beta1" = result2_cll$beta1,
                            "Cadeia3,beta1" = result3_cll$beta1, 
                            "idx" = 1:((niter - nburnin)/nthin))
am_parametros_cll <- am_parametros_cll %>% pivot_longer(cols = 1:6, names_to = "Cadeia", values_to = "Valor") %>% 
  separate("Cadeia", into = c("Cadeia", "Parâmetro"), sep = "\\.")

# Estudo da convergência

mcmc_trace(simullBSgs_conv_cll,pars=pars3)

# Outra forma de ver esse gráfico
pdf("TracePlot_cll.pdf", width = 7, height = 4)
am_parametros_cll %>% ggplot(aes(x = idx, y = Valor, color = Cadeia))+
  geom_line()+
  facet_wrap(~Parâmetro, scales = "free")+
  theme_bw()
dev.off()

pdf("acf_cll.pdf", width = 7, height = 4)
mcmc_acf(simullBSgs_conv_cll,pars=pars3)
dev.off()

#pdf(paste(file.save,sep="","Violin",nome_file,".#pdf"))
mcmc_violin(simullBSgs_conv_cll,pars=pars3)
#dev.off()

rhatMCMC <- rhat(simullBSgs_conv_cll,pars=pars1)
rhatMCMC_aux <- rhatMCMC %>% data.frame()
colnames(rhatMCMC_aux) <- c("r")
pdf("r_cll.pdf", width = 7, height = 4)
ggplot(rhatMCMC_aux)+
  geom_histogram(aes(x = r), color = 1, fill = "darkblue", bins = 10)+
  theme_bw()+
  geom_vline(xintercept = 1, color = "red", size = 1.5)
dev.off()

pdf("dens_cll.pdf", width = 7, height = 4)
mcmc_dens_overlay(simullBSgs_conv_cll,pars=pars3)
dev.off()
#pdf(paste(file.save,sep="","Hist",nome_file,".#pdf"))
mcmc_hist_by_chain(simullBSgs_conv_cll,pars=pars3)
#dev.off()

######### modelo probito
# Parâmetros MCMC
nchains <- 3 # são 3 pra ver se converge
niter <- 20500 # número de iterações total
nburnin <- 500 # número de observações iniciais jogadas fora
nthin <- 10 # espaço entre observações consideradas, ou seja, ele pega de 10 em 10
model_probito<- file.path(folder, "Binomial_Reg_probito_stan.stan")

simullBSgs_conv_prob <-stan(file = model_probito, data = data,
                           pars=pars,
                           chains = nchains, iter = niter,
                           warmup = nburnin, thin = nthin)
simullBSgs_conv_prob
summary(simullBSgs_conv_prob)

# Extraindo as cadeias
mchain_prob <- (rstan::extract(simullBSgs_conv_prob,pars=pars3,
                          permuted=FALSE))
result1_prob <- data.frame(mcmc(mchain_prob[,1,]))
result2_prob <- data.frame(mcmc(mchain_prob[,2,]))
result3_prob <- data.frame(mcmc(mchain_prob[,3,]))

# para usar, se for o caso, outros softwares MCMC
am_parametros_prob <- data.frame("Cadeia1,beta0" = result1_prob$beta0,
                                "Cadeia2,beta0" = result2_prob$beta0,
                                "Cadeia3,beta0" = result3_prob$beta0,
                                "Cadeia1,beta1" = result1_prob$beta1,
                                "Cadeia2,beta1" = result2_prob$beta1,
                                "Cadeia3,beta1" = result3_prob$beta1, 
                                "idx" = 1:((niter - nburnin)/nthin))
am_parametros_prob <- am_parametros_prob %>% pivot_longer(cols = 1:6, names_to = "Cadeia", values_to = "Valor") %>% 
  separate("Cadeia", into = c("Cadeia", "Parâmetro"), sep = "\\.")

# Estudo da convergência
#pdf(paste(file.save,sep="","TracePlot",nome_file,".#pdf"))
mcmc_trace(simullBSgs_conv_prob,pars=pars3)
#dev.off()
# Outra forma de ver esse gráfico

pdf("TracePlot_prob.pdf", width = 7, height = 4)
am_parametros_prob %>% ggplot(aes(x = idx, y = Valor, color = Cadeia))+
  geom_line()+
  facet_wrap(~Parâmetro, scales = "free")+
  theme_bw()
dev.off()

pdf("acf_prob.pdf", width = 7, height = 4)
mcmc_acf(simullBSgs_conv_prob,pars=pars3)
dev.off()

#pdf(paste(file.save,sep="","Violin",nome_file,".#pdf"))
mcmc_violin(simullBSgs_conv_prob,pars=pars3)
#dev.off()

rhatMCMC <- rhat(simullBSgs_conv_prob,pars=pars1)
rhatMCMC_aux <- rhatMCMC %>% data.frame()
colnames(rhatMCMC_aux) <- c("r")

pdf("r_prob.pdf", width = 7, height = 4)
ggplot(rhatMCMC_aux)+
  geom_histogram(aes(x = r), color = 1, fill = "darkblue", bins = 10)+
  theme_bw()+
  geom_vline(xintercept = 1, color = "red", size = 1.5)
dev.off()

pdf("dens_prob.pdf", width = 7, height = 4)
mcmc_dens_overlay(simullBSgs_conv_prob,pars=pars3)
dev.off()
#pdf(paste(file.save,sep="","Hist",nome_file,".#pdf"))
mcmc_hist_by_chain(simullBSgs_conv_prob,pars=pars3)
#dev.off()

