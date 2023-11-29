set.seed(2023)
n <- 1000
  x1 <- sample(c(0,1), size = n, replace = TRUE)
x2 <- round(runif(n, 18, 80))
eta <- -9 + 3.5*x1 + 0.2*x2
pi <- exp(eta)/(1 + exp(eta))
summary(pi)
y <- rbinom(n = n, size = 1, prob = pi)
initial <- matrix(c(-9, 3.5, 0.2), ncol = 1)
dta <- data.frame(x1 = x1, x2 = x2, y = y)
fit.glm <- glm(y ~ x1 + x2, family = binomial(link = "logit"), data = dta)
summary(fit.glm)
BIC(fit.glm)
fit.mlg <- mlg(y ~ x1 + x2, dados = dta, init = initial, distr = "Binomial")


x1 <- sample(c(0,1), size = n, replace = TRUE)
x2 <- round(runif(n, 18, 20))
eta <- -5 + 3.5*x1 + 0.2*x2
pi <- pnorm(eta)
summary(pi)
y <- rbinom(n = n, size = 1, prob = pi)
initial <- matrix(c(-5, 3.5, 0.2), ncol = 1)
dta <- data.frame(x1 = x1, x2 = x2, y = y)
fit.glm2 <- glm(y ~ x1 + x2, family = binomial(link = "probit"), data = dta)
summary(fit.glm2)
fit.mlg2 <- mlg(y ~ x1 + x2, dados = dta, init = initial, distr = "Binomial", link = "probit")
fit.mlg2


x1 <- sample(c(0,1), size = n, replace = TRUE)
x2 <- round(runif(n, 18, 20))
eta <- -5 + 1.5*x1 + 0.2*x2
pi <- 1-exp(-exp(eta))
summary(pi)
y <- rbinom(n = n, size = 1, prob = pi)
initial <- matrix(c(-5, 1.5, 0.2), ncol = 1)
dta <- data.frame(x1 = x1, x2 = x2, y = y)
fit.glm3 <- glm(y ~ x1 + x2, family = binomial(link = "cloglog"), data = dta)
summary(fit.glm3)
fit.mlg3 <- mlg(y ~ x1 + x2, dados = dta, init = initial, distr = "Binomial", link = "cloglog")
fit.mlg3

