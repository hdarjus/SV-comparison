library("rstan")
library("rjags")
library("coda")
library("stochvol")
library("zoo")

source("utils.R")

# Read dataset
everything <- readRDS("phi-05.RDS")
jafit.c <- everything$jafit.c
stfit.c <- everything$stfit.c
svfit.c <- everything$svfit.c
jafit.nc <- everything$jafit.nc
stfit.nc <- everything$stfit.nc
svfit.nc <- everything$svfit.nc
svfit.asis <- everything$svfit.asis
dat <- everything$data
n.draws <- NROW(jafit.nc[[1]])
N <- length(dat$y)
i.chain <- 3

# Stan diagnostics
# Look for dependence between energy and the parameters. They should be independent
stan_pair(stfit.c, i.chain)
stan_pair(stfit.nc, i.chain)

# Parameter plots
pars <- c("phi", "sigma", "mu")
japar.c <- jafit.c[[i.chain]][, pars]
japar.nc <- jafit.nc[[i.chain]][, pars]
stpar.c <- as.matrix(extract(stfit.c, pars, inc_warmup = FALSE, permuted = FALSE)[,i.chain,])[, pars]
stpar.nc <- as.matrix(extract(stfit.nc, pars, inc_warmup = FALSE, permuted = FALSE)[,i.chain,])[, pars]

plot(zoo(japar.c), screen = 1:3)
plot(zoo(japar.nc), screen = 1:3)
plot(stan_mytrace(stfit.c))
plot(stan_mytrace(stfit.nc))
stan_dens(stfit.c, pars = pars)
stan_dens(stfit.nc, pars = pars)
stan_diag(stfit.c, info = "sample", chain = i.chain)
stan_diag(stfit.nc, info = "sample", chain = i.chain)
plot(svfit.c)
plot(svfit.nc)
plot(svfit.asis)

# Parameters' inefficiency factors
as.data.frame(n.draws / rbind(
  JAGS.c = effectiveSize(japar.c),
  JAGS.nc = effectiveSize(japar.nc),
  Stan.c = sapply(c("phi" = "phi", "sigma" = "sigma", "mu" = "mu"), function (x) effectiveSize(unlist(extract(stfit.c, x, inc_warmup = FALSE, permuted = FALSE)[,i.chain,]))),
  Stan.nc = sapply(c("phi" = "phi", "sigma" = "sigma", "mu" = "mu"), function (x) effectiveSize(unlist(extract(stfit.nc, x, inc_warmup = FALSE, permuted = FALSE)[,i.chain,]))),
  stochvol.c = effectiveSize(svfit.c$para),
  stochvol.nc = effectiveSize(svfit.nc$para),
  stochvol.asis = effectiveSize(svfit.asis$para)
))

# Latent vector plots
latent.names <- paste0("h[", 1:N, "]")
jalatent.c <- as.matrix(jafit.c[[i.chain]][, latent.names])
jalatent.nc <- apply(as.matrix(jafit.nc[[i.chain]][, latent.names]), 2, function (x, mu, sigma) { sigma*x + mu }, japar.nc[, "mu"], japar.nc[, "sigma"])
stlatent.c <- as.matrix(as.data.frame(extract(stfit.c, pars = latent.names, permuted = FALSE)[,i.chain,]))
stlatent.nc <- apply(as.matrix(as.data.frame(extract(stfit.nc, pars = latent.names, permuted = FALSE)[,i.chain,])), 2, function (x, mu, sigma) { sigma*x + mu }, stpar.nc[, "mu"], stpar.nc[, "sigma"])

myvolplot(exp(jalatent.c/2), dat)
myvolplot(exp(jalatent.nc/2), dat)
myvolplot(exp(stlatent.c/2), dat)
myvolplot(exp(stlatent.nc/2), dat)
myvolplot(as.matrix(exp(svfit.c$latent/2)), dat)
myvolplot(as.matrix(exp(svfit.nc$latent/2)), dat)
myvolplot(as.matrix(exp(svfit.asis$latent/2)), dat)

# Latent vector's inefficiency factors
as.data.frame(rbind(
  JAGS.c = summary(n.draws / effectiveSize(jalatent.c)),
  JAGS.nc = summary(n.draws / effectiveSize(jalatent.nc)),
  Stan.c = summary(n.draws / apply(stlatent.c, 2, effectiveSize)),
  Stan.nc = summary(n.draws / apply(stlatent.nc, 2, effectiveSize)),
  stochvol.c = summary(n.draws / effectiveSize(svfit.c$latent)),
  stochvol.nc = summary(n.draws / effectiveSize(svfit.nc$latent)),
  stochvol.asis = summary(n.draws / effectiveSize(svfit.asis$latent))
))
