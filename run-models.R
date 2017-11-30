# Run models and save draws

library("rstan")
library("rjags")
library("coda")
library("stochvol")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

set.seed(42)

N <- 1000
phi <- 0.05
n.draws <- 2500  # warmup = adaptation + burnin is hard-coded 2500
sim <- svsim(N, mu = -10, phi = phi, sigma = 1, nu = Inf)
ret <- sim$y

dat <- list(y = ret, N = length(ret))

n.chains <- 3

# JAGS centered
jamodel.c <- jags.model("sv-c.bug", data = dat, n.chains = n.chains, n.adapt = 1000, inits = lapply(sample.int(.Machine$integer.max, n.chains), function (x) list(.RNG.seed = x, .RNG.name = "base::Wichmann-Hill")))
update(jamodel.c, 1500)
jafit.c <- coda.samples(jamodel.c, c("phi", "mu", "sigma", paste0("h[", 1:N, "]")), n.draws)
# JAGS non-centered
jamodel.nc <- jags.model("sv-nc.bug", data = dat, n.chains = n.chains, n.adapt = 1000, inits = lapply(sample.int(.Machine$integer.max, n.chains), function (x) list(.RNG.seed = x, .RNG.name = "base::Wichmann-Hill")))
update(jamodel.nc, 1500)
jafit.nc <- coda.samples(jamodel.nc, c("phi", "mu", "sigma", paste0("h[", 1:N, "]")), n.draws)

# Stan centered
stfit.c <- stan(file = 'sv-c.stan', data = dat, iter = n.draws + 2500, chains = n.chains, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = sample.int(.Machine$integer.max, 1))  # iter = 5000 ==> draws == burnin == 2500
# Stan non-centered
stfit.nc <- stan(file = 'sv-nc.stan', data = dat, iter = n.draws + 2500, chains = n.chains, control = list(adapt_delta = 0.95, max_treedepth = 15), seed = sample.int(.Machine$integer.max, 1))  # iter = 5000 ==> draws == burnin == 2500

# stochvol centered
svfit.c <- svsample(ret, n.draws, 2500, priormu = c(0,100), priorphi = c(1,1), priorsigma = 0.5, expert = list(parameterization = "centered"))
# stochvol non-centered
svfit.nc <- svsample(ret, n.draws, 2500, priormu = c(0,100), priorphi = c(1,1), priorsigma = 0.5, expert = list(parameterization = "noncentered"))
# stochvol ASIS
svfit.asis <- svsample(ret, n.draws, 2500, priormu = c(0,100), priorphi = c(1,1), priorsigma = 0.5, expert = list(parameterization = "GIS_C"))

# save everything
sphi <- as.character(phi)
filename <- paste0("phi-", substr(sphi, 3, nchar(sphi)), ".RDS")
saveRDS(list(
  data = sim,
  jamodel.c = jamodel.c,
  jamodel.nc = jamodel.nc,
  jafit.c = jafit.c,
  jafit.nc = jafit.nc,
  stfit.c = stfit.c,
  stfit.nc = stfit.nc,
  svfit.c = svfit.c,
  svfit.nc = svfit.nc,
  svfit.asis = svfit.asis
), filename)

rm(list = ls())
