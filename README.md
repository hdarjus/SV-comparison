# SV-comparison
Comparison of different implementations of the same stochastic volatility model (`stochvol`, `JAGS`, `Stan`).
The goal is to see when and why one needs to develop a model-specific algorithm even if fast prototyping is easy using existing software packages.

The model in its standard parameterization is
```
y[t] ~ N(0,exp(h[t])),
h[t] ~ N(mu+phi*(h[t-1]-mu),sigma^2),
```
where, in order to be able to do Bayesian inference, one needs to put a prior distribution on the parameters `mu`, `sigma`, and `phi`, and on the initial latent state `h[0]`.

The R package `stochvol` uses several tricks to improve the sampling efficiency of the latent vector, most notably ASIS (ancillarity-sufficiency interweaving strategy) and AWOL (all without a loop).
`Stan` is the state-of-the-art probabilistic programming language, and `JAGS` is a language for graphical models.

The repository is still under development.

References:
- Kastner, Gregor. "Dealing with stochastic volatility in time series using the R package stochvol." Journal of Statistical software 69.5 (2016): 1-30.
- McCausland, William J., Shirley Miller, and Denis Pelletier. "Simulation smoothing for state–space models: A computational efficiency analysis." Computational Statistics & Data Analysis 55.1 (2011): 199-212.
- Yu, Yaming, and Xiao-Li Meng. "To center or not to center: That is not the question -- an Ancillarity–Sufficiency Interweaving Strategy (ASIS) for boosting MCMC efficiency." Journal of Computational and Graphical Statistics 20.3 (2011): 531-570.
