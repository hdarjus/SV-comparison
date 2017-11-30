library(rstan)

panel.hist <- function(x, ...) {  # helper function
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
stan_pair <- function (stfit, i.chain) {
  pairs(cbind(as.matrix(extract(stfit, pars = c("phi", "sigma", "mu"),
                                inc_warmup = TRUE,
                                permuted = FALSE)[,i.chain,]),
              energy = get_sampler_params(stfit)[[i.chain]][, "energy__"]),  # energy has to be independent of model parameters (if I understand it correctly)
        diag.panel = panel.hist)                                             # so we check correlations here
}

stan_mytrace <- function (stfit) {
  sttrace_gg <- ggplot_build(stan_trace(stfit, pars = c("phi", "sigma", "mu")))
  sttrace_gg$layout$panel_layout$COL <- rep(1, 3)
  sttrace_gg$layout$panel_layout$ROW <- 1:3
  return(ggplot_gtable(sttrace_gg))
}

myvolplot <- function (fitmat, dat) {
  vol.quants <- t(apply(fitmat, 2, function (x) quantile(x, probs = c(2.5, 50, 97.5)/100)))
  rownames(vol.quants) <- as.character(1:N)
  plot(zoo(vol.quants), screen = 1, col = rainbow(NCOL(vol.quants)))
  lines(zoo(dat$vol))
}
