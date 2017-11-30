data {
  int<lower=0> N;
  real y[N];
}

transformed data {
  real<lower=0> phi_prior_a = 1;
  real<lower=0> phi_prior_b = 1;
  real mu_prior_mu = 0;
  real<lower=0> mu_prior_sigma = 100;
  real<lower=0> sigma_prior_shape = 0.5;
  real<lower=0> sigma_prior_rate = 0.5;
}

parameters {
  real h[N];
  real<lower=0,upper=1> phi_beta;
  real mu;
  real<lower=0> sigma2;
}

transformed parameters {
  real<lower=-1,upper=1> phi;
  real<lower=0> sigma;
  real exph[N];
  real h_i_mean[N];
  real<lower=0> h_i_sigma[N];
  
  phi = fma(phi_beta, 2, -1);
  sigma = sqrt(sigma2);
  for (i in 1:N) {
    if (i == 1) {
      h_i_mean[i] = mu;
      h_i_sigma[i] = sigma * pow(1-square(phi), -0.5);
    } else {
      h_i_mean[i] = fma(phi, h[i-1], fma(-phi, mu, mu));
      h_i_sigma[i] = sigma;
    }
    exph[i] = exp(fma(h[i], 0.5, 0));
  }
}

model {
  mu ~ normal(mu_prior_mu, mu_prior_sigma);
  phi_beta ~ beta(phi_prior_a, phi_prior_b);
  sigma2 ~ gamma(sigma_prior_shape, sigma_prior_rate);
  
  for (i in 1:N) {
    h[i] ~ normal(h_i_mean[i], h_i_sigma[i]);
    y[i] ~ normal(0, exph[i]);
  }
}
