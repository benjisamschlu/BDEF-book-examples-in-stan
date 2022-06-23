
// Model from BDEF for Infant mortality rates in Sweden, 1995-2015


data {
  int<lower=0> T; // nber of years
  int<lower=0> C; // nber of counties

  int<lower=0> y[C,T]; // deaths
  int<lower=0> w[C,T]; // exposure
  
  int<lower=0> nprojyears; // projection from 2016-2025

}


parameters {
  matrix[C,T] eta;
  
  real beta_0;
  vector[C] beta_c;
  row_vector[T] beta_t;
  real<lower=0> sigma;
  
  real<lower=0> tau_c;
  
  real<lower=0> tau_t;
  vector[T] alpha;
  vector[T] delta;
  real<lower=0> omega_alpha;
  real<lower=0> omega_delta;



}


transformed parameters {
  matrix[C,T] mu;
  matrix[C,T] gamma;
  mu = rep_matrix(rep_vector(beta_0, C), T) + 
        rep_matrix(beta_c, T) + 
        rep_matrix(beta_t, C);
        
  gamma = inv_logit(eta);
}


model {
  
  // Likelihood
  for (t in 1:T) {
      target += binomial_logit_lpmf(y[,t] | w[,t], eta[,t]);
  }

  
  // Priors
  
  // Model does not fit perfectly
  for (t in 1:T) {
      target += normal_lpdf(eta[,t] | mu[,t], sigma);
  }
  // Intercept
  target += normal_lpdf(beta_0 | 0, 10);
  
  // County effect
  target += normal_lpdf(beta_c | 0, tau_c);
  
  // Time effect: local trend model
  target += normal_lpdf(beta_t' | alpha, tau_t);
  
  target += normal_lpdf(alpha[1] | 0, omega_alpha);
  target += normal_lpdf(delta[1] | 0, omega_delta);
  target += normal_lpdf(alpha[2:T] | alpha[1:(T-1)] + delta[2:T], omega_alpha);
  target += normal_lpdf(delta[2:T] | delta[1:(T-1)], omega_delta);
  
  // Standard deviations
  target += student_t_lpdf(sigma | 7, 0, 1);
  target += student_t_lpdf(tau_c | 7, 0, 1);
  target += student_t_lpdf(tau_t | 7, 0, 1);
  target += student_t_lpdf(omega_alpha | 7, 0, 1);
  target += student_t_lpdf(omega_delta | 7, 0, 1);
}



// Forecast 2016-2025
generated quantities {
  vector[nprojyears] alpha_proj;
  vector[nprojyears] delta_proj;
  row_vector[nprojyears] beta_t_proj;
  matrix[C,nprojyears] mu_proj;
  matrix[C,nprojyears] eta_proj;
  matrix[C,nprojyears] gamma_proj;

  
  delta_proj[1] = normal_rng(delta[T], omega_delta);
  alpha_proj[1] = normal_rng(alpha[T] + delta_proj[1], omega_alpha);
  for (p in 2:nprojyears) {
    delta_proj[p] = normal_rng(delta_proj[p-1], omega_delta);
    alpha_proj[p] = normal_rng(alpha_proj[p-1] + delta_proj[p], omega_alpha);
  }
  for (p in 1:nprojyears) {
    beta_t_proj[p] = normal_rng(alpha_proj[p], tau_t);
  }
  
  mu_proj = rep_matrix(rep_vector(beta_0, C), nprojyears) + 
            rep_matrix(beta_c, nprojyears) + 
            rep_matrix(beta_t_proj, C);
  for (p in 1:nprojyears) {
    for (c in 1:C) {
      eta_proj[c,p] = normal_rng(mu_proj[c,p], sigma);
    }
  }
  gamma_proj = inv_logit(eta_proj);
}