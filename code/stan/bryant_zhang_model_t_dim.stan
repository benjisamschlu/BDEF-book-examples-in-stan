data {
  int<lower=0> A; // nber of age groups
  vector<lower=0>[A] n; // exposure
  int<lower=0> d[A]; // deaths
  vector<lower=0>[A] age0; // ind. variable for age 0
}



parameters {
  // Intercept
  real beta_0;
  
  // Age dimension
  vector[A] beta_a;
  vector[A] alpha_a;
  real phi;

  vector[A] gamma_a;

  real<lower=0> sigma_a;
  real<lower=0> sigma_aa;
  real<lower=0> sigma_aaa;
  
  // Do not predict with complete accuracy
  vector[A] e_a;
  real<lower=0> sigma;
}



transformed parameters {
  vector[A] log_lambda;
  log_lambda = beta_0 + beta_a + e_a + log(n);
}



model {
  // Likelihood
  d ~ poisson_log(log_lambda);



  // Priors 
  
  // Intercept
  beta_0 ~ normal(0, 10);
  
  // Age dimension
  
  beta_a ~ normal(alpha_a + phi*age0, sigma_a);
  
  phi ~ student_t(7, 0, 1);
  alpha_a[1] ~ normal(0, 10);
  gamma_a[1] ~ normal(0, 1);
  for (i in 2:A) {
      alpha_a[i] ~ normal(alpha_a[i-1] + gamma_a[i-1], sigma_aa);
      gamma_a[i] ~ normal(gamma_a[i-1], sigma_aaa);
  }
  sigma_a ~ student_t(7, 0, 1);
  sigma_aa ~ student_t(7, 0, 1);
  sigma_aaa ~ student_t(7, 0, 1);
  
  // Do not predict with complete accuracy
  e_a ~ normal(0, sigma);
  sigma ~ student_t(7, 0, 1);
}



generated quantities {
    vector[A] eta = log_lambda - log(n); // linear predictor
}