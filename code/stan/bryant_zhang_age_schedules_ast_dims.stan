data {
  int<lower=0> A;  // nber of age groups
  int<lower=0> T;  // nber of years
  int<lower=0> S;  // nber of genders, s==1 is female

  int D[S,A,T];  // deaths: 2-array of deaths by age group,year
  matrix[A,T]  E[S];  // exposure: 2-array of person-yrs by age group,year

  vector<lower=0>[A] age0; // indicator variable for age 0
  vector[S] sex; //indicator variable for sex
}
parameters {
  
  // Model does not predict with complete accuracy
  matrix[A,T] u[S];
  real<lower=0> sigma;
  
  
  // Intercept
  real beta_0;
  
  
  // Age dimension
  vector[A] beta_a;
  
  vector[A] alpha_a;
  real psi;

  vector[A] gamma_a;

  real<lower=0> sigma_a;
  real<lower=0> sigma_aa;
  real<lower=0> sigma_aaa;
  
  
  // Sex dimension
  real beta_s; // beta_0 for female and beta_s diff for male
  
  
  // Time dimension
  row_vector[T] beta_t;
  
  real alpha_t0;
  
  vector[T] gamma_t;

  real<lower=0> sigma_t;
  real<lower=0> sigma_tt;
  
  
  // Age-sex dimension
  vector[A] beta_as;
  
  vector[A] alpha_as;
  
  real<lower=0> sigma_as;
  
  
  // Age-time dimension
  matrix[A,T] beta_at;
  
  vector[A] alpha_at0;
  real<lower=0.8, upper=1> phi_at; // rescaled to lie in 0.8-1
  
  matrix[A,T] gamma_at;

  real<lower=0> sigma_at;
  real<lower=0> sigma_aatt;
  
  
  // Sex-time dimension
  matrix[S,T] beta_st;
  
  vector[S] alpha_st0;
  real<lower=0.8, upper=1> phi_st; // rescaled to lie in 0.8-1
  
  matrix[S,T] gamma_st;

  real<lower=0> sigma_st;
  real<lower=0> sigma_sstt;
}



transformed parameters {
  matrix[A,T] eta[S] ; // linear predictor
  vector[T] alpha_t;
  matrix[A,T] alpha_at;
  matrix[S,T] alpha_st;

  // eta maybe in generated quantities block ?
  for (s in 1:2) { 
    eta[s] = rep_matrix(rep_vector(beta_0, A), T) + 
             rep_matrix(beta_a, T) + 
             rep_matrix(rep_vector(sex[s]*beta_s, A), T) +  
             rep_matrix(beta_t, A) + 
             rep_matrix(sex[s]*beta_as, T) +
             beta_at +
             rep_matrix(beta_st[s,], A) +
             u[s];
  }
  
  // Time dimension
  alpha_t[1] = alpha_t0;
  for (t in 2:T) {  //Could try to replace all these loops by [2:T] = [1:(T-1)] + [2:T]
    alpha_t[t] = alpha_t[t-1] + gamma_t[t];
  }
  
  // Age-time dimension
  alpha_at[,1] = alpha_at0;
  for (t in 2:T) {
    alpha_at[,t] = alpha_at[,t-1] + gamma_at[,t];
  }
  
  // Sex-time dimension
  alpha_st[,1] = alpha_st0;
  for (t in 2:T) {
    alpha_st[,t] = alpha_st[,t-1] + gamma_st[,t];
  }
  
}



model {
  matrix[A,T] log_lambda[S];
  
  for (s in 1:2) {
    log_lambda[s] = rep_matrix(rep_vector(beta_0, A), T) +
                    rep_matrix(beta_a, T) + 
                    rep_matrix(rep_vector(sex[s]*beta_s, A), T) +  
                    rep_matrix(beta_t, A) + 
                    rep_matrix(sex[s]*beta_as, T) +
                    beta_at +
                    rep_matrix(beta_st[s,], A) +
                    u[s] +
                    log(E[s]);
  }
  
  
  // Likelihood //
  // ---------- //
  
  for (s in 1:2) {
    to_array_1d(D[s]) ~ poisson_log( to_array_1d( log_lambda[s]') ) ;
  }



  // Priors //
  // ------ //
  
  // Model does not predict with complete accuracy
  for (s in 1:2) {
     to_array_1d(u[s]) ~ normal(0, sigma);
  }
  sigma ~ student_t(7, 0, 1);

  
  // Intercept
  beta_0 ~ normal(0, 10);
  
  
  // Age dimension
  beta_a ~ normal(alpha_a + psi*age0, sigma_a);
  
  psi ~ student_t(7, 0, 1);
  alpha_a[1] ~ normal(0, 10);
  gamma_a[1] ~ normal(0, 1);
  for (a in 2:A) {
      alpha_a[a] ~ normal(alpha_a[a-1] + gamma_a[a], sigma_aa);
      gamma_a[a] ~ normal(gamma_a[a-1], sigma_aaa);
  }
  sigma_a ~ student_t(7, 0, 1);
  sigma_aa ~ student_t(7, 0, 1);
  sigma_aaa ~ student_t(7, 0, 1);
  
  
  // Sex dimension
  beta_s ~ normal(0, 1);
  
  
  // Time dimension
  beta_t' ~ normal(alpha_t, sigma_t);
  
  alpha_t0 ~ normal(0, 10);
  gamma_t[1] ~ normal(0, 1);
  for (t in 2:T) {
      gamma_t[t] ~ normal(gamma_t[t-1], sigma_tt);
  }
  sigma_t ~ student_t(7, 0, 1);
  sigma_tt ~ student_t(7, 0, 1);
  
  
  // Age-sex dimension
  beta_as ~ normal(alpha_as, sigma_as);
  
  alpha_as[1] ~ normal(0, 10);
  for (a in 2:A) {
    alpha_as[a] ~ normal(alpha_as[a-1], sigma_as);
  }
  sigma_as ~ student_t(7, 0, 1);
  
  
  // Age-time dimension
  to_array_1d(beta_at) ~ normal(to_array_1d(alpha_at), sigma_at); 
  
  alpha_at0 ~ normal(0, 10);
  gamma_at[, 1] ~ normal(0, 1);
  for (t in 2:T) {
      gamma_at[,t] ~ normal(phi_at*gamma_at[,t-1], sigma_aatt);
      }
  phi_at ~ beta(2,2);
  
  sigma_at ~ student_t(7, 0, 1);
  sigma_aatt ~ student_t(7, 0, 1);
  
  
  // Sex-time dimension
  to_array_1d(beta_st) ~ normal(to_array_1d(alpha_st), sigma_st); 
  
  alpha_st0 ~ normal(0, 10);
  gamma_st[, 1] ~ normal(0, 1);
  for (t in 2:T) {
      gamma_st[,t] ~ normal(phi_st*gamma_st[,t-1], sigma_sstt);
      }
  phi_st ~ beta(2,2);
  
  sigma_st ~ student_t(7, 0, 1);
  sigma_sstt ~ student_t(7, 0, 1);

}

