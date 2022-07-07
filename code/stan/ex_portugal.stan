data {
  int<lower=0> n; // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> T;  // nber of years
  int<lower=0> S;  // nber of genders, s==1 is female
  int<lower=0> Tproj; // nber of projected years 

  int D[n];  // deaths: order is age,year,sex
  vector[n] E;  // exposure: order is age,year,sex

  vector<lower=0>[A] age0; // indicator variable for age 0
  vector[S] sex; //indicator variable for sex
}

parameters {
  //-- MODEL DOES NOT MODEL PERFECTLY --//
  matrix[A,T] u[S];
  real<lower=0> sigma;
  
  real beta0; // Intercept
  
  //-- MAIN EFFECTS --//
  //-- Age
  vector[A] beta_a;
  real<lower=0> sigma_a;
  
  real psi; 
  vector[A] alpha_a;
  vector[A] gamma_a;
  
  real<lower=0> omega_a;
  real<lower=0> omega_aa;

  //-- Time
  vector[T] beta_t;
  real<lower=0> sigma_t;

  real alpha_t0;
  vector[T] gamma_t;
  real<lower=0> omega_t;
  
  //-- Sex 
  real beta_s;
  
  //-- INTERACTIONS EFFECTS --//
  //-- Age-time
  matrix[A,T] beta_at;
  real<lower=0> sigma_at;

  vector[A] alpha_at0;
  real<lower=0.8, upper=1> phi_at;
  matrix[A,T] gamma_at;
  real<lower=0> omega_at;
  
  //-- Sex-time
  vector[T] beta_st;
  real<lower=0> sigma_st;

  real alpha_st0;
  real<lower=0.8, upper=1> phi_st;
  vector[T] gamma_st;
  real<lower=0> omega_st;
  
  //-- Age-sex 
  vector[A] beta_as;
  real<lower=0> sigma_as;

  vector[A] alpha_as;
  real<lower=0> omega_as;
}

transformed parameters {
  vector[T] alpha_t;
  matrix[A,T] alpha_at;
  vector[T] alpha_st;

  matrix[A,T] eta[S];

  //-- MAIN EFFECTS --//
  //-- Age 
  
  //-- Time
  alpha_t[1] = alpha_t0;
for (t in 2:T) {  // Could try to replace all these loops by [2:T] = [1:(T-1)] + [2:T]
    alpha_t[t] = alpha_t[t-1] + gamma_t[t];
  }  
  //-- Age-time
  alpha_at[, 1] = alpha_at0;
  for (t in 2:T) {
    alpha_at[,t] = alpha_at[,t-1] + gamma_at[,t];
  }
  
  //-- Sex-time
  alpha_st[1] = alpha_st0;
  for (t in 2:T) {
    alpha_st[t] = alpha_st[t-1] + gamma_st[t];
  }  

  for (s in 1:S) {
    eta[s] = rep_matrix(rep_vector(beta0, A), T) +
                    rep_matrix(beta_a, T) + 
                    rep_matrix(rep_vector(sex[s]*beta_s, A), T) +  
                    rep_matrix(beta_t', A) + 
                    rep_matrix(sex[s] * beta_as, T) +
                    beta_at +
                    rep_matrix(sex[s] * beta_st', A) +
                    u[s] * sigma; // model does not predict perfectly
  }
}  

model {
  //-- SAMPLE LIKELIHOOD --//
  {
    matrix[A*T, S] eta_ord;
    
    for (s in 1:S) {
      eta_ord[, s] = to_vector(eta[s]); // age,time,sex
    }
  target += poisson_log_lpmf(D | to_vector(eta_ord) + log(E));
  }



  //-- PRIORS --//
  
  for (s in 1:S) {
    to_vector(u[s]) ~ std_normal();
  }
  sigma ~ student_t(7, 0, 1);
  
  //-- Intercept
  beta0 ~ normal(0, 1);
  
  //-- Age 
  beta_a ~ normal(alpha_a + psi * age0, sigma_a);
  psi ~ student_t(7, 0, 1);
  sigma_a ~ student_t(7, 0, 1);

  alpha_a[2:A] ~ normal(alpha_a[1:(A-1)] + gamma_a[2:A], omega_a);
  gamma_a[2:A] ~ normal(gamma_a[1:(A-1)], omega_aa);
  omega_a ~ student_t(7, 0, 1);
  omega_aa ~ student_t(7, 0, 1);
  
  //-- Sex 
  beta_s ~ normal(0, 1);
  
  //-- Time 
  beta_t ~ normal(alpha_t, sigma_t);
  sigma_t ~ student_t(7, 0, 1);
  
  alpha_t0 ~ normal(0, 1);
  gamma_t[2:T] ~ normal(gamma_t[1:(T-1)], omega_t);
  omega_t ~ student_t(7, 0, 1);

  // Age-time
  to_vector(beta_at) ~ normal(to_vector(alpha_at), sigma_at);
  sigma_at ~ student_t(7, 0, 1);
  
  alpha_at0 ~ normal(0, 1);
  phi_at ~ beta(2, 2);
  for (a in 1:A) {
    gamma_at[a, 2:T] ~ normal(phi_at * gamma_at[a, 1:(T-1)], omega_at);
  }
  omega_at ~ student_t(7, 0, 1);
  
  // Sex-time
  beta_st ~ normal(alpha_st, sigma_st);
  sigma_st ~ student_t(7, 0, 1);
  
  alpha_st0 ~ normal(0, 1);
  phi_st ~ beta(2, 2);
  gamma_st[2:T] ~ normal(phi_st * gamma_st[1:(T-1)], omega_st);
  omega_st ~ student_t(7, 0, 1);

  // Age-sex 
  beta_as ~ normal(alpha_as, sigma_as);
  sigma_as ~ student_t(7, 0, 1);
  
  alpha_as[1] ~ normal(0, 1);
  alpha_as[2:A] ~ normal(alpha_as[1:(A-1)], omega_as);
  omega_as ~ student_t(7, 0, 1);
}

// Forecast 1996-2015
generated quantities {
  vector[Tproj] gamma_t_proj;
  vector[Tproj] alpha_t_proj;
  vector[Tproj] beta_t_proj;
  
  matrix[A, Tproj] gamma_at_proj;
  matrix[A, Tproj] alpha_at_proj;
  matrix[A, Tproj] beta_at_proj;
  
  vector[Tproj] gamma_st_proj;
  vector[Tproj] alpha_st_proj;
  vector[Tproj] beta_st_proj;
  
  matrix[A,Tproj] eta_proj[S];
  matrix[A,Tproj] u_proj[S];


  // time
  gamma_t_proj[1] = normal_rng(gamma_t[T], omega_t);
  alpha_t_proj[1] = alpha_t[T] + gamma_t_proj[1];
  for (p in 2:Tproj) {
    gamma_t_proj[p] = normal_rng(gamma_t_proj[p-1], omega_t);
    alpha_t_proj[p] = alpha_t_proj[p-1] + gamma_t_proj[p];
  }
  for (p in 1:Tproj) {
    beta_t_proj[p] = normal_rng(alpha_t_proj[p], sigma_t);
  }
  
  // age-time
  for (a in 1:A) {
    gamma_at_proj[a, 1] = normal_rng(phi_at * gamma_at[a, T], omega_at);
    alpha_at_proj[a, 1] = alpha_at[a, T] + gamma_at_proj[a, 1];
  }
  
  for (p in 2:Tproj) {
    for (a in 1:A) {
      gamma_at_proj[a, p] = normal_rng(phi_at * gamma_at_proj[a, p-1], omega_at);
      alpha_at_proj[a, p] = alpha_at_proj[a, p-1] + gamma_at_proj[a, p];
    }
  }
  for (p in 1:Tproj) {
    for (a in 1:A) {
      beta_at_proj[a, p] = normal_rng(alpha_at_proj[a, p], sigma_at);
    }
  }
  
  // sex-time
  gamma_st_proj[1] = normal_rng(gamma_st[T], omega_st);
  alpha_st_proj[1] = alpha_st[T] + gamma_st_proj[1];
  for (p in 2:Tproj) {
    gamma_st_proj[p] = normal_rng(phi_st * gamma_st_proj[p-1], omega_st);
    alpha_st_proj[p] = alpha_st_proj[p-1] + gamma_st_proj[p];
  }
  for (p in 1:Tproj) {
    beta_st_proj[p] = normal_rng(alpha_st_proj[p], sigma_st);
  }
  for (s in 1:S) {
    for (a in 1:A) {
      for (t in 1:Tproj) {
        u_proj[s][a,t] = normal_rng(0,1);

      }
    }
  }
  for (s in 1:S) {
    eta_proj[s] = rep_matrix(rep_vector(beta0, A), Tproj) +
                    rep_matrix(beta_a, Tproj) + 
                    rep_matrix(rep_vector(sex[s]*beta_s, A), Tproj) +  
                    rep_matrix(beta_t_proj', A) + 
                    rep_matrix(sex[s] * beta_as, Tproj) +
                    beta_at_proj +
                    rep_matrix(sex[s] * beta_st_proj', A) +
                    u_proj[s] * sigma; // model does not predict perfectly
  }
}