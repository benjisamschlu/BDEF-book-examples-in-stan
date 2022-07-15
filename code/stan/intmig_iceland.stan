functions {
  real customProb_lpmf(int x, int y) {
    if (x == y) { 
      
      return(log(1));
      
    } else if (abs(x-y) == 1) {
      
      return(log(0.6666667));
      
    } else {
      
      return(log(0.3333333));
      
    }
  }
}

data {
  int<lower=0> n; // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> T;  // nber of years
  int<lower=0> S;  // nber of genders, s==1 is female
  int<lower=0> R; // nber of regions 

  int X[n];  // confidentialized counts: order is ...
  vector[n] E;  // exposure: order is ...

  vector[S] sex; //indicator variable for sex
  
  vector[n] L; // lower bounds for y
  vector[n] U; // upper bounds for y
  
}

parameters {
  # int<lower=L, upper=U> y[n];
  
  real beta0;
  real beta_s;
  
  
}



model {}