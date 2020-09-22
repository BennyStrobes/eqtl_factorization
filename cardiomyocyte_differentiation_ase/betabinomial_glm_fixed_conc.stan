data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0> conc[N]; 
  real intercept[N]; 
}
parameters {
  vector[P] beta;
}
model {
  vector[N] xb; 
  real a[N];
  real b[N];
  real p[N]; 
  xb <- x * beta;
  for (n in 1:N) {
    p[n] <- inv_logit(xb[n] + intercept[n]); 
    a[n] <- conc[n]*p[n];
    b[n] <- conc[n]*(1.0-p[n]);
  }
  // beta ~ normal(0,5);
  ys ~ beta_binomial(ns, a, b);
}