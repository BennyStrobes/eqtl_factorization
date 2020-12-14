data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0,upper=1> prob1[N];
  real<lower=0,upper=2> prob2[N];
  real intercept[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  vector[P] beta;
  real<lower=0> conc;
}
model {
  vector[N] xb; 
  real p[N]; 
  xb <- x * beta;
  for (n in 1:N) {
    p[n] <- inv_logit(xb[n] + intercept[n]);
    increment_log_prob((prob1[n])*beta_binomial_log(ys[n], ns[n], conc*(p[n]), conc*(1.0-(p[n]))));
    increment_log_prob((prob2[n])*beta_binomial_log(ys[n], ns[n], conc*(1.0-p[n]), conc*(p[n])));

  }
  conc ~ gamma(concShape, concRate);
}
