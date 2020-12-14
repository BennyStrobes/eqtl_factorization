functions {
  real folded_beta_binomial_log(int y, int n, real p, real conc) {
    real log_proby; 
    log_proby = log_sum_exp( beta_binomial_log(y, n, conc*p, conc*(1.0-p)), beta_binomial_log(y, n, conc*(1.0-p), conc*p));
    return log_proby;
  }
}
data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
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
    p[n] <- inv_logit(xb[n]);
    increment_log_prob(folded_beta_binomial_log(ys[n], ns[n], p[n], conc));
  }
  conc ~ gamma(concShape, concRate);
}
