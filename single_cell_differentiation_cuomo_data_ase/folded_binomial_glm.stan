functions {
  real folded_binomial_log(int y, int n, real p) {
    real log_proby; 
    log_proby = lchoose(n, y) + log_sum_exp( y*log(p) + (n-y)*log(1.0-p), (n-y)*log(p) + y*log(1.0-p));
    if ((y+y)==n) {
      log_proby = log_proby + log(0.5);
    }
    return log_proby;
  }
}
data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
}
parameters {
  vector[P] beta;
}
model {
  vector[N] xb; 
  real p[N]; 
  xb <- x * beta;
  for (n in 1:N) {
    p[n] <- inv_logit(xb[n]);
    increment_log_prob(folded_binomial_log(ys[n], ns[n], p[n]));
  }
}