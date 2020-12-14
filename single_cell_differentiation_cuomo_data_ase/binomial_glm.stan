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
  }
  ys ~ binomial(ns, p);
}