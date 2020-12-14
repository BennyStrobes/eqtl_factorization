data {
  int<lower=0> N; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0,upper=1> prob1[N];
  real<lower=0> cauchy_scale;
}
parameters {
  real<lower=0> a1;
}
model {
  for (n in 1:N) {
    increment_log_prob((prob1[n])*beta_binomial_log(ys[n], ns[n], 1.0, 7.0 + a1));
  }
  a1 ~ cauchy(0, cauchy_scale);
}
