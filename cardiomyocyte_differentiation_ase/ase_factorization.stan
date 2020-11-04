data {
  int<lower=0> N; // Number of cells
  int<lower=0> T; // Number of allelic sites
  int<lower=0> K; // Number of latent factors
  int<lower=0> num_cov; // Number of covariates
  matrix[N,num_cov] cov; 
  int<lower=0> ys[T, N];
  int<lower=0> ns[T, N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
}
parameters {
  real<lower=0> conc[T]; // Concentration parameter for each site
  matrix[N, K] U;  // Low dimensional representation of cells
  matrix[K, T] V;  // Low dimensional representation of test effect sizes
  matrix[num_cov, T] C;  // Fixed effects of covariates
}
model {
  matrix[N,T] mu = U*V + cov*C;
  for (t in 1:T) {
    real a[N];
    real b[N];
    real p[N]; 
    for (n in 1:N) {
      p[n] <- inv_logit(mu[n,t]); 
      a[n] <- conc[t]*p[n];
      b[n] <- conc[t]*(1.0-p[n]);
      if (ns[t,n]>0) {
        ys[t,n] ~ beta_binomial(ns[t,n], a[n], b[n]);
      }
    }
    conc[t] ~ gamma(concShape, concRate);
  }
}