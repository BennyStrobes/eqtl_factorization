// Linear Model with Normal Errors
data {
  // number of observations
  int N;
  
  // response
  vector[N] y;
  
  // number of columns in the design matrix X
  int K;
  
  // design matrix X
  // should not include an intercept
  matrix [N, K] X;
  
  // Prior on beta
  real scale_beta;
}
parameters {
  // regression coefficient vector
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;

  mu = alpha + X * beta;
}
model {
  // priors
  beta ~ normal(0., scale_beta);
  // likelihood
  y ~ normal(mu, sigma);
}