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

  vector[N] sigma;
}
parameters {
  // regression coefficient vector
  vector<lower=0>[K] beta;
}
transformed parameters {
  vector[N] mu;

  mu = X * beta;
}
model {
  // priors
  beta ~ normal(0., scale_beta);
  // half normal
  // likelihood
  y ~ normal(mu, sigma);
}