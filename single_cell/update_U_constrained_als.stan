data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  simplex[K] beta;
  real<lower=0> sigma;  // error scale
}
model {
  y ~ normal(x * beta, sigma);  // likelihood
}

