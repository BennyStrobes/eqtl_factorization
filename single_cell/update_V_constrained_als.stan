data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  real<lower=0> scale_beta;
  real<lower=0> scale_beta_genotype;
  matrix[N, K] x;   // predictor matrix
  matrix[N, 1] x_genotype; // genotype predictor
  vector[N] y;      // outcome vector
}
parameters {
  real alpha;           // intercept
  vector<lower=-.4, upper=.4>[K] beta;       // coefficients for predictors
  vector[1] beta_genotype;
  real<lower=0> sigma;  // error scale
}
model {
  beta ~ normal(0., scale_beta);
  beta_genotype ~ normal(0., scale_beta_genotype);
  y ~ normal(x * beta + x_genotype * beta_genotype + alpha, sigma);  // likelihood
}

