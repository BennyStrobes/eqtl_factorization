data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
  real tau;  // temperature
}
parameters {
  simplex[K] pi;
  simplex[K] beta;
  real<lower=0> sigma;  // error scale
}
model {
  real temp = 0.0;
  for (k in 1:K) {
    target += log(pi[k]) - ((tau+1)*log(beta[k]));
    temp += (pi[k])/(pow(beta[k], tau));
   }
  target += (1.0/K)*log(temp);
  y ~ normal(x * beta, sigma);  // likelihood
}

