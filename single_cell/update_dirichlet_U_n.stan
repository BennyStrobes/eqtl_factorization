data {
  int<lower=0> T;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[T, K] V;   // predictor matrix
  matrix[T, K] V_squared;   // predictor matrix
  vector[T] y;      // outcome vector
  vector[T] g;
  vector[T] intercept;
  vector[T] F;
  vector[T] precisions; 
  real<lower=0> alpha_zero;
}
parameters {
  simplex[K] beta;
}
model {
  beta ~ dirichlet(rep_vector(alpha_zero, K)); 
  for (t in 1:T) {
    target += g[t]*precisions[t]*(y[t] - intercept[t] - (g[t]*F[t]))*(V[t,:]*beta);
    target += -precisions[t]*square(g[t])*(V_squared[t,:]*square(beta))/2.0;
    target += -precisions[t]*square(g[t])*((V[t,:]*beta)*(V[t,:]*beta))/2.0;
    for (k in 1:K) {
      target += precisions[t]*square(g[t])*square(V[t,k]*beta[k])/2.0;
    }
  }
}


