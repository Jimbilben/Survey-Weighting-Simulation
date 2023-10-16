// generated with brms 2.17.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] weights;  // model weights
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(sigma | 3, 0, 15.6)
    - 1 * student_t_lccdf(0 | 3, 0, 15.6);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = X * b;
    for (n in 1:N) {
      target += weights[n] * (normal_lpdf(Y[n] | mu[n], sigma));
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
}
