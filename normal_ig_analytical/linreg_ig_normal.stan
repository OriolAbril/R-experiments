data {
  int n;
  int m;
  vector[n] y;
  matrix[n, m] Z;
  matrix[m, m] V0;
  real g;
  real a;
  real b;
}

transformed data {
  vector[m] zeros_vect;

  for (j in 1:m) {
    zeros_vect[j] = 0;
  }
}

parameters {
  vector[m] theta;
  real<lower=0> phi;
}

model {
  // prior
  phi ~ inv_gamma(a/2, b/2);
  theta ~ multi_normal(zeros_vect, g*phi*V0);

  // likelihood
  y ~ normal(Z * theta, sqrt(phi));
}

