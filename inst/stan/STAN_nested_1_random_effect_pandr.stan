data {
  int<lower=1> J; // Number of isotopes
  int<lower=1> N; // Number of observations per group
  int<lower=1> K; // Dimension of p and q
  int<lower=1> L1; // N cols of X_fixed
  int<lower=1> L2; // N cols of X_random
  matrix[N, J] y; // Data matrix
  matrix[K, J] q; // q matrix
  matrix[K, J] s_mean; // s_mean matrix
  matrix[K, J] c_mean; // c_mean matrix
  matrix[K, J] s_sd; // s_sd matrix
  matrix[K, J] c_sd; // c_sd matrix
  vector[J] sigma_shape; // Prior shape for sigma
  vector[J] sigma_rate; // Prior rate for sigma
  real<lower=0> not_solo; // Adjustment factor for sigma
    matrix[N, L1] X_fixed; // Intercept
  matrix[N, L2] X_random; // Covariates matrix for pack/inside/nested one
  real<lower=0.001> cauchy_scale;
}

parameters {
  vector<lower=0>[J] sigma_raw; // log raw sigma values
  vector<lower=0>[K] omega; // log raw sigma values
    matrix[L1, K] beta0;
  matrix[L2, K] beta1; // Matrix of coefficients for inside values (also have intercept here I guess??)
}

transformed parameters {
  vector[J] sigma = 0.001 + not_solo * sigma_raw; // Transform sigma_raw
  matrix[N, K] p; // Main parameter
  matrix[N, J] var_y; // Variance for each group
  matrix[N, K] f; // f matrix for CLR prior on p

  for (i in 1:N) {
    for (k in 1:K) {
          f[i,k] = dot_product(X_fixed[i,:], beta0[:, k]) + dot_product(X_random[i,:], beta1[:,k]);

    }
  }

  for (i in 1:N) {
    p[i, :] = to_row_vector(softmax(to_vector(f[i, :])));
  }

  for (i in 1:N) {
    for (j in 1:J) {
      var_y[i,j] = dot_product(square(to_vector(p[i,:]) .* q[:, j]), square(s_sd[:, j]) + square(c_sd[:, j]))
                / square(dot_product(to_vector(p[i,:]), q[:, j])) + square(sigma[j]);
    }
  }

}

 model {
  // Prior on betas

       to_vector(beta0) ~ normal(0,1);

  for (k in 1:K){

  beta1[,k] ~ normal(0, omega[k]); // Prior for beta
  omega[k] ~ cauchy(0, cauchy_scale);
}

  // Prior on sigma_raw
  sigma_raw ~ gamma(sigma_shape, sigma_rate);

  // Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      real mu_ij = dot_product(to_vector(p[i, :]) .* q[:, j], s_mean[:, j] + c_mean[:, j]) / dot_product(to_vector(p[i, :]), q[:, j]);
      y[i, j] ~ normal(mu_ij, sqrt(var_y[i,j]));
    }
  }
}
