data {
  int<lower=1> J; // Number of isotopes
  int<lower=1> N; // Number of observations per group
  int<lower=1> K; // Dimension of p and q
  int<lower=1> L1; // Number of covariates in inner
  int<lower=1> L2; //Number of covariates in outer
  matrix[N, J] y; // Data matrix
  matrix[K, J] q; // q matrix
  matrix[K, J] s_mean; // s_mean matrix
  matrix[K, J] c_mean; // c_mean matrix
  matrix[K, J] s_sd; // s_sd matrix
  matrix[K, J] c_sd; // c_sd matrix
  vector[J] sigma_shape; // Prior shape for sigma
  vector[J] sigma_rate; // Prior rate for sigma
  real<lower=0> not_solo; // Adjustment factor for sigma
  matrix[N, 1] X_intercept; // Covariates matrix for pack/inside/nested one
  matrix[N, L1] X_inner; // Covariates matrix for pack/inside/nested one
  matrix[N, L2] X_outer; // Covariates matrix for pack/inside/nested one
 int<lower=1> region_for_pack[L1];
  vector[J] omicron_shape; // Prior shape for sigma
  vector[J] omicron_rate; // Prior rate for sigma
}

parameters {
  vector<lower=0>[J] sigma_raw; // log raw sigma values
 matrix<lower=0>[K, L1] sigma_pack; // log raw sigma values
 vector<lower=0>[K] sigma_region; // log raw sigma values
  vector<lower=0>[K] mu_region; // mean for beta1
  matrix[K, L1] beta1; // Matrix of coefficients for inside values (also have intercept here I guess??)
   matrix[K, L2] beta2; // Matrix of coefficients for outside values

}

transformed parameters {
  vector[J] sigma = 0.001 + not_solo * sigma_raw; // Transform sigma_raw
 // vector[N] alpha;
 vector[J] omicron;
  matrix[N, K] p; // Main parameter
  matrix[N, J] var_y; // Variance for each group
  matrix[N, K] f; // f matrix for CLR prior on p

  for (i in 1:N) {
    for (k in 1:K) {
      f[i,k] = dot_product(X_inner[i,:], beta1[k,:]) + dot_product(X_outer[i,:], beta2[k,:]);
         // f[i,k] = dot_product(X_inner[i,:], beta1[k,:]);

    }
  }

  for (i in 1:N) {
    p[i, :] = to_row_vector(softmax(to_vector(f[i, :])));
  }

  for (i in 1:N) {
    for (j in 1:J) {
      var_y[i,j] = dot_product(square(to_vector(p[i,:]) .* q[:, j]), square(s_sd[:, j]) + square(c_sd[:, j]))
                / square(dot_product(to_vector(p[i,:]), q[:, j])) .*square(omicron[j]) + square(sigma[j]);
    }
  }

}

 model {
  // Prior on beta region
  for (k in 1:K) {
    for (l in 1:L2) {

      beta2[k,l] ~ normal(mu_region[k], square(sigma_region[k])); // Prior for beta

    }
    mu_region[k] ~ normal(0,1);

  }


  // Prior on beta pack
  // Need to match pack number to region number

  for(k in 1:K){
    for(l in 1:L1){
      sigma_pack[k,l] ~ gamma(1,1);
      beta1[k,l] ~ normal(beta2[k,region_for_pack[l]], square(sigma_pack[k,l]));
    }
  }



   //To match mixsiar




  // Prior on sigma_raw
  sigma_region ~ gamma(1,1);

  sigma_raw ~ gamma(1, 1);
  omicron ~ gamma(omicron_shape, omicron_rate);
//  alpha ~ normal(0,1);

  // Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      real mu_ij = dot_product(to_vector(p[i, :]) .* q[:, j], s_mean[:, j] + c_mean[:, j]) / dot_product(to_vector(p[i, :]), q[:, j]);
      y[i, j] ~ normal(mu_ij, sqrt(var_y[i,j]));
    }
  }
}
