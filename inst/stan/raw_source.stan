data {
  int<lower=1> J; // Number of isotopes
  int<lower=1> K; // Dimension of p and q (no. sources)
  int<lower=1> N; //N samples total in y
  matrix[N, J] Y; // Data matrix of food
  int<lower=1, upper=K> source[N];
  real<lower=1> shape_sig;
}

parameters {
 matrix[J,K] mu_jk; // Matrix of mean values
 corr_matrix[J] Omega[K];
 vector<lower=0>[K] tau; // Separate tau for each source
}

transformed parameters{
    matrix[J,J] Sigma[K];

     for (k in 1:K) {
    Sigma[k] = quad_form_diag(Omega[k], tau[k] * ones_vector(J)); // Covariance matrix for each source
  }

}

 model {
   tau ~ cauchy(0, 2.5);
 // Prior on mu
  for(j in 1:J){
    for(k in 1:K){
      mu_jk[j,k] ~ normal(0,100);
    }
  }

for(k in 1:K){
  Omega[k] ~ lkj_corr(shape_sig);
}



  // Likelihood
  for (i in 1:N) {
      Y[i] ~ multi_normal(mu_jk[,source[i]], Sigma[source[i]]);
    }
}
