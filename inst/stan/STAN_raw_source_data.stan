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
 corr_matrix[J] Sigma_k[K]; 

}



 model {
 // Prior on mu
  for(j in 1:J){
    for(k in 1:K){
      mu_jk[j,k] ~ normal(0,100);
    }
  }
  
for(k in 1:K){
  Sigma_k[k] ~ lkj_corr(shape_sig);
}



  // Likelihood
  for (i in 1:N) {
      Y[i] ~ multi_normal(mu_jk[,source[i]], Sigma_k[source[i]]);
    }
}
