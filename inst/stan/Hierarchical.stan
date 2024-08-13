data {
    int<lower=1> J;  // number of tracers
    int<lower=1> K; //number of sources
  int<lower=1> n[K];           // sample size for each source

  matrix[K, J] source_mean; // mean for each source 1 to K for isotope 1 to J
  matrix[K,J] source_sd; //sd for each source 1 to K for isotope 1 to J 1 to J
}

parameters {
   matrix[K, J] tmp_X; // Temporary matrix for chi-square random variables

  matrix[K,J] mu; //true mean to be estimated
}

transformed parameters {
   matrix[K, J] src_tau; // The final tau values

  // Generate src_tau values
  for (k in 1:K) {
    for (j in 1:J) {
      src_tau[k, j] = tmp_X[k, j] / (source_sd[k, j] * (n[k] - 1));
    }
  }
}

model {

for(k in 1:K){
for(j in 1:J){
mu[k, j] ~ normal(source_mean[k,j], n[k]/source_sd[k,j]);
tmp_X[k,j] ~ chi_square(n[k]);
}
}



}

generated quantities {
  matrix[K,J] sigma;
  sigma = 1/sqrt(src_tau);
}
