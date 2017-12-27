data {
  int<lower=1> N;  // sample number
  int<lower=1> K;  // category number
  int<lower=1> D;  // covariate number
  row_vector[D] cov[N];  // covariate data
  vector[N] geno;  // genotype data
  int<lower=0,upper=K> pheno[N];  // phenotype data
  int<lower=0,upper=1> prior_non; // non-informative prior
  real prior_inf; // informative prior
}

parameters {
  real a;  // intercept
  real p;  // variant effect
  vector[D] beta;  // covariate effect
  simplex[K-1] cut1;  // cut points primitive
  real xx;
}

transformed parameters {
  ordered[K-1] cut; // cut points
  real cutsum; // cut points

  cut[1] = 0; // solve identity issue with the intercept
  cutsum = 0; // make the ever-growing vector
  for (cutp in 2:(K-1)) {
    cutsum = cutsum + cut1[cutp-1];
    cut[cutp] = 10 * cutsum;
  }
}

model {
  xx ~ normal(0, 1e-4);
  if (prior_non == 1) { // otherwise use flat prior for all parameters
    a ~ normal(0, 1); // intercept
    beta ~ normal(0, 1); // covariates

    if (prior_inf == 0) { // non-informative prior for variant effect
      p ~ normal(0, 1); // variant effect
    }
  } // end of if flow for specifying priors

  if (K > 2) { // ordered categorical model
    cut1 ~ dirichlet(rep_vector(1, K-1));

    for (n in 1:N)
      pheno[n] ~ ordered_logistic(a + p * geno[n] + cov[n] * beta, cut);
  }
}
