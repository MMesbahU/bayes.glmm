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
  real<lower=1e-4> sigma_gau; // variance of gaussian models
  simplex[K-1] cut1;  // cut points primitive
  real t; // standardized effect size
  real<lower=1e-4> sigma_inf; // prior variance
}

transformed parameters {
    ordered[K-1] cut; // cut points
    real cutsum; // cut points
    real mu; // prior effect size

    cut[1] = 0; // solve identity issue with the intercept
    cutsum = 0; // make the ever-growing vector
    for (cutp in 2:(K-1)) {
      cutsum = cutsum + cut1[cutp-1];
      cut[cutp] = 10 * cutsum;
    }

    mu = sigma_inf * t; // informative prior
}

model {
  t ~ normal(prior_inf, 1); // standarized size of variant effect
  sigma_inf ~ inv_gamma(2, 1); // variance of variant effect

  sigma_gau ~ inv_gamma(2, 1);
  cut1 ~ dirichlet(rep_vector(1, K-1));

  if (prior_non == 1) { // otherwise use flat prior for all parameters
    a ~ normal(0, 1); // intercept
    beta ~ normal(0, 1); // covariates

    if (prior_inf == 0) { // non-informative prior for variant effect
      p ~ normal(0, 1); // variant effect
    }
    if (prior_inf != 0) { // informative prior for variant effect
      p ~ normal(mu, sigma_inf); // variant effect
    }
  } // end of if flow for specifying priors

  if (K == 1) { // Gaussian model
    for (n in 1:N)
      pheno[n] ~ normal(a + p * geno[n] + cov[n] * beta, sigma_gau);
  }

  if (K == 2) { // binary logistic model
    for (n in 1:N)
      pheno[n] ~ bernoulli_logit(a + p * geno[n] + cov[n] * beta);
  }

  if (K > 2) { // ordered categorical model
    for (n in 1:N)
      pheno[n] ~ ordered_logistic(a + p * geno[n] + cov[n] * beta, cut);
  }
}
