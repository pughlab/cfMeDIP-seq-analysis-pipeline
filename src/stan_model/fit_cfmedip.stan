data {
  int<lower=1> N;          // number of data points
  int<lower=0> coverage[N];               // observations
  
  real<lower=0> ns_mu_a;				// Nonspecific binding logistic coefficients
  real<lower=0> ns_mu_b;
  real<lower=0> ns_mu_c;
  real<lower=0> ns_theta_a;
  real<lower=0> ns_theta_b;
  real<lower=0> ns_theta_c;
  
  vector<lower=0>[N] cpg_count;
  vector<lower=0, upper=1>[N] gc_content;
}

parameters {
  real<lower=0, upper=1> theta;		// proportion of methylated component
  real<lower=0> cpg_max;			// Logistic coefficients for regression of mu
  real<lower=0> cpg_growth;
  real<lower=0> cpg_mean;
  real gc_coef;
  real intercept;
  real methylated_mu_sigma[N];
  vector[N] methylated_mu;
  vector[N] methylated_theta;
}

transformed parameters {
  vector[N] methylated_mu_mean;
  vector[N] unmethylated_mu;
  vector[N] unmethylated_theta;
  
  methylated_mu_mean = intercept + (cpg_max ./ (1 + exp(-1 * cpg_growth * ( cpg_count - cpg_mean )))) + (gc_coef * gc_content);
  unmethylated_mu = ns_mu_a ./ (1 + exp(-1 * ns_mu_b * (gc_content - ns_mu_c)));
  unmethylated_theta = ns_theta_a ./ (1 + exp(- ns_theta_b * (gc_content - ns_theta_c)));
}

model {  
  cpg_growth ~ gamma(1, 1/10);
  cpg_mean ~ gamma(3, 1/10);
  cpg_max ~ gamma(5, 1/20);
  gc_coef ~ normal(0, 5);
  intercept ~ normal(0, 5);
  methylated_mu_sigma ~ lognormal(0, 2);
  methylated_mu ~ normal(methylated_mu_mean, methylated_mu_sigma);
  methylated_theta ~ lognormal(0, 2);
    
  for (n in 1:N) {
    real unmethylated = log(1-theta) + neg_binomial_lpmf(coverage[n] | unmethylated_mu[n], unmethylated_theta[n]);
	real methylated = log(theta) + neg_binomial_lpmf(coverage[n] | methylated_mu[n], methylated_theta[n]);
    target += log_sum_exp(methylated, unmethylated);
  }
}
