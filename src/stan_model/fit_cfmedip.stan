data {
  int<lower=1> N;               // Number of data points
  int<lower=0> coverage[N];     // Coverage at each bin
   
  vector<lower=0>[N] cpg_count; // Number of CpGs in bin
  vector<lower=0, upper=1>[N] gc_content; // GC content in bin
 
  real<lower=0> ns_mu_a;        // Nonspecific binding logistic coefficients
  real<lower=0> ns_mu_b;
  real<lower=0> ns_mu_c;
  real<lower=0> ns_theta_a;
  real<lower=0> ns_theta_b;
  real<lower=0> ns_theta_c;

  int<lower=0, upper=1> run_estimation; // a switch to evaluate likelihood. Use 0 for simulating data.
}

parameters {
  real<lower=0, upper=1> theta; // proportion of methylated component
  //real<lower=0> cpg_max;        // Logistic coefficients for regression of mu
  //real<lower=0> cpg_growth;
  //real<lower=0> cpg_mean;
  real<lower=0> cpg_slope;
  //real<lower=0> gc_coef;
  real<lower=0> gc_slope;
  real<lower=0> intercept;
  real methylated_mu_sigma;
  real methylated_phi_intercept;
  real methylated_phi_cpg;
  real methylated_phi_gc;
}

model {  
  vector[N] methylated_mu;
  vector[N] methylated_phi;
  vector[N] unmethylated_mu;
  vector[N] unmethylated_phi;

  // Priors
  theta ~ beta(1, 1);
  //cpg_max ~ gamma(5, 1.0/20.0);
  //cpg_growth ~ gamma(1, 1.0/10.0);
  //cpg_mean ~ gamma(3, 1.0/10.0);
  //gc_coef ~ gamma(1, 1.0/10.0);
  cpg_slope ~ gamma(1, 1);
  gc_slope ~ gamma(1, 1);
  intercept ~ gamma(1, 1);
  methylated_phi_intercept ~ normal(0, 1);
  methylated_phi_cpg ~ normal(0, 1);
  methylated_phi_gc ~ normal(0, 1);

  unmethylated_mu = ns_mu_a ./ (1 + exp(-1 * ns_mu_b * (gc_content - ns_mu_c)));
  unmethylated_phi = ns_theta_a ./ (1 + exp(- ns_theta_b * (gc_content - ns_theta_c)));

  // methylated_mu = intercept + (cpg_max ./ (1 + exp(-1 * cpg_growth * ( cpg_count - cpg_mean )))) + (gc_coef * gc_content);
  methylated_mu = intercept + cpg_slope * cpg_count + gc_slope * gc_content;
  methylated_phi = exp(methylated_phi_intercept + methylated_phi_cpg * cpg_count + methylated_phi_gc * gc_content);

  if (run_estimation == 1) {
      for (n in 1:N) {
        target += log_mix(
          theta,
          neg_binomial_2_lpmf(coverage[n] | methylated_mu[n], methylated_phi[n]),
          neg_binomial_2_lpmf(coverage[n] | unmethylated_mu[n], unmethylated_phi[n])
        );
      }
  }
}

generated quantities {
  vector[N] methylated_mu;
  vector[N] methylated_phi;
  vector[N] unmethylated_mu;
  vector[N] unmethylated_phi;
  vector[N] coverage_sim;
  int<lower=0, upper=1> is_methylated;

  unmethylated_mu = ns_mu_a ./ (1 + exp(-1 * ns_mu_b * (gc_content - ns_mu_c)));
  unmethylated_phi = ns_theta_a ./ (1 + exp(- ns_theta_b * (gc_content - ns_theta_c)));
  methylated_mu = intercept + (cpg_max ./ (1 + exp(-1 * cpg_growth * ( cpg_count - cpg_mean )))) + (gc_coef * gc_content);
  methylated_phi = exp(methylated_phi_intercept + methylated_phi_cpg * cpg_count + methylated_phi_gc * gc_content);

  for (n in 1:N) {
    is_methylated = bernoulli_rng(theta);
    if (is_methylated == 1) {
        coverage_sim[n] = neg_binomial_2_rng(methylated_mu[n], methylated_phi[n]);
    } else {
        coverage_sim[n] = neg_binomial_2_rng(unmethylated_mu[n], unmethylated_phi[n]);
    }
  }
}
