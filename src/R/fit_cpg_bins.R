' fit_cpg_bins.R v1.0
Fit bins using Stan

Usage:
    fit_cpg_bins.R -i INPUT -o OUTPUT [ --method METHOD ]

Options:
    -i --input INPUT    Path to bin data output by bin_stats.R
    -o --output OUTPUT  Output path
    --method METHOD     Default is MCMC. Can specify LBFGS, BFGS, Newton, or VB.
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Fitting Bins vs. 1.0')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(GenomicRanges)
library(tidyverse)
library(flexmix)
library(countreg)
library(rstan)

DOWNSAMPLING_GROUP_SIZE = 200

bins <- read_tsv(args[['input']], comment = '#', col_types='ciiiddddi') %>%
  mutate(coverage_int = mean_coverage %>% round %>% as.integer) %>%
  filter(mean_coverage > 0 | gc_content > 0 | cpg_count > 0) %>%
  mutate(
    cpg_bin = factor(round(round(cpg_count * 20 / max(cpg_count)) * max(cpg_count) / 20)),
    gc_bin = factor(round(gc_content * 4, 1) / 4)
  )

##
## Determining the profile of non-specific binding at zero-CPG regions ##
##

zero_profile_gc_model_output <- bins %>%
  filter(cpg_count == 0) %>%
  mutate(
    gc_bin = factor(round(gc_content * 4, 1) / 4)
  ) %>%
  group_by(gc_bin) %>%
  filter(n() > 50) %>%
  ungroup() %>%
  plyr::ddply(c('gc_bin'), function(z) {
    print(sprintf('Running for gc_bin = %s', z$gc_bin %>% unique))
    if (all(z$coverage_int == 0)) {
      return(tibble(
        log_mu = -Inf,
        theta = Inf,
        log_likelihood = NA,
        count = nrow(z)
      ))
    } else {
      tryCatch({
        zero_model <- flexmix(coverage_int ~ 1, data = z, k = 1, model = FLXMRnegbin())
        return(tibble(
          log_mu = parameters(zero_model, component=1)[[1]],
          theta = parameters(zero_model, component=1)[[2]],
          log_likelihood = logLik(zero_model) %>% as.numeric,
          count = nrow(z)
        ))
      }, error = function(e) {
        message('Error: skipping this GC bin')
        NULL
      })
    }
  }) %>%
  mutate(mu = exp(log_mu))

##
## Fitting to logistic functions.
##

mu_fit <- nls(
  mu ~ a / (1 + exp(-b * (gc_bin - c))), start=list(a=1,b=10,c=0.6),
  data = zero_profile_gc_model_output %>% mutate(gc_bin = gc_bin %>% as.character %>% as.numeric),
  control = nls.control(warnOnly = TRUE, maxiter=1000)
)
mu_fit_coef = summary(mu_fit)$coef[, 1]

theta_fit <- nls(
  theta ~ a / (1 + exp(-b * (gc_bin - c))), start=list(a=1,b=10,c=0.6),
  data = zero_profile_gc_model_output %>% mutate(gc_bin = gc_bin %>% as.character %>% as.numeric) %>% filter(!is.infinite(theta)),
  control = nls.control(warnOnly = TRUE, maxiter=1000)
)
theta_fit_coef = summary(theta_fit)$coef[, 1]

bins_downsampled <- bind_rows(
    bins %>% group_by(cpg_bin, gc_bin) %>% filter(n() <= DOWNSAMPLING_GROUP_SIZE),
    bins %>% group_by(cpg_bin, gc_bin) %>% filter(n() > DOWNSAMPLING_GROUP_SIZE) %>% filter(row_number() %in% sample(1:n(), DOWNSAMPLING_GROUP_SIZE, replace=F))
)

stan_input <- with(
    bins_downsampled,
    list(
      N = nrow(bins_downsampled),
      coverage = coverage_int,
      ns_mu_a = mu_fit_coef[['a']],
      ns_mu_b = mu_fit_coef[['b']],
      ns_mu_c = mu_fit_coef[['c']],
      ns_theta_a = theta_fit_coef[['a']],
      ns_theta_b = theta_fit_coef[['b']],
      ns_theta_c = theta_fit_coef[['c']],
      cpg_count = cpg_count,
      gc_content = gc_content
    )
)

message('Compiling Stan model...')

cfmedip_stan_model <- stan_model(
  file = '/cluster/home/zhaoe/git/cfMeDIP-seq-analysis-pipeline/src/stan_model/fit_cfmedip.stan',
  verbose = TRUE
)

stop('testing')

if (is.null(args[['method']])) {
    method = 'MCMC'
} else {
    method = args[['method']]
}

message(sprintf('Running Stan model using method %s', method))

if (method %in% c('LBFGS', 'BFGS', 'Newton')) {
    stan_output <- optimizing(
        object = cfmedip_stan_model,
        data = stan_input,
    )
} else if (method == 'VB') {
    stan_output <- vb(
        object = cfmedip_stan_model,
        data = stan_input
    ) 
} else if (method == 'MCMC') {
    stan_output <- sampling(
        object = cfmedip_stan_model,
        data = stan_input,
        chains = 1,
        warmup = 100,
        iter = 300,
        cores = 1,
        verbose = TRUE
    )
} else {
    stop("Invalid method supplied. Must be MCMC, VB, LBFGS, BFGS, or Newton")
}

stan_output$input_data <- stan_input

saveRDS(stan_output, args[['output']])
