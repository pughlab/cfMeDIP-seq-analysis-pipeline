' cfmedip_nbglm.R
Fit bins from cfMeDIP-seq data using negative binomial regression.

Usage:
    cfmedip_nbglm.R -i INPUT -o OUTPUT [ --method METHOD --modelout MODELOUT ]

Options:
    -i --input INPUT            Path to bin data output by bin_stats.R
    -o --output OUTPUT          Output path for the table of

    --method METHOD             Default is MCMC. Can specify LBFGS, BFGS, Newton, or VB.
    --modelout MODELOUT         Path into which to dump a .Rds file containing the final
                                    model specifications and the coefficients used in
                                    each iteration.
' -> doc

if (! interactive()) {
  library(docopt)
  args <- docopt(doc, version='cfMeDIP-seq negative binomial GLM v1.0')
  print(args)
} else {
  message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)
library(MASS)
library(flexmix)
library(countreg)

MAX_ITER = 20
PERCENT_CHANGE_THRESHOLD = 0.1

message('- Importing data.')

bins <- read_tsv(args[['input']], comment = '#', col_types='ciiddddddi') %>%
  mutate(coverage_int = mean_coverage %>% round %>% as.integer) %>%
  filter(mean_coverage > 0 | gc_content > 0 | cpg_count > 0) %>%
  mutate(
    cpg_bin = factor(round(round(cpg_count * 20 / max(cpg_count)) * max(cpg_count) / 20)),
    gc_bin = factor(round(gc_content * 4, 1) / 4)
  )

##
## Determining the profile of non-specific binding at zero-CPG regions ##
##

message('- Determining the profile of bins with no CpGs')

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
## The zero profile to logistic functions, depending on gc_content.
##

zero_mu_fit <- lm(
  log(mu) ~ log(gc_content),
  data = zero_profile_gc_model_output %>% mutate(gc_content = gc_bin %>% as.character %>% as.numeric) %>% filter(mu > 0, gc_content > 0)
)

zero_theta_fit <- lm(
  log(theta) ~ log(gc_content),
  data = zero_profile_gc_model_output %>% mutate(gc_content = gc_bin %>% as.character %>% as.numeric) %>% filter(!is.infinite(theta), theta < 100, gc_content > 0)
)

##
## Begin fitting of CpG by making an initial estimate based on
##      mean_coverage = exp(b_0 + b_1 * cpg_count + b_2 * gc_content)
## setting b_1 to 0.5, b_2 to 0, and allowing theta to vary as 0.5 * cpg_count.
##

message('- Fitting CpG profile iteratively.\n')

bin_methylation <- bins %>%
  mutate(
    unmethylated_mu = exp(predict(zero_mu_fit, newdata = .)),
    unmethylated_theta = exp(predict(zero_theta_fit, newdata = .)),
    unmethylated_likelihood = dnbinom(coverage_int, mu = unmethylated_mu, size = unmethylated_theta),
    methylated_mu = ifelse(cpg_count == 0, NA, cpg_count),
    methylated_theta = ifelse(cpg_count == 0, NA, 0.5 * cpg_count),
    methylated_likelihood = ifelse(is.na(methylated_mu), 0, dnbinom(coverage_int, mu = methylated_mu, size = methylated_theta)),
    unmethylated_posterior = unmethylated_likelihood / (methylated_likelihood + unmethylated_likelihood),
    methylated_posterior = methylated_likelihood / (methylated_likelihood + unmethylated_likelihood)
  ) %>%
  mutate(methylation_status = ifelse(methylated_posterior > unmethylated_posterior, 'methylated', 'unmethylated'))

##
## Iterate: For each iteration...
##      1. Separate out the presumed methylated bins as determined by the current model.
##      2. Fit a new negative binomial GLM regression model
##              coverage = exp(b_0 + b_1 * cpg_count + b_2 * gc_content + b_3 * cpg_count * gc_content)
##      3. Re-estimate the posterior probabilities of each point being methylated or unmethylated
##         based on the new negative binomial regression model.
##      4. Check for convergence, defined as < 1% change in all coefficients from the previous iteration.
##         Stop if convergence or the maximum number of allowable iterations has been reached.
##

methylated_fits <- list()
for (i in 1:MAX_ITER) {
  ##
  ## 1. Separate out the presumed methylated bins
  ##

  message(sprintf('Running iteration %s', i))
  bin_methylation_subset <- bin_methylation %>%
    filter(methylated_posterior > 0.5)
  
  ##
  ## 2. Fit NB regression to the methylated bins
  ##

  methylated_bins_nbfit <- glm.nb(
    coverage_int ~ cpg_count * gc_content,
    data = bin_methylation_subset
  )
  methylated_fits[[i]] = methylated_bins_nbfit

  ##
  ## 3. Refit the points based on the regression
  ##
  
  refitted <- bins %>%
    mutate(
      unmethylated_mu = exp(predict(zero_mu_fit, newdata = .)),
      unmethylated_theta = exp(predict(zero_theta_fit, newdata = .)),
      unmethylated_likelihood = dnbinom(coverage_int, mu = unmethylated_mu, size = unmethylated_theta),
      methylated_mu = ifelse(cpg_count == 0, NA, exp(predict(methylated_bins_nbfit, newdata = .))),
      methylated_theta = methylated_bins_nbfit$theta,
      methylated_likelihood = ifelse(is.na(methylated_mu), 0, dnbinom(coverage_int, mu = methylated_mu, size = methylated_theta)),
      unmethylated_posterior = unmethylated_likelihood / (methylated_likelihood + unmethylated_likelihood),
      methylated_posterior = methylated_likelihood / (methylated_likelihood + unmethylated_likelihood)
    )
  
  if (any(is.nan(refitted$unmethylated_posterior)) || any(is.nan(refitted$methylated_posterior))) {
    nan_rows = refitted %>% filter(is.nan(unmethylated_posterior) | is.nan(methylated_posterior))
    message(sprintf('%s bins yielded NaNs - these may be outliers and have been removed. They are printed below.', nrow(nan_rows)))
    print(nan_rows)
    
    refitted <- refitted %>% filter(!is.nan(methylated_posterior), !is.nan(unmethylated_posterior))
  }
  
  bin_methylation <- refitted %>%
    mutate(
      methylation_status = ifelse(methylated_posterior > unmethylated_posterior, 'methylated', 'unmethylated')
    )

  ##
  ## 4. Check for convergence
  ##
  
  if (i > 1) {
    percent_changes = abs((methylated_fits[[i]]$coefficients - methylated_fits[[i-1]]$coefficients) / methylated_fits[[i-1]]$coefficients) * 100
    message('Percent Changes:')
    message(sprintf('%s  %s  %s\n', names(percent_changes), signif(percent_changes, 3), ifelse(percent_changes < PERCENT_CHANGE_THRESHOLD, 'converged', '--')))
    if (all(percent_changes < PERCENT_CHANGE_THRESHOLD)) {
      message('All coefficients converged.')
      methylated_fit <- methylated_bins_nbfit
      break
    } else if (i == MAX_ITER) {
      message('Maximum iterations hit. Some coefficients did not converge.')
    }
  }
}

##
## Write output
##

message('- Writing output to file.')

bin_methylation %>%
  dplyr::select(
    bin_chr,
    bin_start,
    bin_end,
    methylation_status,
    coverage = coverage_int,
    cpg_count,
    gc_content,
    methylated_posterior,
    methylated_mu,
    mean_fragment_length
  ) %>%
  write_tsv(args[['output']])
message(sprintf('Output data written to %s', args[['output']]))

if (!is.null(args[['modelout']])) {
  message('- Serializing model data to file.')
  model_data <- list(
    final_model = methylated_fit,
    iteration_models = methylated_fits,
    zero_model = list(
        model_output = zero_profile_gc_model_output,
        theta_fit = zero_theta_fit,
        mu_fit = zero_mu_fit
    )
  )
  saveRDS(model_data, args[['modelout']])
  message(sprintf('Model data written to %s', args[['modelout']]))
}
