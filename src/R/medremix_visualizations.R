plot_scatter_bin_methylation <- function(bin_methylation_data, nbfit) {
  scatterplot_data <- bin_methylation_data %>%
    group_by(gc_bin) %>%
    summarise(max_cpg_count = max(cpg_count)) %>%
    ungroup() %>%
    filter(max_cpg_count > 0) %>%
    expand_grid(tibble(cpg_count = 1:50)) %>%
    filter(cpg_count <= max_cpg_count) %>%
    mutate(
      gc_content = gc_bin %>% as.character %>% as.numeric
    ) %>%
    mutate(
      predicted_coverage = exp(predict(nbfit, newdata=.)),
      gc_label = sprintf('%s%% GC', gc_bin %>% as.character() %>% as.numeric * 100)
    )

  bin_methylation_data %>%
    mutate(
      methylation_status = factor(methylation_status, levels = c('methylated', 'unmethylated')),
      gc_label = sprintf('%s%% GC', gc_bin %>% as.character() %>% as.numeric * 100)
    ) %>%
    ggplot(aes(
      x = cpg_count,
      y = coverage,
    )) +
    geom_point(aes(
      color = methylation_status
      ), alpha=0.1, shape=20
    ) +
    geom_line(aes(
      x = cpg_count,
      y = predicted_coverage
    ), data = scatterplot_data) +
    facet_wrap(~ gc_label, scales='free') +
    labs(
      x = 'CpG count',
      y = 'Coverage'
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)
    ) +
    scale_color_brewer(palette = 'Set1')
}

plot_bin_methylation_fit <- function(bin_methylation_data, nbfit) {
  bin_methylation_fit_plot_data_raw <- bin_methylation_data %>%
    group_by(cpg_bin, gc_bin) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    expand_grid(tibble(x = 1:50)) %>%
    mutate(
      gc_content = gc_bin %>% as.character %>% as.numeric,
      cpg_count = cpg_bin %>% as.character %>% as.integer
    )
  
  bin_methylation_fit_plot_data <- bin_methylation_fit_plot_data_raw %>%
    mutate(
      unmethylated_mu = exp(predict(nbfit$zero_model$mu_fit, newdata=.)),
      unmethylated_theta = exp(predict(nbfit$zero_model$theta_fit, newdata=.)),
      unmethylated_fit = dnbinom(x, mu = unmethylated_mu, size = unmethylated_theta) * count,
      methylated_mu = ifelse(cpg_count == 0, NA, exp(predict(nbfit$final_model, newdata = . ))),
      methylated_theta = nbfit$final_model$theta,
      methylated_fit = dnbinom(x, mu = methylated_mu, size =  methylated_theta) * count
    ) %>%
    filter(!is.nan(methylated_fit))
  
  bin_methylation_fit_plot_data_reshaped <- bin_methylation_fit_plot_data %>%
    dplyr::select(cpg_bin, gc_bin, x, methylated_fit, unmethylated_fit) %>%
    gather(fit, value, -cpg_bin, -gc_bin, -x) %>%
    arrange(cpg_bin, gc_bin) %>%
    mutate(facet_label = sprintf(
      '%s CpGs, GC: %s%%',
      cpg_bin %>% as.character,
      gc_bin %>% as.character %>% as.numeric * 100 %>% round
    )) %>%
    mutate(facet_label = factor(facet_label, levels = unique(facet_label)))
  
  return(
    bin_methylation_data %>%
      arrange(cpg_bin, gc_bin) %>%
      mutate(facet_label = sprintf(
        '%s CpGs, GC: %s%%',
        cpg_bin %>% as.character,
        gc_bin %>% as.character %>% as.numeric * 100 %>% round
      )) %>%
      mutate(facet_label = factor(facet_label, levels = unique(facet_label))) %>%
      ggplot(aes(
        x = coverage,
      )) +
      geom_histogram(aes(
        fill = methylation_status
      ), alpha=0.3) +
      geom_line(aes(
        x = x,
        y = value,
        group = fit,
        color = fit
      ), size = 1, data = bin_methylation_fit_plot_data_reshaped) +
      facet_wrap(~ facet_label, scales = 'free') +
      labs(
        x = 'Coverage',
        y = 'Number of bins',
        color = 'Model fit',
        fill = 'Classification'
      ) +
      scale_color_brewer(palette = 'Set1') +
      scale_fill_brewer(palette = 'Set1')
  )
}
