library(testthat)
library(maaslin3)

output_tmp <- tempfile()
data_in <- data.frame('a' = c(0, 0, 0, 0, 0, 0, 0, 0, 7, 5, 8, 4, 7, 8, 5, 8), 
                      'b' = c(2, 0, 0, 5, 0, 5, 3, 6, 5, 2, 6, 9, 0, 8, 3, 7),
                      'c' = c(3, 4, 0, 0, 2, 7, 0, 3, 3, 7, 7, 2, 7, 4, 7, 0))
rownames(data_in) <- paste0("sample", seq(nrow(data_in)))

metadata <- data.frame('var1' = c(rep(0, 8), rep(1, 7), 2),
                       'var2' = c(0, 1, 0, 1, 0, 1, 0, 1, 
                                  0, 1, 0, 1, 0, 1, 0, 1))
rownames(metadata) <- paste0("sample", seq(nrow(data_in)))

data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0] <- NA
data_in_tss_log <- log2(data_in_tss)

results <- maaslin_fit(data_in_tss,
                       data_in_tss_log,
                       metadata,
                       formula = formula('expr ~ var1 + var2'),
                       random_effects_formula = NULL, 
                       min_abundance = 0, 
                       min_prevalence = 0, 
                       min_variance = 0,
                       data = data_in, 
                       median_comparison_abundance = FALSE)

plot_out <- maaslin_plot_results(output = output_tmp,
                     transformed_data = data_in_tss_log,
                     unstandardized_metadata = metadata,
                     fit_data_abundance = results$fit_data_abundance, 
                     fit_data_prevalence = results$fit_data_prevalence, 
                     normalization = 'TSS',
                     transform = 'LOG',
                     median_comparison_abundance = FALSE,
                     max_significance = 0.1)

expect_is(plot_out$assocation_plots$var1$a$logistic, 'ggplot')

expect_equal(list.files(file.path(output_tmp, 'figures', 
                                  'association_plots', 'var1', 'logistic')),
             'var1_a_logistic.png')

expect_equal(sort(list.files(file.path(output_tmp, 'figures'))),
             sort(c("association_plots", "summary_plot_gg.RDS", "summary_plot.pdf", "summary_plot.png")))

unlink(output_tmp, recursive = T)






