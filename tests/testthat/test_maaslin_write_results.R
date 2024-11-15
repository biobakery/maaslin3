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

maaslin_write_results(output = output_tmp,
                      fit_data_abundance = results$fit_data_abundance, 
                      fit_data_prevalence = results$fit_data_prevalence, 
                      random_effects_formula = NULL, 
                      max_significance = 0.1)

results_combined <- rbind(results$fit_data_abundance$results,
                          results$fit_data_prevalence$results)
results_combined <- results_combined[order(results_combined$qval_individual, 
                                           na.last = T),]
results_combined$model <- ifelse(results_combined$model == 'linear',
                                 'abundance',
                                 'prevalence')
rownames(results_combined) <- NULL
all_results <- read.csv(file.path(output_tmp, 'all_results.tsv'), sep = '\t')
rownames(all_results) <- NULL
all_results$error <- as.character(all_results$error)
results_combined$error <- as.character(results_combined$error)

expect_equal(results_combined, all_results)

signif_results <- read.csv(file.path(output_tmp, 'significant_results.tsv'), sep = '\t')
rownames(signif_results) <- NULL
results_combined <- results_combined[results_combined$qval_joint < 0.1 | 
                                    results_combined$qval_individual < 0.1,]
rownames(results_combined) <- NULL

results_combined$error <- NULL
expect_equal(results_combined, signif_results)

unlink(output_tmp, recursive = T)










