library(testthat)
library(maaslin3)
library(dplyr)

# The idea of these checks is that the contrast test should be the same as
# refactoring with categorical data
set.seed(1)

data_in <- data.frame('a' = c(0, 0, 3, 0, 5, 7, 0, 8, 4, 7), 
                      'b' = c(2, 0, 0, 5, 0, 5, 2, 6, 9, 0),
                      'c' = c(3, 4, 0, 0, 7, 2, 7, 4, 7, 0))
rownames(data_in) <- paste0("sample", c(1:10))

metadata <- data.frame('var1' = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       'var2' = c('a', 'a', 'a', 'a', 'b', 
                                  'b', 'b', 'c', 'c', 'c'))
metadata$var2 <- factor(metadata$var2)
rownames(metadata) <- paste0("sample", c(1:10))

data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0] <- NA
data_in_tss_log <- log2(data_in_tss)

out_dir = tempdir()

results <- maaslin_fit(data_in_tss,
                       data_in_tss_log,
                       metadata,
                       formula = formula('expr ~ var1 + var2'),
                       random_effects_formula = NULL, 
                       min_abundance = 0, 
                       min_prevalence = 0, 
                       min_variance = 0,
                       data = data_in, 
                       save_models = TRUE)

contrast_mat <- matrix(c(0, -1, 1), 
                       ncol = 3, nrow = 1, byrow = TRUE)

colnames(contrast_mat) <- c("var1",
                            "var2b",
                            "var2c")

contrast_test_out <- maaslin_contrast_test(results, 
                      contrast_mat)

metadata$var2 <- factor(metadata$var2, levels = c('b', 'a', 'c'))

results2 <- maaslin_fit(data_in_tss,
                       data_in_tss_log,
                       metadata,
                       formula = formula('expr ~ var1 + var2'),
                       random_effects_formula = NULL, 
                       min_abundance = 0, 
                       min_prevalence = 0, 
                       min_variance = 0,
                       data = data_in, 
                       save_models = TRUE)

new_mod_results <- results2$fit_data_abundance$results[
    results2$fit_data_abundance$results$value == 'c',]
new_mod_results <- new_mod_results[order(new_mod_results$pval_individual),]

expect_that(new_mod_results$coef, equals(contrast_test_out$fit_data_abundance$results$coef))
expect_that(new_mod_results$stderr, equals(contrast_test_out$fit_data_abundance$results$stderr))
expect_equal(new_mod_results$pval_individual, 
             contrast_test_out$fit_data_abundance$results$pval_individual, tolerance = 0.01)

new_mod_results2 <- results2$fit_data_prevalence$results[
    results2$fit_data_prevalence$results$value == 'c',]
new_mod_results2 <- new_mod_results2[order(new_mod_results2$pval_individual),]

contrast_test_out$fit_data_prevalence$results <- 
    contrast_test_out$fit_data_prevalence$results[
        order(contrast_test_out$fit_data_prevalence$results$pval_individual),]

expect_that(new_mod_results2$coef, equals(contrast_test_out$fit_data_prevalence$results$coef))
expect_that(new_mod_results2$stderr, equals(contrast_test_out$fit_data_prevalence$results$stderr))
expect_equal(new_mod_results2$pval_individual, 
             contrast_test_out$fit_data_prevalence$results$pval_individual, tolerance = 0.01)

unlink(out_dir, recursive = TRUE)

