library(testthat)
library(maaslin3)

# When evaluate_only='abundance', maaslin_fit must rename pval/qval to
# pval_individual/qval_individual on fit_data_abundance$results.
# A previous bug referenced a nonexistent variable "outputs$results" instead
# of "fit_data_abundance$results", which would crash this code path.

test_that("evaluate_only='abundance' renames pval columns correctly", {
    data_in <- data.frame('a' = c(0, 0, 3, 0, 5, 7, 0, 8, 4, 7), 
                          'b' = c(2, 0, 0, 5, 0, 5, 2, 6, 9, 0))
    rownames(data_in) <- paste0("sample", 1:10)
    
    metadata <- data.frame('var1' = 1:10,
                           'var2' = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
    rownames(metadata) <- paste0("sample", 1:10)
    
    data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                      FUN = function(x) x / sum(x))))
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
                           median_comparison_abundance = FALSE,
                           warn_prevalence = FALSE,
                           evaluate_only = 'abundance')
    
    res_cols <- colnames(results$fit_data_abundance$results)
    expect_true("pval_individual" %in% res_cols)
    expect_true("qval_individual" %in% res_cols)
    expect_false("pval" %in% res_cols)
    expect_false("qval" %in% res_cols)
    expect_true("pval_joint" %in% res_cols)
    expect_true("qval_joint" %in% res_cols)
})
