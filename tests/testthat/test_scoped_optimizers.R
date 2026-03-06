library(testthat)
library(maaslin3)

# The optimizers and optCtrlList vectors are passed to mirai daemons via
# everywhere(). Verify they exist in the package namespace.

test_that("optimizers and optCtrlList are accessible via maaslin3 namespace", {
    expect_true(is.character(maaslin3:::optimizers))
    expect_equal(length(maaslin3:::optimizers), 4)
    expect_true(is.list(maaslin3:::optCtrlList))
    expect_equal(length(maaslin3:::optCtrlList), 4)
})

test_that("logistic random effects model completes without error", {
    set.seed(42)
    n <- 60
    data_in <- data.frame(
        a = rpois(n, lambda = 5),
        b = rpois(n, lambda = 3)
    )
    rownames(data_in) <- paste0("sample", 1:n)
    
    metadata <- data.frame(
        var1 = rnorm(n),
        group = factor(rep(paste0("g", 1:6), each = 10))
    )
    rownames(metadata) <- paste0("sample", 1:n)
    
    data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1,
                                      FUN = function(x) x / sum(x))))
    data_in_tss[data_in_tss == 0] <- NA
    data_in_tss_log <- log2(data_in_tss)
    
    results <- maaslin_fit(data_in_tss,
                           data_in_tss_log,
                           metadata,
                           formula = formula('expr ~ var1 + (1|group)'),
                           random_effects_formula = formula('expr ~ (1|group)'),
                           min_abundance = 0,
                           min_prevalence = 0,
                           min_variance = 0,
                           data = data_in,
                           median_comparison_abundance = FALSE,
                           median_comparison_prevalence = FALSE)
    
    expect_true(!is.null(results$fit_data_abundance))
    expect_true(!is.null(results$fit_data_prevalence))
    expect_true(nrow(results$fit_data_abundance$results) > 0)
    expect_true(nrow(results$fit_data_prevalence$results) > 0)
})
