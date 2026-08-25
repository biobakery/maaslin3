library(testthat)
library(maaslin3)

test_that("two crossed random intercepts return a stacked ranef matrix", {
    set.seed(42)
    n <- 60
    data_in <- data.frame(
        a = rpois(n, lambda = 8),
        b = rpois(n, lambda = 5)
    )
    rownames(data_in) <- paste0("sample", 1:n)

    metadata <- data.frame(
        x_cont = rnorm(n),
        cluster_a = factor(rep(paste0("a", 1:5), length.out = n)),
        cluster_b = factor(rep(paste0("b", 1:4), each = n / 4))
    )
    rownames(metadata) <- paste0("sample", 1:n)

    data_in_tss <- data.frame(t(apply(data_in, 1, function(x) x / sum(x))))
    data_in_tss[data_in_tss == 0] <- NA
    data_in_tss_log <- log2(data_in_tss)

    results <- maaslin_fit(
        data_in_tss,
        data_in_tss_log,
        metadata,
        formula = formula("expr ~ x_cont + (1 | cluster_a) + (1 | cluster_b)"),
        random_effects_formula = formula("expr ~ (1 | cluster_a) + (1 | cluster_b)"),
        min_abundance = 0,
        min_prevalence = 0,
        min_variance = 0,
        data = data_in,
        median_comparison_abundance = FALSE,
        median_comparison_prevalence = FALSE,
        bypass_small_group_warning = TRUE
    )

    ranef_ab <- results$fit_data_abundance$ranef
    expect_true(is.matrix(ranef_ab))
    expect_equal(nrow(ranef_ab), 2)
    expect_equal(rownames(ranef_ab), c("a", "b"))
    expect_equal(ncol(ranef_ab), 9)
    expect_true(all(grepl("^cluster_[ab]::", colnames(ranef_ab))))
    expect_true(is.numeric(ranef_ab))

    ranef_pr <- results$fit_data_prevalence$ranef
    expect_true(is.matrix(ranef_pr))
    expect_equal(nrow(ranef_pr), 2)
    expect_equal(colnames(ranef_pr), colnames(ranef_ab))
})
