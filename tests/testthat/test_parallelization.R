library(testthat)
library(maaslin3)

setup_test_data <- function() {
    set.seed(123)
    n <- 30
    data_in <- data.frame(
        a = rpois(n, lambda = 5),
        b = rpois(n, lambda = 3),
        c = rpois(n, lambda = 4)
    )
    rownames(data_in) <- paste0("sample", 1:n)

    metadata <- data.frame(
        var1 = rnorm(n),
        var2 = c(rep(0, 15), rep(1, 15))
    )
    rownames(metadata) <- paste0("sample", 1:n)

    data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1,
                                      FUN = function(x) x / sum(x))))
    data_in_tss[data_in_tss == 0] <- NA
    data_in_tss_log <- log2(data_in_tss)

    list(data_in = data_in, metadata = metadata,
         data_in_tss = data_in_tss, data_in_tss_log = data_in_tss_log)
}

test_that("serial execution (cores=1, no daemons) produces results", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_test_data()
    results <- maaslin_fit(td$data_in_tss,
                           td$data_in_tss_log,
                           td$metadata,
                           formula = formula('expr ~ var1 + var2'),
                           random_effects_formula = NULL,
                           min_abundance = 0,
                           min_prevalence = 0,
                           min_variance = 0,
                           data = td$data_in,
                           cores = 1,
                           median_comparison_abundance = FALSE,
                           median_comparison_prevalence = FALSE)

    expect_false(mirai::daemons_set())
    expect_true(!is.null(results$fit_data_abundance))
    expect_true(!is.null(results$fit_data_prevalence))
    expect_true(nrow(results$fit_data_abundance$results) > 0)
    expect_true(nrow(results$fit_data_prevalence$results) > 0)
})

test_that("cores parameter starts daemons and cleans up", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_test_data()
    results <- maaslin_fit(td$data_in_tss,
                           td$data_in_tss_log,
                           td$metadata,
                           formula = formula('expr ~ var1 + var2'),
                           random_effects_formula = NULL,
                           min_abundance = 0,
                           min_prevalence = 0,
                           min_variance = 0,
                           data = td$data_in,
                           cores = 2,
                           median_comparison_abundance = FALSE,
                           median_comparison_prevalence = FALSE)

    expect_false(mirai::daemons_set())
    expect_true(!is.null(results$fit_data_abundance))
    expect_true(!is.null(results$fit_data_prevalence))
    expect_true(nrow(results$fit_data_abundance$results) > 0)
    expect_true(nrow(results$fit_data_prevalence$results) > 0)
})

test_that("pre-set daemons are used and preserved", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)
    mirai::daemons(2)
    Sys.sleep(0.5)

    expect_true(mirai::daemons_set())

    td <- setup_test_data()
    results <- maaslin_fit(td$data_in_tss,
                           td$data_in_tss_log,
                           td$metadata,
                           formula = formula('expr ~ var1 + var2'),
                           random_effects_formula = NULL,
                           min_abundance = 0,
                           min_prevalence = 0,
                           min_variance = 0,
                           data = td$data_in,
                           cores = 1,
                           median_comparison_abundance = FALSE,
                           median_comparison_prevalence = FALSE)

    expect_true(mirai::daemons_set())
    expect_true(!is.null(results$fit_data_abundance))
    expect_true(!is.null(results$fit_data_prevalence))
    expect_true(nrow(results$fit_data_abundance$results) > 0)
    expect_true(nrow(results$fit_data_prevalence$results) > 0)

    mirai::daemons(0)
    Sys.sleep(0.5)
})

test_that("cores>1 warns when daemons already set", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)
    mirai::daemons(2)
    Sys.sleep(0.5)

    td <- setup_test_data()
    expect_warning(
        maaslin_fit(td$data_in_tss,
                    td$data_in_tss_log,
                    td$metadata,
                    formula = formula('expr ~ var1 + var2'),
                    random_effects_formula = NULL,
                    min_abundance = 0,
                    min_prevalence = 0,
                    min_variance = 0,
                    data = td$data_in,
                    cores = 2,
                    median_comparison_abundance = FALSE,
                    median_comparison_prevalence = FALSE),
        "Daemons already set"
    )

    expect_true(mirai::daemons_set())
    mirai::daemons(0)
    Sys.sleep(0.5)
})

test_that("maaslin3 function with cores parameter works end-to-end", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_test_data()
    tmp_output <- tempfile("maaslin3_parallel_test_")
    dir.create(tmp_output)

    results <- maaslin3(td$data_in,
                        td$metadata,
                        output = tmp_output,
                        formula = 'expr ~ var1 + var2',
                        normalization = 'TSS',
                        transform = 'LOG',
                        cores = 2,
                        plot_summary_plot = FALSE,
                        plot_associations = FALSE,
                        median_comparison_abundance = FALSE,
                        warn_prevalence = FALSE)

    expect_false(mirai::daemons_set())
    expect_true(!is.null(results$fit_data_abundance))
    expect_true(!is.null(results$fit_data_prevalence))
    expect_true(nrow(results$fit_data_abundance$results) > 0)
    expect_true(nrow(results$fit_data_prevalence$results) > 0)

    unlink(tmp_output, recursive = TRUE)
})

test_that("serial and parallel results are consistent", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_test_data()

    results_serial <- maaslin_fit(td$data_in_tss,
                                  td$data_in_tss_log,
                                  td$metadata,
                                  formula = formula('expr ~ var1 + var2'),
                                  random_effects_formula = NULL,
                                  min_abundance = 0,
                                  min_prevalence = 0,
                                  min_variance = 0,
                                  data = td$data_in,
                                  cores = 1,
                                  median_comparison_abundance = FALSE,
                                  median_comparison_prevalence = FALSE)

    results_parallel <- maaslin_fit(td$data_in_tss,
                                    td$data_in_tss_log,
                                    td$metadata,
                                    formula = formula('expr ~ var1 + var2'),
                                    random_effects_formula = NULL,
                                    min_abundance = 0,
                                    min_prevalence = 0,
                                    min_variance = 0,
                                    data = td$data_in,
                                    cores = 2,
                                    median_comparison_abundance = FALSE,
                                    median_comparison_prevalence = FALSE)

    expect_false(mirai::daemons_set())

    serial_res <- results_serial$fit_data_abundance$results
    parallel_res <- results_parallel$fit_data_abundance$results

    serial_res <- serial_res[order(serial_res$feature, serial_res$metadata, serial_res$value), ]
    parallel_res <- parallel_res[order(parallel_res$feature, parallel_res$metadata, parallel_res$value), ]

    expect_equal(nrow(serial_res), nrow(parallel_res))
    expect_equal(serial_res$feature, parallel_res$feature)
    expect_equal(serial_res$metadata, parallel_res$metadata)
    expect_equal(serial_res$coef, parallel_res$coef, tolerance = 1e-6)
    expect_equal(serial_res$pval, parallel_res$pval, tolerance = 1e-6)
})
