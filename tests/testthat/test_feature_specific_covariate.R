library(testthat)
library(maaslin3)

setup_fsc_data <- function() {
    set.seed(42)
    n <- 60
    data_in <- data.frame(
        a = rpois(n, lambda = 5),
        b = rpois(n, lambda = 3),
        c = rpois(n, lambda = 4)
    )
    rownames(data_in) <- paste0("sample", 1:n)

    metadata <- data.frame(
        var1 = rnorm(n),
        var2 = c(rep(0, 30), rep(1, 30)),
        var_grp = factor(rep(c("A", "B", "C"), each = 20),
                         levels = c("A", "B", "C")),
        var_ord = factor(rep(c("low", "med", "high"), each = 20),
                         levels = c("low", "med", "high"),
                         ordered = TRUE),
        var_re = factor(rep(paste0("g", 1:6), each = 10))
    )
    rownames(metadata) <- paste0("sample", 1:n)

    data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1,
                                      FUN = function(x) x / sum(x))))
    data_in_tss[data_in_tss == 0] <- NA
    data_in_tss_log <- log2(data_in_tss)

    set.seed(7)
    dna <- data.frame(
        a = rnorm(n),
        b = rnorm(n),
        c = rnorm(n)
    )
    rownames(dna) <- paste0("sample", 1:n)

    list(data_in = data_in, metadata = metadata,
         data_in_tss = data_in_tss, data_in_tss_log = data_in_tss_log,
         dna = dna)
}

fit_fsc <- function(td, formula, cores, random_effects_formula = NULL) {
    maaslin_fit(td$data_in_tss,
                td$data_in_tss_log,
                td$metadata,
                formula = formula,
                random_effects_formula = random_effects_formula,
                min_abundance = 0,
                min_prevalence = 0,
                min_variance = 0,
                data = td$data_in,
                cores = cores,
                feature_specific_covariate = td$dna,
                feature_specific_covariate_name = "DNA",
                feature_specific_covariate_record = TRUE,
                median_comparison_abundance = FALSE,
                median_comparison_prevalence = FALSE)
}

sort_results <- function(res) {
    res[order(res$feature, res$metadata, res$value), ]
}

test_that("feature_specific_covariate runs serially and matches manual lm", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    results <- fit_fsc(td, formula("expr ~ var1 + var2 + DNA"), cores = 1)

    abundance <- results$fit_data_abundance$results
    expect_true(!is.null(abundance))
    expect_true("DNA" %in% abundance$metadata)
    expect_true(all(is.na(abundance$error)))

    for (feature in c("a", "b", "c")) {
        manual_df <- data.frame(
            expr = td$data_in_tss_log[[feature]],
            var1 = td$metadata$var1,
            var2 = td$metadata$var2,
            DNA = td$dna[[feature]]
        )
        manual_coef <- summary(
            lm(expr ~ var1 + var2 + DNA, manual_df))$coef["DNA", "Estimate"]

        fitted_coef <- abundance$coef[abundance$feature == feature &
                                          abundance$metadata == "DNA"]
        expect_equal(length(fitted_coef), 1)
        expect_equal(fitted_coef, manual_coef, tolerance = 1e-6)
    }
})

test_that("feature_specific_covariate parallel results match serial", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    serial <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + var2 + DNA"),
                cores = 1)$fit_data_abundance$results)
    parallel <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + var2 + DNA"),
                cores = 2)$fit_data_abundance$results)

    expect_false(mirai::daemons_set())
    expect_true("DNA" %in% parallel$metadata)
    expect_true(all(is.na(parallel$error)))
    expect_equal(nrow(serial), nrow(parallel))
    expect_equal(serial$metadata, parallel$metadata)
    expect_equal(serial$coef, parallel$coef, tolerance = 1e-6)
})

test_that("feature_specific_covariate works with random effects", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    serial <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + var2 + DNA + (1|var_re)"),
                cores = 1,
                random_effects_formula =
                    formula("expr ~ (1|var_re)"))$fit_data_abundance$results)
    parallel <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + var2 + DNA + (1|var_re)"),
                cores = 2,
                random_effects_formula =
                    formula("expr ~ (1|var_re)"))$fit_data_abundance$results)

    expect_false(mirai::daemons_set())
    expect_true("DNA" %in% serial$metadata)
    expect_true("DNA" %in% parallel$metadata)
    expect_true(all(is.na(serial$error)))
    expect_equal(nrow(serial), nrow(parallel))
    expect_equal(serial$metadata, parallel$metadata)
    expect_equal(serial$coef, parallel$coef, tolerance = 1e-6)
})

test_that("feature_specific_covariate works with a group effect", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    serial <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + group(var_grp)"),
                cores = 1)$fit_data_abundance$results)
    parallel <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + group(var_grp)"),
                cores = 2)$fit_data_abundance$results)

    expect_false(mirai::daemons_set())
    expect_true("DNA" %in% serial$metadata)
    expect_true("var_grp" %in% serial$metadata)
    expect_true("DNA" %in% parallel$metadata)
    expect_true("var_grp" %in% parallel$metadata)
    expect_equal(nrow(serial), nrow(parallel))
    expect_equal(serial$metadata, parallel$metadata)
    expect_equal(serial$coef, parallel$coef, tolerance = 1e-6)
})

test_that("feature_specific_covariate works with an ordered effect", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    serial <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + ordered(var_ord)"),
                cores = 1)$fit_data_abundance$results)
    parallel <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + ordered(var_ord)"),
                cores = 2)$fit_data_abundance$results)

    expect_false(mirai::daemons_set())
    expect_true("DNA" %in% serial$metadata)
    expect_true("var_ord" %in% serial$metadata)
    expect_true("DNA" %in% parallel$metadata)
    expect_true("var_ord" %in% parallel$metadata)
    expect_equal(nrow(serial), nrow(parallel))
    expect_equal(serial$metadata, parallel$metadata)
    expect_equal(serial$coef, parallel$coef, tolerance = 1e-6)
})

test_that("feature_specific_covariate works with a strata effect", {
    if (mirai::daemons_set()) mirai::daemons(0)
    Sys.sleep(0.5)

    td <- setup_fsc_data()

    serial <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + strata(var_grp)"),
                cores = 1)$fit_data_abundance$results)
    parallel <- sort_results(
        fit_fsc(td, formula("expr ~ var1 + DNA + strata(var_grp)"),
                cores = 2)$fit_data_abundance$results)

    expect_false(mirai::daemons_set())
    expect_true("DNA" %in% serial$metadata)
    expect_true("DNA" %in% parallel$metadata)
    expect_equal(nrow(serial), nrow(parallel))
    expect_equal(serial$metadata, parallel$metadata)
    expect_equal(serial$coef, parallel$coef, tolerance = 1e-6)
})
