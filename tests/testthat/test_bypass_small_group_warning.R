library(testthat)
library(maaslin3)

set.seed(42)

data_in <- data.frame(
    'a' = c(0, 3, 5, 7, 0, 8, 4, 7, 2, 6, 3, 5),
    'b' = c(2, 0, 5, 0, 5, 2, 6, 9, 0, 4, 7, 1),
    'c' = c(3, 4, 0, 7, 2, 7, 4, 7, 0, 5, 8, 3)
)
rownames(data_in) <- paste0("sample", 1:12)

metadata <- data.frame(
    'var1' = rnorm(12),
    'group_var' = rep(paste0("g", 1:4), each = 3)
)
rownames(metadata) <- paste0("sample", 1:12)

data_in_tss <- data.frame(t(apply(
    data_in, MARGIN = 1,
    FUN = function(x) { x / sum(x) }
)))
data_in_tss[data_in_tss == 0] <- NA
data_in_tss_log <- log2(data_in_tss)

formula_obj <- formula('expr ~ var1 + (1 | group_var)')
random_effects_formula <- formula('expr ~ (1 | group_var)')

small_group_msg <- paste0(
    "<4 average observations per random effect group often inflates ",
    "coefficients and deflates p-values: consider setting ",
    "small_random_effects=TRUE and see tutorial")

test_that("small group warning appears when bypass is FALSE", {
    results <- maaslin_fit(
        data_in_tss,
        data_in_tss_log,
        metadata,
        formula = formula_obj,
        random_effects_formula = random_effects_formula,
        min_abundance = 0,
        min_prevalence = 0,
        min_variance = 0,
        data = data_in,
        median_comparison_abundance = FALSE,
        bypass_small_group_warning = FALSE
    )

    expect_true(
        any(grepl("<4 average observations per random effect",
                  results$fit_data_prevalence$results$error))
    )
})

test_that("small group warning is bypassed when bypass is TRUE", {
    results <- maaslin_fit(
        data_in_tss,
        data_in_tss_log,
        metadata,
        formula = formula_obj,
        random_effects_formula = random_effects_formula,
        min_abundance = 0,
        min_prevalence = 0,
        min_variance = 0,
        data = data_in,
        median_comparison_abundance = FALSE,
        bypass_small_group_warning = TRUE
    )

    has_small_group_warning <- any(grepl(
        "<4 average observations per random effect",
        results$fit_data_prevalence$results$error
    ))
    expect_false(has_small_group_warning)
})
