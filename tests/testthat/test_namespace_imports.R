library(testthat)
library(maaslin3)

# Verify NAMESPACE declares all imports that R CMD check requires.

test_that("tail is imported from utils", {
    ns_file <- system.file("NAMESPACE", package = "maaslin3")
    ns_lines <- readLines(ns_file)
    tail_import <- grep("importFrom.*utils.*tail", ns_lines, value = TRUE)
    expect_true(length(tail_import) > 0,
                info = "NAMESPACE should importFrom('utils', 'tail')")
})

test_that("dplyr is declared in Suggests for vignette use", {
    desc <- packageDescription("maaslin3")
    expect_true(grepl("\\bdplyr\\b", desc$Suggests),
                info = "dplyr should be in Suggests (used in vignettes)")
})

test_that("internal functions passed to daemons are accessible", {
    fns <- c("check_for_zero_one_obs", "check_missing_first_factor_level",
             "fit_augmented_logistic", "non_augmented", "run_group_models",
             "run_ordered_models", "get_character_cols", "fitting_wrap_up",
             "get_fixed_effects", "make_lm_plot", "make_logistic_plot",
             "optimizers", "optCtrlList")
    
    for (fn_name in fns) {
        obj <- get(fn_name, envir = asNamespace("maaslin3"))
        expect_true(!is.null(obj),
                    info = paste(fn_name, "should exist in maaslin3 namespace"))
    }
})
