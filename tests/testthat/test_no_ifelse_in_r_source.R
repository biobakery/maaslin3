library(testthat)
library(maaslin3)

# After migration to data.table::fifelse and if/else, the only remaining
# base::ifelse calls should be on matrices (where fifelse doesn't work).
# This test checks the installed package's decompiled code for stray ifelse.

get_function_body_text <- function(fn) {
    paste(deparse(body(fn)), collapse = "\n")
}

test_that("fit.model does not use base ifelse", {
    fn_body <- get_function_body_text(maaslin3:::fit.model)
    ifelse_matches <- gregexpr("\\bifelse\\(", fn_body)[[1]]
    n_matches <- if (ifelse_matches[1] == -1) 0 else length(ifelse_matches)
    expect_equal(n_matches, 0,
                 info = "fit.model should not use base ifelse()")
})

test_that("key viz functions do not use base ifelse", {
    for (fn_name in c("maaslin_plot_results_from_output", 
                       "maaslin_plot_results")) {
        fn <- getFromNamespace(fn_name, "maaslin3")
        fn_body <- get_function_body_text(fn)
        ifelse_matches <- gregexpr("\\bifelse\\(", fn_body)[[1]]
        n_matches <- if (ifelse_matches[1] == -1) 0 else length(ifelse_matches)
        expect_equal(n_matches, 0,
                     info = paste(fn_name, "should not use base ifelse()"))
    }
})
