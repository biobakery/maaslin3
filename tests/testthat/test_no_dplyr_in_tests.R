library(testthat)
library(maaslin3)

test_that("test files do not use dplyr in active code", {
    test_dir <- system.file("tests", "testthat", package = "maaslin3")

    if (nchar(test_dir) == 0) {
        test_dir <- file.path(getwd(), "tests", "testthat")
    }

    test_files <- list.files(test_dir, pattern = "\\.R$", full.names = TRUE)
    test_files <- test_files[basename(test_files) != "test_no_dplyr_in_tests.R"]

    for (tf in test_files) {
        lines <- readLines(tf, warn = FALSE)
        code_lines <- lines[!grepl("^\\s*#", lines)]
        pattern <- paste0("dplyr", "::", "|library\\(dplyr\\)")
        dplyr_calls <- grep(pattern, code_lines, value = TRUE)
        expect_equal(length(dplyr_calls), 0,
                     info = paste("File", basename(tf),
                                  "should not reference dplyr in code:",
                                  paste(dplyr_calls, collapse = "; ")))
    }
})
