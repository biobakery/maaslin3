library(testthat)
library(maaslin3)

# R CMD check warns when a package uses ::: to access its own objects.
# All internal functions should be called directly or passed to daemons
# via mirai::everywhere().

test_that("no maaslin3::: self-calls in R source files", {
    r_dir <- system.file("R", package = "maaslin3")
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    
    all_violations <- character(0)
    for (f in r_files) {
        lines <- readLines(f)
        hits <- grep("maaslin3:::", lines, value = TRUE)
        hits <- hits[!grepl("^\\s*#", hits)]
        if (length(hits) > 0) {
            all_violations <- c(all_violations,
                                paste0(basename(f), ": ", trimws(hits)))
        }
    }
    expect_equal(length(all_violations), 0,
                 info = paste("Found maaslin3::: self-calls:",
                              paste(all_violations, collapse = "\n")))
})

test_that("no library() or require() calls in R source files", {
    r_dir <- system.file("R", package = "maaslin3")
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    
    all_violations <- character(0)
    for (f in r_files) {
        lines <- readLines(f)
        hits <- grep("\\b(library|require)\\s*\\(", lines, value = TRUE)
        hits <- hits[!grepl("^\\s*#", hits)]
        hits <- hits[!grepl("requireNamespace", hits)]
        if (length(hits) > 0) {
            all_violations <- c(all_violations,
                                paste0(basename(f), ": ", trimws(hits)))
        }
    }
    expect_equal(length(all_violations), 0,
                 info = paste("Found library/require calls:",
                              paste(all_violations, collapse = "\n")))
})
