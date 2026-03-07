library(testthat)
library(maaslin3)

# The inherits() check for BioC DataFrame must test input_metadata, not the
# already-converted local variable 'metadata'. The old code used
# inherits(metadata, 'DataFrame') which would always be FALSE after the
# data.frame conversion above it, causing DataFrame inputs to hit the
# stop() branch.

test_that("maaslin_read_data handles a regular data.frame for metadata", {
    data_in <- data.frame(a = c(1, 2, 3, 4, 5),
                          b = c(5, 4, 3, 2, 1))
    rownames(data_in) <- paste0("s", 1:5)

    meta <- data.frame(var1 = rnorm(5))
    rownames(meta) <- paste0("s", 1:5)

    result <- maaslin_read_data(data_in, meta, NULL, NULL)
    expect_true(is.data.frame(result$metadata))
    expect_equal(nrow(result$metadata), 5)
})

test_that("maaslin_read_data rejects non-data-frame non-file metadata", {
    data_in <- data.frame(a = c(1, 2, 3), b = c(3, 2, 1))
    rownames(data_in) <- paste0("s", 1:3)

    expect_error(
        maaslin_read_data(data_in, list(x = 1), NULL, NULL),
        "neither a file nor a data frame"
    )
})
