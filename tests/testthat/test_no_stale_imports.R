library(testthat)
library(maaslin3)

# The NAMESPACE previously imported dplyr::%>% and plyr::mapvalues even though
# neither is used. The DESCRIPTION previously listed pbapply and parallel in
# Imports even though neither is used in active code paths.

test_that("dplyr pipe is not imported", {
    ns_file <- system.file("NAMESPACE", package = "maaslin3")
    ns_lines <- readLines(ns_file)
    dplyr_pipe_lines <- grep('importFrom.*dplyr.*%>%', ns_lines, value = TRUE)
    expect_equal(length(dplyr_pipe_lines), 0,
                 info = "NAMESPACE should not import dplyr pipe")
})

test_that("plyr mapvalues is not imported", {
    ns_file <- system.file("NAMESPACE", package = "maaslin3")
    ns_lines <- readLines(ns_file)
    plyr_lines <- grep('importFrom.*plyr.*mapvalues', ns_lines, value = TRUE)
    expect_equal(length(plyr_lines), 0,
                 info = "NAMESPACE should not import plyr::mapvalues")
})

test_that("pbapply is not in Imports", {
    desc <- packageDescription("maaslin3")
    expect_false(grepl("pbapply", desc$Imports),
                 info = "pbapply should not be in Imports")
})

test_that("parallel is not in Imports", {
    desc <- packageDescription("maaslin3")
    expect_false(grepl("\\bparallel\\b", desc$Imports),
                 info = "parallel should not be in Imports")
})
