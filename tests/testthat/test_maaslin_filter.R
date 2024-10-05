library(testthat)
library(maaslin3)

data_in <- data.frame('a' = c(1, 2, 3, 4, 5), 
                      'b' = c(2, 3, 4, 5, 6),
                      'c' = c(3, 4, 5, 6, 7))
rownames(data_in) <- paste0("sample", c(1:5))

metadata <- data.frame('a' = c(1, 2, 3, 4, 5),
                       'b' = c(0, 1, 0, 1, 0),
                       'c' = c('a', 'b', 'c', 'a', 'b'))
rownames(metadata) <- paste0("sample", c(1:5))

output_tmp <- tempfile()
expect_that(maaslin_filter(data_in, output_tmp),
            equals(data_in))

expect_that(maaslin_filter(data_in, output_tmp, 
                           min_abundance = 3, min_prevalence = 1/2),
            equals(data_in[,-1]))

expect_that(maaslin_filter(data_in, output_tmp, zero_threshold = 3,
                           min_prevalence = 1/2),
            equals(data_in[,-1]))

data_in <- data.frame('a' = c(1, 1, 1, 1, 1), 
                      'b' = c(2, 3, 4, 5, 6),
                      'c' = c(3, 4, 5, 6, 7))
rownames(data_in) <- paste0("sample", c(1:5))

expect_that(maaslin_filter(data_in, output_tmp),
            equals(data_in[,-1]))

unlink(output_tmp, recursive = T)