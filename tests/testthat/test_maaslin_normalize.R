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

data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))

output_tmp <- tempfile()
expect_that(maaslin_normalize(data_in, output_tmp),
            equals(data_in_tss))

expect_that(maaslin_normalize(data_in, output_tmp, normalization = 'NONE'),
            equals(data_in))

unscaled <- data.frame('total' = c(5, 4, 3, 2, 1))
rownames(unscaled) <- paste0("sample", c(1:5))
scaled_data <- data.frame(apply(data_in_tss, MARGIN = 2, 
                    FUN = function(x){x * unscaled$total}))
expect_that(maaslin_normalize(data_in, output_tmp, 
                              normalization = 'TSS', 
                              unscaled_abundance = unscaled),
            equals(scaled_data))

unscaled <- data.frame('a' = c(5, 4, 3, 2, 1))
rownames(unscaled) <- paste0("sample", c(1:5))
scaled_data <- data.frame(apply(data_in_tss, MARGIN = 2, 
                        FUN = function(x){x * unscaled$a / data_in_tss$a}))
scaled_data$a <- NULL
expect_that(maaslin_normalize(data_in, output_tmp, 
                              normalization = 'TSS', 
                              unscaled_abundance = unscaled),
            equals(scaled_data))

data_in <- data.frame('a' = c(0, 0, 3, 4, 5), 
                      'b' = c(2, 0, 0, 5, 6),
                      'c' = c(3, 4, 0, 0, 7))
rownames(data_in) <- paste0("sample", c(1:5))
data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0] <- NA
expect_that(maaslin_normalize(data_in, output_tmp),
            equals(data_in_tss))

unlink(output_tmp, recursive = T)
