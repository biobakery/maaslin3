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
data_in_tss_log <- log2(data_in_tss)
output_tmp <- tempfile()
expect_that(maaslin_transform(data_in_tss, output_tmp),
            equals(data_in_tss_log))

data_in <- data.frame('a' = c(0, 0, 3, 4, 5), 
                      'b' = c(2, 0, 0, 5, 6),
                      'c' = c(3, 4, 0, 0, 7))
rownames(data_in) <- paste0("sample", c(1:5))
data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0] <- NA
data_in_tss_log <- log2(data_in_tss)
expect_that(maaslin_transform(data_in_tss, output_tmp),
            equals(data_in_tss_log))

data_in_tss[is.na(data_in_tss)] <- min(data_in_tss[!is.na(data_in_tss)]) / 2
data_in_tss_log <- log2(data_in_tss)
data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
expect_that(maaslin_transform(data_in_tss, output_tmp, transform = 'PLOG'),
            equals(data_in_tss_log))

expect_that(maaslin_transform(data_in_tss, output_tmp, transform = 'NONE'),
            equals(data_in_tss))

unlink(output_tmp, recursive = T)
