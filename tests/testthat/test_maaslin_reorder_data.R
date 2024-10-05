library(testthat)
library(maaslin3)

data_in <- data.frame('a' = c(1, 2, 3, 4, 5), 
                      'b' = c(2, 3, 4, 5, 6),
                      'c' = c(3, 4, 5, 6, 7))
rownames(data_in) <- paste0("sample", c(1:5))

covar_data <- data.frame('a' = c(5, 4, 3, 2, 1), 
                         'b' = c(6, 5, 4, 3, 2),
                         'c' = c(7, 6, 5, 4, 3))
rownames(covar_data) <- paste0("sample", c(1:5))

unscaled <- data.frame('total' = c(5, 4, 3, 2, 1))
rownames(unscaled) <- paste0("sample", c(1:5))

metadata <- data.frame('var1' = c(1, 2, 3, 4, 5),
                       'var2' = c(0, 1, 0, 1, 0),
                       'var3' = c('a', 'b', 'c', 'a', 'b'))
rownames(metadata) <- paste0("sample", c(1:5))

data_in_new <- maaslin_reorder_data(data_in, metadata, covar_data, unscaled)

expect_that(data_in_new$data, equals(data_in))
expect_that(data_in_new$metadata, equals(metadata))
expect_that(data_in_new$feature_specific_covariate, equals(covar_data))
expect_that(data_in_new$unscaled_abundance, equals(unscaled))

data_in_new <- maaslin_reorder_data(t(data_in), metadata, 
                                    t(covar_data), unscaled)

expect_that(data_in_new$data, equals(data_in))
expect_that(data_in_new$metadata, equals(metadata))
expect_that(data_in_new$feature_specific_covariate, equals(covar_data))
expect_that(data_in_new$unscaled_abundance, equals(unscaled))

data_in_new <- maaslin_reorder_data(data_in[-1,], metadata, 
                                    t(covar_data), unscaled)

expect_that(data_in_new$data, equals(data_in[-1,]))
expect_that(data_in_new$metadata, equals(metadata[-1,]))
expect_that(data_in_new$feature_specific_covariate, equals(covar_data[-1,]))
expect_that(data_in_new$unscaled_abundance, equals(unscaled[-1, ,drop=F]))



