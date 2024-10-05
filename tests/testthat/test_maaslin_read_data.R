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

unscaled <- data.frame('a' = c(5, 4, 3, 2, 1))
rownames(unscaled) <- paste0("sample", c(1:5))

metadata <- data.frame('a' = c(1, 2, 3, 4, 5),
                       'b' = c(0, 1, 0, 1, 0),
                       'c' = c('a', 'b', 'c', 'a', 'b'))
rownames(metadata) <- paste0("sample", c(1:5))

output_tmp <- tempfile()
dir.create(output_tmp)
write.table(data_in, file.path(output_tmp, 'data.tsv'), 
            row.names = T, sep = '\t')
write.table(metadata, file.path(output_tmp, 'metadata.tsv'), 
            row.names = T, sep = '\t')
write.table(covar_data, file.path(output_tmp, 'covar.tsv'), 
            row.names = T, sep = '\t')
write.table(unscaled, file.path(output_tmp, 'unscaled.tsv'), 
            row.names = T, sep = '\t')

data_in_new <- maaslin_read_data(input_data = file.path(output_tmp, 'data.tsv'), 
                input_metadata = file.path(output_tmp, 'metadata.tsv'), 
                feature_specific_covariate = file.path(output_tmp, 'covar.tsv'), 
                unscaled_abundance = file.path(output_tmp, 'unscaled.tsv'))

expect_that(data_in_new$data, equals(data_in))
expect_that(data_in_new$metadata, equals(metadata))
expect_that(data_in_new$feature_specific_covariate, equals(covar_data))
expect_that(data_in_new$unscaled_abundance, equals(unscaled))

