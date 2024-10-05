library(testthat)
library(maaslin3)

metadata <- data.frame('a' = c(1, 2, 3, 4, 5),
                       'b' = c(0, 1, 0, 1, 0),
                       'c' = c('a', 'b', 'c', 'a', 'b'))
rownames(metadata) <- paste0("sample", c(1:5))

expect_error(maaslin_process_metadata(metadata = metadata, 
                         formula = formula('~ a + b + c'), 
                         standardize = T))

metadata_out <- metadata
metadata_out$a <- scale(metadata_out$a)
metadata_out$b <- scale(metadata_out$b)
metadata_out$c <- factor(metadata_out$c)
expect_that(maaslin_process_metadata(metadata = metadata, 
                                      formula = formula('~ a + b + c'), 
                                      reference = 'c,a',
                                      standardize = T),
            equals(metadata_out))

metadata$c <- factor(metadata$c)
expect_that(maaslin_process_metadata(metadata = metadata, 
                                     formula = formula('~ a + b + c'), 
                                     standardize = F),
            equals(metadata))

expect_that(maaslin_process_metadata(metadata = metadata, 
                             formula = formula('~ a + (1|b) + strata(c)'), 
                             standardize = F),
            equals(metadata))



