library(testthat)
library(maaslin3)

mgx_in <- data.frame('a' = c(1, 2, 0, 4, 5), 
                      'b' = c(2, 3, 4, 5, 6),
                      'c' = c(3, 4, 5, 6, 0))
rownames(mgx_in) <- paste0("sample", c(1:5))

mtx_in <- data.frame('a' = c(1, 2, 3, 4, 5), 
                     'b' = c(2, 3, 4, 5, 0),
                     'c' = c(3, 4, 5, 6, 0))
rownames(mtx_in) <- paste0("sample", c(1:5))

data_in_tss <- data.frame(t(apply(mgx_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0 & mtx_in == 0] <- NA
data_in_tss[data_in_tss == 0] <- min(data_in_tss[data_in_tss > 0], na.rm=T) / 2
data_in_tss_log <- log2(data_in_tss)

expect_that(preprocess_dna_mtx(mgx_in, mtx_in)$dna_table, 
            equals(data_in_tss_log))

data_in_tss <- data.frame(t(apply(mtx_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))

expect_that(preprocess_dna_mtx(mgx_in, mtx_in)$rna_table, 
            equals(data_in_tss))



