library(testthat)
library(maaslin3)

data_in <- data.frame('a' = c(0, 0, 3, 0, 5, 7, 0, 8, 4, 7), 
                      'b' = c(2, 0, 0, 5, 0, 5, 2, 6, 9, 0),
                      'c' = c(3, 4, 0, 0, 7, 2, 7, 4, 7, 0))
rownames(data_in) <- paste0("sample", c(1:10))

metadata <- data.frame('var1' = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       'var2' = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
rownames(metadata) <- paste0("sample", c(1:10))

data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                  FUN = function(x){x / sum(x)})))
data_in_tss[data_in_tss == 0] <- NA
data_in_tss_log <- log2(data_in_tss)

results <- maaslin_fit(data_in_tss,
            data_in_tss_log,
            metadata,
            formula = formula('expr ~ var1 + var2'),
            random_effects_formula = NULL, 
            min_abundance = 0, 
            min_prevalence = 0, 
            min_variance = 0,
            data = data_in, 
            median_comparison_abundance = FALSE)

data_in_tss_log$sample <- rownames(data_in_tss_log)
metadata$sample <- rownames(metadata)
data_merged <- merge(data_in_tss_log, metadata, by = "sample", all = TRUE)
rownames(data_merged) <- data_merged$sample
data_merged$sample <- NULL

growing_df <- data.frame()
for (feature in c('a', 'b', 'c')) {
    summary_out <- summary(lm(formula(paste0(feature, '~ var1 + var2')), data_merged))$coef
    summary_out <- data.frame(summary_out, check.names = F)
    summary_out$variable <- rownames(summary_out)
    summary_out <- summary_out[-1,] # Drop intercept
    summary_out$`t value` <- NULL
    colnames(summary_out) <- c("coef", "stderr", "pval_individual", 'metadata')
    summary_out$name <- summary_out$value <- summary_out$metadata
    summary_out$feature <- feature
    summary_out$error <- NA_character_
    summary_out$model <- 'linear'
    summary_out$N <- nrow(data_merged)
    summary_out$N_not_zero <- sum(!is.na(data_merged[,feature]))
    growing_df <- rbind(growing_df, summary_out)
}

for (feature in c('a', 'b', 'c')) {
    data_subset <- data_merged[,c(feature, 'var1', 'var2')]
    data_subset[,1] <- ifelse(is.na(data_subset[,1]), 0, 1)
    orig_length <- nrow(data_subset)
    data_subset <- rbind(data_subset, data_subset, data_subset)
    data_subset[(orig_length + 1):(2 * orig_length),1] <- 1
    data_subset[(2 * orig_length + 1):(3 * orig_length),1] <- 0
    summary_out <- tryCatch({
        withCallingHandlers({
                summary(glm(formula(paste0(feature, ' ~ var1 + var2')), 
                            data_subset, family = 'binomial', 
                            weights = c(rep(1, orig_length),
                                        rep(2 / 2 / orig_length, 
                                            2 * orig_length))))$coef
            },
            warning = function(w) {
                invokeRestart("muffleWarning")
            }
        )
    })
    
    summary_out <- data.frame(summary_out, check.names = F)
    summary_out$variable <- rownames(summary_out)
    summary_out <- summary_out[-1,] # Drop intercept
    summary_out$`z value` <- NULL
    colnames(summary_out) <- c("coef", "stderr", "pval_individual", 'metadata')
    summary_out$name <- summary_out$value <- summary_out$metadata
    summary_out$feature <- feature
    summary_out$error <- NA_character_
    summary_out$model <- 'logistic'
    summary_out$N <- nrow(data_merged)
    summary_out$N_not_zero <- sum(!is.na(data_merged[,feature]))
    growing_df <- rbind(growing_df, summary_out)
}
growing_df$qval_individual <- p.adjust(growing_df$pval_individual, 'BH')

merged_df <- dplyr::full_join(growing_df[growing_df$model == 'linear',],
                       growing_df[growing_df$model == 'logistic',], 
                       by = c('metadata', 'value', 'name', 'feature', 'N', 
                              'N_not_zero'))
pval_joint <- pbeta(pmin(merged_df$pval_individual.x, 
                        merged_df$pval_individual.y), 1, 2)
qval_joint <- p.adjust(pval_joint, 'BH')
growing_df$pval_joint <- rep(pval_joint, 2)
growing_df$qval_joint <- rep(qval_joint, 2)

col_order <- c("feature", "metadata", "value", "name", "coef", "stderr", 
               "pval_individual", "qval_individual", "pval_joint", 
               "qval_joint", "error", "model", "N", "N_not_zero")
growing_df <- growing_df[,col_order]
growing_df <- growing_df[order(growing_df$qval_individual),]

linear_piece <- growing_df[growing_df$model == 'linear',]
logistic_piece <- growing_df[growing_df$model == 'logistic',]
rownames(linear_piece) <- NULL
rownames(logistic_piece) <- NULL
rownames(results$fit_data_abundance$results) <- NULL
rownames(results$fit_data_prevalence$results) <- NULL

expect_that(results$fit_data_abundance$results,
            equals(linear_piece))

expect_that(results$fit_data_prevalence$results,
            equals(logistic_piece))
