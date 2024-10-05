library(testthat)
library(maaslin3)

expect_error(
    maaslin_check_arguments(unscaled_abundance = data.frame('total' = c(1,2,3)), 
                            median_comparison_abundance = TRUE)
)

expect_error(
    maaslin_check_arguments(feature_specific_covariate = NULL,
                            feature_specific_covariate_name = "something", 
                            feature_specific_covariate_record = NULL)
)

expect_error(
    maaslin_check_arguments(normalization = 'nonsense')
)

expect_error(
    maaslin_check_arguments(normalization = 'NONE',
                            unscaled_abundance = data.frame('total' = c(1,2,3)),
                            median_comparison_abundance = FALSE)
)

expect_error(
    maaslin_check_arguments(transform = 'nonsense')
)

expect_error(
    maaslin_check_arguments(correction = 'nonsense')
)

expect_error(
    maaslin_check_arguments(normalization = 'CLR',
                            transform = 'LOG')
)

expect_error(
    maaslin_check_arguments(transform = 'PLOG',
                            zero_threshold = 0.1)
)

expect_error(
    maaslin_check_arguments(normalization = 'TSS',
                            transform = 'NONE')
)

expect_error(
    maaslin_check_arguments(evaluate_only = 'abundance')
)

expect_error(
    maaslin_check_arguments(evaluate_only = 'nonsense')
)

expect_error(
    maaslin_check_arguments(transform = 'PLOG')
)

expect_that(maaslin_check_arguments(), equals(NULL))

expect_that(maaslin_check_arguments(normalization = 'TSS', 
                                    transform = 'LOG', 
                                    zero_threshold = 0.5,
                                    median_comparison_abundance = TRUE), 
            equals(NULL))
