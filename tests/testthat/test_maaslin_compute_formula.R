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

expect_that(maaslin_compute_formula(data,
              metadata,
              fixed_effects = c('a', 'b', 'c'),
              random_effects = NULL,
              group_effects = NULL,
              ordered_effects = NULL,
              strata_effects = NULL,
              feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + c)))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = NULL,
                                    random_effects = NULL,
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + c)))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a', 'b'),
                                    random_effects = 'c',
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + (1|c))))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a', 'b'),
                                    random_effects = NULL,
                                    group_effects = 'c',
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + group(c))))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a', 'b'),
                                    random_effects = NULL,
                                    group_effects = NULL,
                                    ordered_effects = 'c',
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + ordered(c))))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a', 'b'),
                                    random_effects = NULL,
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = 'c',
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + strata(c))))

expect_that(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a,b,c'),
                                    random_effects = NULL,
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL)$formula, 
            equals(formula(expr ~ a + b + c)))

expect_error(maaslin_compute_formula(data,
                                    metadata,
                                    fixed_effects = c('a,b,c'),
                                    random_effects = "d",
                                    group_effects = NULL,
                                    ordered_effects = NULL,
                                    strata_effects = NULL,
                                    feature_specific_covariate_name = NULL))

expect_error(maaslin_compute_formula(data,
                                     metadata,
                                     fixed_effects = character(0),
                                     random_effects = character(0),
                                     group_effects = character(0),
                                     ordered_effects = character(0),
                                     strata_effects = character(0),
                                     feature_specific_covariate_name = NULL))

