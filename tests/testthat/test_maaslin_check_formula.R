library(testthat)
library(maaslin3)

data_in <- data.frame('a' = c(1, 2, 3, 4, 5), 
                    'b' = c(2, 3, 4, 5, 6),
                    'c' = c(3, 4, 5, 6, 7))
rownames(data_in) <- paste0("sample", c(1:5))

metadata <- data.frame('a' = c(1, 2, 3, 4, 5),
                       'b' = c(0, 1, 0, 1, 0),
                       'c' = c('a', 'b', 'c', 'a', 'b'))

expect_that(maaslin_check_formula(data,
                    metadata,
                    input_formula = 'a + b + c')$formula, 
            equals(formula(expr ~ a + b + c)))

expect_error(maaslin_check_formula(data,
                                  metadata,
                                  input_formula = 'strata(b)')$formula)

expect_error(maaslin_check_formula(data,
                                   metadata,
                                   input_formula = ''))

expect_error(maaslin_check_formula(data,
                                metadata,
                                input_formula = '~ a + b + d'))

expect_error(maaslin_check_formula(data,
                                   metadata,
                                   input_formula = NULL))

expect_that(maaslin_check_formula(data,
                                metadata,
                                input_formula = '~ a + b + (1|c)')$formula,
            equals(formula(expr ~ a + b + (1|c))))

expect_that(maaslin_check_formula(data,
                    metadata,
                    input_formula = '~ a + b + (1|c)')$random_effects_formula,
            equals(formula(expr ~ a + b + (1|c))))

expect_that(maaslin_check_formula(data,
                                metadata,
                                input_formula = '~ a + b + (1|c)',
                                feature_specific_covariate_name = 'd')$formula,
            equals(formula(expr ~ d + a + b + (1|c))))

expect_that(maaslin_check_formula(data,
                                  metadata,
                                  input_formula = formula('~ a + b + (1|c)'),
                                  feature_specific_covariate_name = 'd')$formula,
            equals(formula(expr ~ d + a + b + (1|c))))




