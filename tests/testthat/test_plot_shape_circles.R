library(testthat)
library(maaslin3)

# Association plots must use circular points (shape 21). A previous change
# switched to squares (shape 22) for rendering speed, but this was reverted
# to preserve the original visual appearance.

test_that("association plots use circular points (shape 21)", {
    output_tmp <- tempfile()
    
    data_in <- data.frame('a' = c(0, 0, 0, 0, 0, 0, 0, 0, 7, 5, 8, 4, 7, 8, 5, 8), 
                          'b' = c(2, 0, 0, 5, 0, 5, 3, 6, 5, 2, 6, 9, 0, 8, 3, 7))
    rownames(data_in) <- paste0("sample", seq(nrow(data_in)))
    
    metadata <- data.frame('var1' = c(rep(0, 8), rep(1, 7), 2))
    rownames(metadata) <- paste0("sample", seq(nrow(data_in)))
    
    data_in_tss <- data.frame(t(apply(data_in, MARGIN = 1, 
                                      FUN = function(x) x / sum(x))))
    data_in_tss[data_in_tss == 0] <- NA
    data_in_tss_log <- log2(data_in_tss)
    
    results <- maaslin_fit(data_in_tss,
                           data_in_tss_log,
                           metadata,
                           formula = formula('expr ~ var1'),
                           random_effects_formula = NULL, 
                           min_abundance = 0, 
                           min_prevalence = 0, 
                           min_variance = 0,
                           data = data_in, 
                           median_comparison_abundance = FALSE)
    
    plot_out <- maaslin_plot_results(output = output_tmp,
                                    transformed_data = data_in_tss_log,
                                    unstandardized_metadata = metadata,
                                    fit_data_abundance = results$fit_data_abundance, 
                                    fit_data_prevalence = results$fit_data_prevalence, 
                                    normalization = 'TSS',
                                    transform = 'LOG',
                                    median_comparison_abundance = FALSE,
                                    max_significance = 0.1,
                                    save_plots_rds = TRUE)
    
    gg <- plot_out$assocation_plots[[1]]
    expect_is(gg, 'ggplot')
    
    point_layers <- Filter(function(layer) {
        inherits(layer$geom, "GeomPoint")
    }, gg$layers)
    
    shapes_used <- vapply(point_layers, function(layer) {
        as.numeric(layer$aes_params$shape %||% NA_real_)
    }, numeric(1))
    
    expect_true(all(shapes_used[!is.na(shapes_used)] == 21),
                info = "All point geoms should use shape 21 (circles)")
    
    unlink(output_tmp, recursive = TRUE)
})
