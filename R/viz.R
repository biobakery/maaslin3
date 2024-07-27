#!/usr/bin/env Rscript
###############################################################################
# MaAsLin3 visualizations

# Copyright (c) 2018 Harvard School of Public Health

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# Load libararies
for (lib in c('dplyr', 'ggplot2', 'viridis', "grid", 'RColorBrewer', 'patchwork', 'scales')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# MaAsLin3 theme based on Nature journal requirements
nature_theme <- function(x_axis_labels, y_label) {
    # set default text format based on categorical and length
    angle = NULL
    hjust = NULL
    size = 8
    if (max(nchar(x_axis_labels), na.rm=TRUE) > 5) {
        angle = 45
        hjust = 1
        size = 6
    }
    axis_title_size = 10
    if (nchar(y_label) > 15) {
        axis_title_size = 8
    }
    if (nchar(y_label) > 25) {
        axis_title_size = 6
    }
    return ( ggplot2::theme_bw() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = size, vjust = 1, hjust = hjust, angle = angle),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title = ggplot2::element_text(size = axis_title_size),
        plot.title = ggplot2::element_text(size = 7, face = 'bold'),
        legend.title = ggplot2::element_text(size = 6, face = 'bold'),
        legend.text = ggplot2::element_text(size = 6),
        axis.line = ggplot2::element_line(colour = 'black', linewidth = .25),
        axis.line.x = ggplot2::element_line(colour = 'black', linewidth = .25),
        axis.line.y = ggplot2::element_line(colour = 'black', linewidth = .25),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())
   )
}

# MaAsLin3 summary_plot function for overall view of associations
maaslin3_summary_plot <-
  function(
    merged_results,
    summary_plot_file,
    figures_folder,
    first_n = 30,
    max_significance = 0.1,
    coef_plot_vars = NULL,
    heatmap_vars = NULL,
    median_comparison_abundance = FALSE,
    median_comparison_prevalence = FALSE) {
    
    if (first_n > 200) {
      logging::logerror(
        paste("At most 200 features can be plotted in the heatmap. Please choose a smaller first_n."))
      return()
    }
    
    # Preprocessing
    merged_results <- merged_results[is.na(merged_results$error) & 
                                       !is.na(merged_results$qval_individual) &
                                       !is.na(merged_results$coef),]
    if (nrow(merged_results) == 0) {
      logging::loginfo(
        paste("No associtions were without errors. No summary plot generated."))
      return(NULL)
    }
    merged_results$model <- ifelse(merged_results$model == 'LM', 'Abundance', 'Prevalence')
    merged_results$full_metadata_name <- 
      ifelse(merged_results$metadata == merged_results$value,
             merged_results$metadata,
             paste0(merged_results$metadata, ' ', merged_results$value))
    
    median_df <- merged_results %>%
      dplyr::group_by(.data$full_metadata_name, .data$model) %>%
      dplyr::summarize(median_val = median(.data$coef), .groups = 'drop')
    
    if (!median_comparison_abundance) {
      median_df$median_val[median_df$model == 'Abundance'] <- 0
    }
    if (!median_comparison_prevalence) {
      median_df$median_val[median_df$model == 'Prevalence'] <- 0
    }
    
    if (!is.null(coef_plot_vars) | !is.null(heatmap_vars)) {
      if (any(!c(coef_plot_vars, heatmap_vars) %in% unique(merged_results$full_metadata_name))) {
        logging::loginfo(
          paste0("The following specified variables were not found in the associations: ", 
                paste0(setdiff(c(coef_plot_vars, heatmap_vars), unique(merged_results$full_metadata_name)), collapse = ', '), collapse = ''))
        logging::loginfo(
          paste0("Available associations: ", 
                 paste0(unique(merged_results$full_metadata_name), collapse = ', '), collapse = ''))
        return(NULL)
      }
      merged_results <- merged_results[merged_results$full_metadata_name %in% 
                                         c(coef_plot_vars, heatmap_vars),]
    }
    
    merged_results_joint_only <- unique(merged_results[,c('feature', 'qval_joint')])
    merged_results_joint_only <- merged_results_joint_only[order(merged_results_joint_only$qval_joint),]
    if (length(unique(merged_results_joint_only$feature)) < first_n) {
      first_n <- length(unique(merged_results_joint_only$feature))
    }
    signif_taxa <- unique(merged_results_joint_only$feature)[1:first_n]

    merged_results_sig <- merged_results %>%
      dplyr::filter(.data$feature %in% signif_taxa)
    
    # order feature
    ord_feature <- with(merged_results_sig, reorder(feature, coef))
    ord_feature <- levels(ord_feature)
    
    merged_results_sig$feature <- factor(merged_results_sig$feature, levels = ord_feature)
    
    if(is.null(coef_plot_vars)) {
      mean_log_qval <- merged_results_sig %>%
        dplyr::group_by(.data$full_metadata_name) %>%
        dplyr::summarise(mean_value = mean(log(.data$qval_joint), na.rm = TRUE))
      
      coef_plot_vars <- mean_log_qval$full_metadata_name[order(mean_log_qval$mean_value)]
      coef_plot_vars <- setdiff(coef_plot_vars, heatmap_vars)
      if (length(coef_plot_vars) > 0) {
        coef_plot_vars <- coef_plot_vars[1:min(2, length(coef_plot_vars))]
      }
    }
    
    if(is.null(heatmap_vars)) {
      mean_log_qval <- merged_results_sig %>%
        dplyr::group_by(.data$full_metadata_name) %>%
        dplyr::summarise(mean_value = mean(log(.data$qval_joint), na.rm = TRUE))
      
      heatmap_vars <- mean_log_qval$full_metadata_name[order(mean_log_qval$mean_value)]
      heatmap_vars <- setdiff(heatmap_vars, coef_plot_vars)
    }
    
    if (length(coef_plot_vars) > 0 & 
        sum(merged_results_sig$full_metadata_name %in% coef_plot_vars) >= 1) {
      coef_plot_data <- merged_results_sig[merged_results_sig$full_metadata_name %in% coef_plot_vars,]
      
      quantile_df <- coef_plot_data %>%
        dplyr::group_by(.data$full_metadata_name) %>%
        dplyr::summarise(lower_q = median(.data$coef) - 10 * (median(.data$coef) - quantile(.data$coef, 0.25)), 
                         upper_q = median(.data$coef) + 10 * (quantile(.data$coef, 0.75) - median(.data$coef))) %>%
        data.frame()
      rownames(quantile_df) <- quantile_df$full_metadata_name
      
      # Make sure insignificant coefficients don't distort the plot
      coef_plot_data <- coef_plot_data[coef_plot_data$qval_individual < max_significance | 
                                         (coef_plot_data$coef > quantile_df[coef_plot_data$full_metadata_name, 'lower_q'] & 
                                            coef_plot_data$coef < quantile_df[coef_plot_data$full_metadata_name, 'upper_q']),]
      
      custom_break_fun <- function(n) {
        return(function(x) {
          extended_breaks <- scales::breaks_extended(n)(x)
          if (max(x) > 0) {
            extended_breaks <- extended_breaks[extended_breaks <= max(x) * 0.9]
          } else {
            extended_breaks <- extended_breaks[extended_breaks <= max(x) * 1.1]
          }
          if (min(x) > 0) {
            extended_breaks <- extended_breaks[extended_breaks >= min(x) * 1.1]
          } else {
            extended_breaks <- extended_breaks[extended_breaks >= min(x) * 0.9]
          }
          extended_breaks
        })
      }
      
      # Create plot
      p1 <- ggplot2::ggplot(coef_plot_data, ggplot2::aes(x=.data$coef, y=.data$feature))
      
      if (median_comparison_prevalence | median_comparison_abundance) {
        p1 <- p1 + 
          ggplot2::guides(
            linetype = ggplot2::guide_legend(title = 'Null hypothesis', order = 1),
          ) + 
          ggplot2::geom_vline(data = median_df[median_df$full_metadata_name %in% coef_plot_vars,], 
                              ggplot2::aes(xintercept = .data$median_val, linetype = .data$model), color = "darkgray") +
          ggplot2::scale_linetype_manual(values = c("Prevalence" = "dashed", "Abundance" = "solid"))
      } else {
        p1 <- p1 + 
          ggplot2::geom_vline(ggplot2::aes(xintercept = 0), color = "darkgray", linetype = 'dashed')
      }
      
      scale_fill_gradient_limits <- c(min(max_significance, 10^floor(log10(min(coef_plot_data$qval_individual)))), 1)
      if (min(coef_plot_data$qval_individual) < max_significance) {
        scale_fill_gradient_breaks <- c(10^floor(log10(min(coef_plot_data$qval_individual))), max_significance, 1)
      } else {
        scale_fill_gradient_breaks <- c(max_significance, 1)
      }
      if (min(coef_plot_data$qval_individual) < max_significance) {
        scale_fill_gradient_labels <- c(paste0("1e", floor(log10(min(coef_plot_data$qval_individual)))), 
                                        paste0(max_significance),
                                        "1")
      } else {
        scale_fill_gradient_labels <- c(paste0(max_significance),
                                        "1")
      }
      
      
      p1 <- p1 +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$coef - .data$stderr, xmax = .data$coef + .data$stderr), width = 0.2) + 
        ggplot2::geom_point(data = coef_plot_data[coef_plot_data$model == 'Prevalence',], 
                            ggplot2::aes(shape = .data$model, fill = .data$qval_individual), size = 4.5, color = "black")+
        ggplot2::scale_fill_gradient(low="darkgreen", high="white",
                              limits = scale_fill_gradient_limits,
                              breaks = scale_fill_gradient_breaks,
                              labels = scale_fill_gradient_labels,
                              transform = scales::pseudo_log_trans(sigma = 0.001),
                              name = bquote("Prevalence" ~ P["FDR"])) +
        ggnewscale::new_scale_fill() + 
        ggplot2::geom_point(data = coef_plot_data[coef_plot_data$model == 'Abundance',], 
                            ggplot2::aes(shape = .data$model, fill = .data$qval_individual), size = 4.5, color = "black")+
        ggplot2::scale_fill_gradient(low="purple4", high="white",
                                     limits = scale_fill_gradient_limits,
                                     breaks = scale_fill_gradient_breaks,
                                     labels = scale_fill_gradient_labels,
                            transform = scales::pseudo_log_trans(sigma = 0.001),
                            name = bquote("Abundance" ~ P["FDR"])) +
        ggplot2::scale_x_continuous(breaks = custom_break_fun(n = 6), 
                           limits = c(min(coef_plot_data$coef) - quantile(coef_plot_data$stderr, 0.8), 
                                      max(coef_plot_data$coef) + quantile(coef_plot_data$stderr, 0.8))) +
        ggplot2::scale_shape_manual(name = "Association", values=c(21, 24))+
        ggplot2::guides(
          shape = ggplot2::guide_legend(order = 2),
        ) + 
        ggplot2::labs(x =expression(paste(beta, " coefficient")),  y = "Feature") +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.title = ggplot2::element_text(size = 16),
              axis.text.y = ggplot2::element_text(size = 14),
              axis.text.x =ggplot2::element_text(size = 14),
              legend.title = ggplot2::element_text(size = 16),
              legend.text = ggplot2::element_text(size = 14, face = "plain"),
              legend.position = "right",
              legend.background = ggplot2::element_rect(fill = "transparent"),
              panel.spacing=ggplot2::unit(0, "lines"),
              panel.grid.minor = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(size=14),
              strip.background = ggplot2::element_rect(fill = "transparent")) + 
        ggplot2::facet_wrap(~ factor(full_metadata_name, levels = unique(coef_plot_vars)), scales = 'free_x', ncol = length(coef_plot_vars))
      
    } else {
      p1 <- NULL
    }
    
    # Create column for significance star annotation
    merged_results_sig$sig_star <- cut(merged_results_sig$qval_individual, breaks=c(-Inf, max_significance / 10, max_significance, Inf), label=c("**", "*", ""))  

    # Bin coefficients into categories
    coefficient_thresh <- round(max(abs(quantile(merged_results_sig$coef, c(0.1, 0.9)))) / 10, 1) * 5
    coef_breaks <- c(-coefficient_thresh, -coefficient_thresh / 2, 0, coefficient_thresh / 2, coefficient_thresh, Inf)
    threshold_set <- c(paste0("(-Inf,", -1 * coefficient_thresh, "]"),
                       paste0("(", -1 * coefficient_thresh, ",", -1/2 * coefficient_thresh,"]"),
                       paste0("(", -1/2 * coefficient_thresh, ",0]"),
                       paste0("(0,", 1/2 * coefficient_thresh,"]"),
                       paste0("(", 1/2 * coefficient_thresh, ",", 1 * coefficient_thresh,"]"),
                       paste0("(", 1 * coefficient_thresh, ",Inf)"))
    
    threshold_indices <- sapply(merged_results_sig$coef, function(value) {
      which(value < coef_breaks)[1]
    })
    
    merged_results_sig <- merged_results_sig %>%
      dplyr::mutate(coef_cat = threshold_set[threshold_indices])
    merged_results_sig$coef_cat <- factor(merged_results_sig$coef_cat, levels = threshold_set)
    
    scale_fill_values <- rev((RColorBrewer::brewer.pal(n = 6, name = "RdBu")))
    names(scale_fill_values) <- threshold_set
    
    if (length(heatmap_vars) > 0 & 
      sum(merged_results_sig$full_metadata_name %in% heatmap_vars) >= 1) {
      heatmap_data <- merged_results_sig[merged_results_sig$full_metadata_name %in% heatmap_vars,]
      
      grid <- expand.grid(
        feature = unique(heatmap_data$feature),
        full_metadata_name = unique(heatmap_data$full_metadata_name),
        model = unique(heatmap_data$model)
      )
      heatmap_data <- merge(grid, heatmap_data, by = c("feature", "full_metadata_name", "model"), all.x = TRUE)
      heatmap_data$coef[is.na(heatmap_data$coef)] <- NA

      p2 <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = factor(.data$full_metadata_name, unique(heatmap_vars)), y = .data$feature)) +
        ggplot2::geom_tile(data = heatmap_data, ggplot2::aes(fill = .data$coef_cat), colour="white", linewidth=0.2) +
        ggplot2::scale_fill_manual(name = "Beta coefficient", na.value="#EEEEEE",
                          values = scale_fill_values) + 
        ggplot2::geom_text(ggplot2::aes(label = .data$sig_star, color = .data$sig_star), size = 6, vjust = 0.75, hjust = 0.5, key_glyph = ggplot2::draw_key_blank) +
        ggplot2::scale_color_manual(name = bquote("Covariates" ~ P["FDR"]),
                           breaks = c("*", "**", ""),
                           values = c("black", "black", "black"),
                           labels = c(paste0("* < ", round(max_significance, 3)), paste0("** < ", round(max_significance / 10, 5)), "")) +
        ggplot2::labs(x ='',  y = "Feature", caption = "") +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.title = ggplot2::element_text(size = 16),
              axis.text.x =ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
              legend.title = ggplot2::element_text(size = 16),
              legend.text = ggplot2::element_text(size = 14, face = "plain"),
              legend.position = "right",
              legend.background = ggplot2::element_rect(fill = "transparent"),
              panel.spacing=ggplot2::unit(0, "lines"),
              panel.grid.minor = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(size=14),
              strip.background = ggplot2::element_rect(fill = "transparent")) + 
        ggplot2::guides(
          fill = ggplot2::guide_legend(order = 1),
          color = ggplot2::guide_legend(order = 2),
        ) + 
        ggplot2::facet_grid(~ model, labeller = ggplot2::labeller(model = c("abundance" = "Abundance", "prevalence" = "Prevalence")))
      
      if (!is.null(p1)) {
        p2 <- p2 + ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
        )             
      }
      
    } else {
      p2 <- NULL
    }
    
    if (!is.null(p1) & !is.null(p2)) {
      final_plot <- patchwork::wrap_plots(p1, p2, ncol = 3, 
                                          widths = c(max(0, length(coef_plot_vars) * (max(15, max(nchar(as.character(coef_plot_vars))))) / 15 - 2) + 2,
                                                     max(0, length(heatmap_vars) / 4 - 2) + 2, 
                                                     0.5), guides = 'collect')
    } else if (is.null(p1) & !is.null(p2)) {
      final_plot <- p2
    } else if (!is.null(p1) & is.null(p2)) {
      final_plot <- p1
    } else {
      final_plot <- NULL
    }
    
    if (!is.null(final_plot)) {
      height_out <- 9.5 + max(first_n / 5 - 5, 0) + max(nchar(c(as.character(coef_plot_vars), as.character(heatmap_vars)))) / 10
      width_out <-  5 + max(nchar(merged_results$feature)) / 12 + 
        (length(coef_plot_vars) * (max(20, max(nchar(as.character(coef_plot_vars))))) / 20) * 2.5 + 
        length(heatmap_vars) * 0.25
      
      ggplot2::ggsave(summary_plot_file, plot = final_plot, height = height_out, width = width_out)
      png_file <- file.path(figures_folder, "summary_plot.png")
      ggplot2::ggsave(png_file, plot = final_plot, height = height_out, width = width_out)
    }
}

maaslin3_association_plots <-
    function(
      merged_results,
      metadata,
      features,
      max_significance = 0.1,
      figures_folder,
      max_pngs = 10,
      normalization,
      transform,
      feature_specific_covariate = NULL,
      feature_specific_covariate_name = NULL,
      feature_specific_covariate_record = NULL) {
      
      new_name_normalization <- c('Total sum scaling', 'Center log ratio', 'Cumulative sum scaling', 'None', 'Trimmed means of M values')
      names(new_name_normalization) <- c("TSS", "CLR", "CSS", "NONE", "TMM")
      normalization <- new_name_normalization[normalization]
      
      new_name_transformation <- c('Log base 2', 'Logit', 'Arcsin square root', 'None')
      names(new_name_transformation) <- c("LOG", "LOGIT", "AST", "NONE")
      transformation <- new_name_transformation[transform]
      
      merged_results <- merged_results[is.na(merged_results$error) & 
                                         !is.na(merged_results$qval_individual) & 
                                         merged_results$qval_individual < max_significance,]
      if (nrow(merged_results) == 0) {
        logging::loginfo(
          paste("All associations had errors or were insignificant."))
        return(NULL)
      }
      
      merged_results <- merged_results[order(merged_results$qval_individual),]
      
      logging::loginfo(
          paste("Plotting associations from most",
              "to least significant,",
              "grouped by metadata"))
      
      saved_plots <- list()
      
      features_by_metadata <- unique(merged_results[,c('feature', 'metadata', 'model')])
      
      for (row_num in 1:min(nrow(features_by_metadata), max_pngs)) {
        feature_name <- features_by_metadata[row_num, 'feature']
        feature_abun <- data.frame(sample = rownames(features),
                                   feature_abun = features[,feature_name])
        
        metadata_name <- features_by_metadata[row_num, 'metadata']
        if (!is.null(feature_specific_covariate_name)) {
          if (metadata_name == feature_specific_covariate_name) {
            metadata_sub <- data.frame(sample = rownames(feature_specific_covariate),
                                       metadata = feature_specific_covariate[,feature_name])
          } else {
            metadata_sub <- data.frame(sample = rownames(metadata),
                                       metadata = metadata[,metadata_name])
          }
        } else {
          metadata_sub <- data.frame(sample = rownames(metadata),
                                     metadata = metadata[,metadata_name])
        }
        joined_features_metadata <- dplyr::inner_join(feature_abun, metadata_sub, by = c('sample'))
        
        model_name = features_by_metadata[row_num, 'model']
        
        this_signif_association <- merged_results[merged_results$feature == feature_name & 
                                                    merged_results$metadata == metadata_name &
                                                    merged_results$model == model_name,]
        
        if ('LM' == model_name) {
          coef_val <- this_signif_association[this_signif_association$model == 'LM',]$coef
          qval <- this_signif_association[this_signif_association$model == 'LM',]$qval_individual
          N_nonzero <- this_signif_association[this_signif_association$model == 'LM',]$N.not.zero
          N_total <- this_signif_association[this_signif_association$model == 'LM',]$N
          results_value <- this_signif_association[this_signif_association$model == 'LM',]$value
          
          joined_features_metadata_abun <- joined_features_metadata[!is.na(joined_features_metadata$feature_abun),]
          if (is.numeric(joined_features_metadata_abun$metadata) & 
              length(unique(joined_features_metadata_abun$metadata)) > 1) {
            logging::loginfo(
              "Creating scatter plot for continuous data (linear), %s vs %s",
              metadata_name,
              feature_name)
            temp_plot <- ggplot2::ggplot(
              data = joined_features_metadata_abun,
              ggplot2::aes(
                as.numeric(.data$metadata),
                as.numeric(.data$feature_abun)
              )) +
              ggplot2::geom_point(
                fill = 'darkolivegreen4',
                color = 'black',
                alpha = .5,
                shape = 21,
                size = 1,
                stroke = 0.15
              ) + 
              ggplot2::scale_x_continuous(
                limits = c(min(joined_features_metadata_abun['metadata']), 
                           max(joined_features_metadata_abun['metadata']))) +
              ggplot2::scale_y_continuous(
                limits = c(min(joined_features_metadata_abun['feature_abun']), 
                           max(joined_features_metadata_abun['feature_abun'])),
                expand = ggplot2::expansion(mult = c(0, 0.2))) +
              ggplot2::stat_smooth(
                method = "glm",
                formula = 'y ~ x',
                linewidth = 0.5,
                color = 'blue',
                na.rm = TRUE
              ) + 
              ggplot2::guides(alpha = 'none') + 
              ggplot2::labs("") +
              ggplot2::xlab(metadata_name) + 
              ggplot2::ylab(paste0(feature_name, '\n(Normalization: ', normalization, ', Transformation: ', transformation, ')')) + 
              nature_theme(joined_features_metadata_abun['metadata'], 
                           paste0(feature_name, '\n(Normalization: ', normalization, ', Transformation: ', transformation, ')')) + 
            ggplot2::annotate(
              geom = "text",
              x = Inf,
              y = Inf,
              hjust = 1,
              vjust = 1,
              label = sprintf(
                "FDR: %s\nCoefficient (in full model): %s\nN: %s\nN (not zero): %s",
                formatC(qval, format = "e", digits = 1),
                formatC(coef_val, format = "e", digits = 1),
                formatC(N_total, format = 'f', digits = 0),
                formatC(N_nonzero, format = 'f', digits = 0)
              ) ,
              color = "black",
              size = 2,
              fontface = "italic"
            )
          } else {
            x_axis_label_names <- unique(joined_features_metadata_abun$metadata)
            
            sorted_fixed_order <- order(match(results_value, levels(x_axis_label_names)))
            coef_val <- coef_val[sorted_fixed_order]
            qval <- qval[sorted_fixed_order]
            N_nonzero <- N_nonzero[sorted_fixed_order]
            N_total <- N_total[sorted_fixed_order]
            results_value <- results_value[sorted_fixed_order]
            
            if (!is.null(feature_specific_covariate_name)) {
              if (metadata_name == feature_specific_covariate_name) {
                renamed_levels <- as.character(levels(feature_specific_covariate[,feature_name]))
              } else {
                renamed_levels <- as.character(levels(metadata[,metadata_name]))
              }
            } else {
              renamed_levels <- as.character(levels(metadata[,metadata_name]))
            }
            
            if (length(renamed_levels) == 0) {
              renamed_levels <- x_axis_label_names
            }
            for (name in x_axis_label_names) {
              total <- length(which(joined_features_metadata_abun$metadata == name))
              new_n <- paste(name, " (n=", total, ")", sep="")
              levels(joined_features_metadata_abun[,'metadata'])[
                levels(joined_features_metadata_abun[,'metadata']) == name] <- new_n
              renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n)
            }

            logging::loginfo(
              "Creating box plot for categorical data (linear), %s vs %s",
              metadata_name,
              feature_name)
            
            temp_plot <-
              ggplot2::ggplot(
                data = joined_features_metadata_abun, ggplot2::aes(.data$metadata, .data$feature_abun)) +
              ggplot2::geom_boxplot(
                ggplot2::aes(fill = .data$metadata),
                outlier.alpha = 0.0,
                na.rm = TRUE,
                alpha = .5,
                show.legend = FALSE
              ) +
              ggplot2::geom_point(
                ggplot2::aes(fill = .data$metadata),
                alpha = 0.75 ,
                size = 1,
                shape = 21,
                stroke = 0.15,
                color = 'black',
                position = ggplot2::position_jitterdodge()
              ) +
              ggplot2::scale_fill_brewer(palette = "Spectral") + 
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2)))
            
            temp_plot <- temp_plot + 
              nature_theme(.data$metadata, 
                           paste0(feature_name, '\n(Normalization: ', normalization, ', Transformation: ', transformation, ')')) +
              ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black")
              ) +
              ggplot2::xlab(metadata_name) +
              ggplot2::ylab(paste0(feature_name, '\n(Normalization: ', normalization, ', Transformation: ', transformation, ')')) +
              ggplot2::theme(legend.position = "none") +
              ggplot2::annotate(
                geom = "text",
                x = Inf,
                y = Inf,
                hjust = 1,
                vjust = 1,
                label = sprintf(
                  "Value: %s\nFDR: %s\nCoefficient (in full model): %s",
                  paste0(results_value, collapse = ', '),
                  paste0(formatC(qval, format = "e", digits = 1), collapse = ', '),
                  paste0(formatC(coef_val, format = "e", digits = 1), collapse = ', ')
                ) ,
                color = "black",
                size = 2,
                fontface = "italic"
              )
          }
        }
        
        if ('logistic' == model_name) {
          coef_val <- this_signif_association[this_signif_association$model == 'logistic',]$coef
          qval <- this_signif_association[this_signif_association$model == 'logistic',]$qval_individual
          N_nonzero <- this_signif_association[this_signif_association$model == 'logistic',]$N.not.zero
          N_total <- this_signif_association[this_signif_association$model == 'logistic',]$N
          results_value <- this_signif_association[this_signif_association$model == 'logistic',]$value
          
          joined_features_metadata_prev <- joined_features_metadata
          joined_features_metadata_prev$feature_abun <- 
            ifelse(is.na(joined_features_metadata_prev$feature_abun), 'Absent', 'Present')
          
          joined_features_metadata_prev$feature_abun <- 
            factor(joined_features_metadata_prev$feature_abun, levels = c('Present', 'Absent'))
          
          if (is.numeric(joined_features_metadata_prev$metadata) & 
              length(unique(joined_features_metadata_prev$metadata)) > 1) {
            logging::loginfo(
              "Creating boxplot for continuous data (logistic), %s vs %s",
              metadata_name,
              feature_name)
            
            temp_plot <-
              ggplot2::ggplot(
                data = joined_features_metadata_prev, ggplot2::aes(.data$feature_abun, .data$metadata)) +
              ggplot2::geom_boxplot(
                ggplot2::aes(fill = .data$feature_abun),
                outlier.alpha = 0.0,
                na.rm = TRUE,
                alpha = .5,
                show.legend = FALSE
              ) +
              ggplot2::geom_point(
                ggplot2::aes(fill = .data$feature_abun),
                alpha = 0.75 ,
                size = 1,
                shape = 21,
                stroke = 0.15,
                color = 'black',
                position = ggplot2::position_jitterdodge()
              ) +
              ggplot2::scale_fill_brewer(palette = "Spectral") + 
              ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0.7)))
            
            temp_plot <- temp_plot + 
              nature_theme(metadata_name, joined_features_metadata_prev['feature_abun']) +
              ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black")
              ) +
              ggplot2::xlab(feature_name) +
              ggplot2::ylab(metadata_name) +
              ggplot2::theme(legend.position = "none") +
              ggplot2::annotate(
                geom = "text",
                x = Inf,
                y = Inf,
                hjust = 1,
                vjust = 1,
                label = sprintf(
                  "FDR: %s\nCoefficient (in full model): %s\nN: %s\nN (not zero): %s",
                  formatC(qval, format = "e", digits = 1),
                  formatC(coef_val, format = "e", digits = 1),
                  formatC(N_total, format = 'f', digits = 0),
                  formatC(N_nonzero, format = 'f', digits = 0)
                ) ,
                color = "black",
                size = 2,
                fontface = "italic"
              ) + 
              ggplot2::coord_flip()
            
          } else {
            joined_features_metadata_prev$feature_abun <- 
              factor(joined_features_metadata_prev$feature_abun, levels = c('Absent', 'Present'))
            
            x_axis_label_names <- unique(joined_features_metadata_prev$metadata)
            
            sorted_fixed_order <- order(match(results_value, levels(x_axis_label_names)))
            coef_val <- coef_val[sorted_fixed_order]
            qval <- qval[sorted_fixed_order]
            N_nonzero <- N_nonzero[sorted_fixed_order]
            N_total <- N_total[sorted_fixed_order]
            results_value <- results_value[sorted_fixed_order]
            
            if (!is.null(feature_specific_covariate_name)) {
              if (metadata_name == feature_specific_covariate_name) {
                renamed_levels <- as.character(levels(feature_specific_covariate[,feature_name]))
              } else {
                renamed_levels <- as.character(levels(metadata[,metadata_name]))
              }
            } else {
              renamed_levels <- as.character(levels(metadata[,metadata_name]))
            }
            
            if (length(renamed_levels) == 0) {
              renamed_levels <- x_axis_label_names
            }
            for (name in x_axis_label_names) {
              mean_abun <- mean(joined_features_metadata_prev$feature_abun[
                joined_features_metadata_prev$metadata == name] == 'Present', na.rm=T)
              new_n <- paste(name, " (p = ", round(mean_abun, 2)*100, "%)", sep="")
              levels(joined_features_metadata_prev[,'metadata'])[
                levels(joined_features_metadata_prev[,'metadata']) == name] <- new_n
              renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n)
            }
            
            logging::loginfo(
              "Creating tile plot for categorical data (logistic), %s vs %s",
              metadata_name,
              feature_name)
            
            count_df <- joined_features_metadata_prev %>%
              dplyr::group_by(.data$feature_abun, .data$metadata) %>%
              dplyr::summarise(count = dplyr::n(), .groups = 'drop')
            
            x_vals <- unique(joined_features_metadata_prev$feature_abun)
            y_vals <- unique(joined_features_metadata_prev$metadata)
            complete_grid <- expand.grid(feature_abun = x_vals, 
                                         metadata = y_vals)
            
            table_df <- complete_grid %>%
              dplyr::left_join(count_df, by = c("feature_abun", "metadata")) %>%
              dplyr::mutate(count = ifelse(is.na(.data$count), 0, .data$count))

            temp_plot <- ggplot2::ggplot(table_df, ggplot2::aes(x = .data$metadata, y = .data$feature_abun)) +
              ggplot2::geom_tile(ggplot2::aes(fill = .data$count), color = "white",
                        lwd = 1.5,
                        linetype = 1) +
              ggplot2::geom_text(ggplot2::aes(label = .data$count), color = "black", size = 32 / nrow(table_df)) + 
              ggplot2::coord_fixed(ratio = 0.5)+
              ggplot2::scale_fill_gradient2(low = "#075AFF",
                                   mid = "#FFFFCC",
                                   high = "#FF0000") +
              ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 1 + nrow(table_df) / 6))) + 
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                    legend.position = "none")
              
            temp_plot <- temp_plot + 
              nature_theme(as.character(table_df$metadata), joined_features_metadata_prev['feature_abun']) + 
              ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black")
              ) +
              ggplot2::xlab(metadata_name) +
              ggplot2::ylab(feature_name) +
              ggplot2::theme(legend.position = "none") +
              ggplot2::annotate(
                geom = "text",
                x = Inf,
                y = Inf,
                hjust = 1,
                vjust = 1,
                label = sprintf(
                  "Value: %s\nFDR: %s\nCoefficient (in full model): %s",
                  paste0(results_value, collapse = ', '),
                  paste0(formatC(qval, format = "e", digits = 1), collapse = ', '),
                  paste0(formatC(coef_val, format = "e", digits = 1), collapse = ', ')
                ) ,
                color = "black",
                size = 2,
                fontface = "italic"
              )
          }
        }
        
        saved_plots[[metadata_name]][[feature_name]][[model_name]] <- temp_plot
      }
      
      association_plots_folder <- file.path(figures_folder, 'association_plots')
      if (!file.exists(association_plots_folder)) {
        dir.create(association_plots_folder)
      }
      
      for (metadata_variable in names(saved_plots)) {
        for (feature in names(saved_plots[[metadata_variable]])) {
          for (model_name in names(saved_plots[[metadata_variable]][[feature]])) {
            this_plot <- saved_plots[[metadata_variable]][[feature]][[model_name]]

            png_file <- file.path(association_plots_folder,
                                  paste0(metadata_variable, '_', feature, "_", model_name, ".png"))
            height <- max(960, 15 * max(nchar(unlist(strsplit(this_plot$labels$y, '\n')))))
            ggplot2::ggsave(filename = png_file, plot = this_plot, dpi = 300, width = 960/300, height = height/300)
          }
        }
      }
      
      return(saved_plots)
}
