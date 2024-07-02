#!/usr/bin/env Rscript
###############################################################################
# MaAsLin2 visualizations

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

# MaAsLin2 theme based on Nature journal requirements
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

# MaAsLin2 summary_plot function for overall view of associations
maaslin3_summary_plot <-
  function(
    merged_results,
    summary_plot_file,
    figures_folder,
    first_n = 25,
    max_significance = 0.1,
    coef_plot_vars = NULL,
    heatmap_vars = NULL,
    median_comparison_abundance = FALSE,
    median_comparison_prevalence = FALSE) {
    
    # Preprocessing
    merged_results <- merged_results[is.na(merged_results$error),]
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

      p1 <- ggplot2::ggplot(coef_plot_data, ggplot2::aes(x=.data$coef, y=.data$feature)) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$coef - .data$stderr, xmax = .data$coef + .data$stderr), width = 0.2) + 
        ggplot2::geom_point(ggplot2::aes(shape = .data$model, fill = .data$qval_individual), size = 4.5, color = "black")+
        ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = 5), 
                           limits = c(min(coef_plot_data$coef) - quantile(coef_plot_data$stderr, 0.8), 
                                      max(coef_plot_data$coef) + quantile(coef_plot_data$stderr, 0.8))) +
        ggplot2::scale_shape_manual(name = "Association", values=c(21, 24))+
        viridis::scale_fill_viridis(option = "viridis", 
                           limits=c(10^floor(log10(min(coef_plot_data$qval_individual))), 1), 
                           breaks=c(10^floor(log10(min(coef_plot_data$qval_individual))), max_significance, 1), 
                           labels = c(paste0("1e", floor(log10(min(coef_plot_data$qval_individual)))), 
                                      paste0("1e", floor(log10(max_significance))),
                                      "1"),
                           transform = scales::pseudo_log_trans(sigma = 0.001),
                           name = expression(P["FDR"]), direction = -1) +
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
        ggplot2::facet_wrap(~ full_metadata_name, scales = 'free_x', ncol = length(coef_plot_vars))
      
      if (median_comparison_prevalence | median_comparison_abundance) {
        p1 <- p1 + 
          ggplot2::guides(
            shape = ggplot2::guide_legend(override.aes = list(color = "black")),
            linetype = ggplot2::guide_legend(title = 'Association median'),
          ) + 
          ggplot2::geom_vline(data = median_df[median_df$full_metadata_name %in% coef_plot_vars,], 
                     ggplot2::aes(xintercept = .data$median_val, linetype = .data$model), color = "black") +
          ggplot2::scale_linetype_manual(values = c("Prevalence" = "dashed", "Abundance" = "solid"))
      } else {
        p1 <- p1 + 
          ggplot2::guides(
            shape = ggplot2::guide_legend(override.aes = list(color = "black")),
          ) + 
          ggplot2::geom_vline(ggplot2::aes(xintercept = 0), color = "black", linetype = 'dashed')
      }
      
    } else {
      p1 <- NULL
    }
    
    # Create column for significance star annotation
    merged_results_sig$sig_star <- cut(merged_results_sig$qval_individual, breaks=c(-Inf, max_significance / 10, max_significance, Inf), label=c("**", "*", ""))  

    # Bin coefficients into categories
    coefficient_thresh <- round(max(abs(quantile(merged_results_sig$coef, c(0.1, 0.9)))) / 10, 1) * 5
    coef_breaks <- c(-Inf, -coefficient_thresh, -coefficient_thresh / 2, coefficient_thresh / 2, coefficient_thresh, Inf)
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

      p2 <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$full_metadata_name, y = .data$feature)) +
        ggplot2::geom_tile(data = heatmap_data, ggplot2::aes(fill = .data$coef_cat), colour="white", linewidth=0.2) +
        ggplot2::scale_fill_manual(name = "Beta coefficient", na.value="#EEEEEE",
                          values = scale_fill_values) + 
        ggplot2::geom_text(label = heatmap_data$sig_star, color = "black", size=6, vjust = 0.75, hjust = 0.5) + 
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
      final_plot <- patchwork::wrap_plots(p1, p2, ncol = 3, widths = c(max(0, length(coef_plot_vars) - 2) + 2, 
                                         max(0, length(heatmap_vars) / 4 - 2) + 2, 0.5), guides = 'collect')
    } else if (is.null(p1) & !is.null(p2)) {
      final_plot <- p2
    } else if (!is.null(p1) & is.null(p2)) {
      final_plot <- p1
    } else {
      final_plot <- NULL
    }
    
    return(final_plot)
}

save_summary_plot <-
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

      # generate a heatmap and save it to a pdf and as a png
      summary_plot <-
          maaslin3_summary_plot(
            merged_results,
            summary_plot_file,
            figures_folder,
            first_n,
            max_significance,
            coef_plot_vars,
            heatmap_vars,
            median_comparison_abundance,
            median_comparison_prevalence)

      if (!is.null(summary_plot)) {
        height_out <- 8 + max(first_n / 5 - 5, 0)
        width_out <-  4 + max(nchar(merged_results$feature)) / 10 + 
          ifelse(is.null(coef_plot_vars), 3, length(coef_plot_vars) * 2.5) + 
          ifelse(is.null(heatmap_vars), 1.5, length(heatmap_vars) * 0.2)
        
        ggplot2::ggsave(summary_plot_file, plot = summary_plot, height = height_out, width = width_out)
        png_file <- file.path(figures_folder, "summary_plot.png")
        ggplot2::ggsave(png_file, plot = summary_plot, height = height_out, width = width_out)
      }
}

maaslin3_association_plots <-
    function(
      merged_results,
      metadata,
      features,
      max_significance = 0.1,
      figures_folder,
      max_pngs = 10) {
      
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
      
      features_by_metadata <- unique(merged_results[,c('feature', 'metadata')])
      
      for (row_num in 1:min(nrow(features_by_metadata), max_pngs)) {
        feature_name <- features_by_metadata[row_num, 'feature']
        feature_abun <- data.frame(sample = rownames(features),
                                   feature_abun = features[,feature_name])
        
        metadata_name <- features_by_metadata[row_num, 'metadata']
        metadata_sub <- data.frame(sample = rownames(metadata),
                               metadata = metadata[,metadata_name])
        joined_features_metadata <- dplyr::inner_join(feature_abun, metadata_sub, by = c('sample'))
        
        this_signif_association <- merged_results[merged_results$feature == feature_name & 
                                                    merged_results$metadata == metadata_name,]
        
        if ('LM' %in% this_signif_association$model) {
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
              ggplot2::ylab(feature_name) + 
              nature_theme(joined_features_metadata_abun['metadata'], feature_name) + 
            ggplot2::annotate(
              geom = "text",
              x = Inf,
              y = Inf,
              hjust = 1,
              vjust = 1,
              label = sprintf(
                "FDR: %s\nCoefficient (in full model): %s\nN: %s\nN (not zero): %s",
                formatC(qval, format = "e", digits = 3),
                formatC(coef_val, format = "e", digits = 2),
                formatC(N_total, format = 'f', digits = 0),
                formatC(N_nonzero, format = 'f', digits = 0)
              ) ,
              color = "black",
              size = 2,
              fontface = "italic"
            )
          } else {
            x_axis_label_names <- unique(joined_features_metadata_abun$metadata)
            renamed_levels <- as.character(levels(metadata[,metadata_name]))
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
              nature_theme(metadata_name, joined_features_metadata_abun['feature_abun']) +
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
                  "FDR: %s\nCoefficient (in full model): %s\nValue: %s",
                  paste0(formatC(qval, format = "e", digits = 3), collapse = ', '),
                  paste0(formatC(coef_val, format = "e", digits = 2), collapse = ', '),
                  paste0(results_value, collapse = ', ')
                ) ,
                color = "black",
                size = 2,
                fontface = "italic"
              )
          }
        }
        
        if ('logistic' %in% this_signif_association$model) {
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
                  formatC(qval, format = "e", digits = 3),
                  formatC(coef_val, format = "e", digits = 2),
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
            renamed_levels <- as.character(levels(metadata[,metadata_name]))
            if (length(renamed_levels) == 0) {
              renamed_levels <- x_axis_label_names
            }
            for (name in x_axis_label_names) {
              mean_abun <- mean(joined_features_metadata_prev$feature_abun[
                joined_features_metadata_prev$metadata == name] == 'Present')
              new_n <- paste(name, " (p=", round(mean_abun, 2), ")", sep="")
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
              ggplot2::coord_fixed(ratio = 0.5)+
              ggplot2::geom_text(ggplot2::aes(label = .data$count), color = "black", size = 4) + 
              ggplot2::scale_fill_gradient2(low = "#075AFF",
                                   mid = "#FFFFCC",
                                   high = "#FF0000") +
              ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0, 1.7))) + 
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                    legend.position = "none")
              
            temp_plot <- temp_plot + 
              nature_theme(metadata_name, joined_features_metadata_prev['feature_abun']) + 
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
                  "FDR: %s\nCoefficient (in full model): %s\nValue: %s",
                  paste0(formatC(qval, format = "e", digits = 3), collapse = ', '),
                  paste0(formatC(coef_val, format = "e", digits = 2), collapse = ', '),
                  paste0(results_value, collapse = ', ')
                ) ,
                color = "black",
                size = 2,
                fontface = "italic"
              )
          }
        }
        
        saved_plots[[metadata_name]][[feature_name]] <- temp_plot
      }
      
      association_plots_folder <- file.path(figures_folder, 'association_plots')
      if (!file.exists(association_plots_folder)) {
        dir.create(association_plots_folder)
      }
      
      for (metadata_variable in names(saved_plots)) {
        for (feature in names(saved_plots[[metadata_variable]])) {
          this_plot <- saved_plots[[metadata_variable]][[feature]]
          
          png_file <- file.path(association_plots_folder,
                                paste0(metadata_variable, '_', feature, ".png"))
          height <- max(960, 15 * nchar(feature))
          ggplot2::ggsave(filename = png_file, plot = this_plot, dpi = 300, width = 960/300, height = height/300)
        }
      }
      
      return(saved_plots)
}
