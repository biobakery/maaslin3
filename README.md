# MaAsLin 3 #

MaAsLin 3 is the next generation of MaAsLin (Microbiome Multivariable Association with Linear Models). This repository contains the MaAsLin 3 code as an early shell of a Bioconductor package.

The following packages are required dependencies:
```
optparse
logging
data.table
dplyr
pbapply
lmerTest
parallel
lme4
plyr
TcGSA
ggplot2
grid
pheatmap
gridExtra
```

To install MaAsLin 3, clone this repository and navigate to the cloned folder. To load the Maaslin3 function, run the following after setting `Maaslin3_path` to be the path to `../Maaslin3/R/`:
```
for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  source(file.path(Maaslin3_path, R_file))
}
```

To run MaAsLin 3 on example data from HMP2, run the following from the MaAsLin 3 directory:
```
taxa_table <- read.csv('data/HMP2_taxonomy.tsv', sep = '\t')
rownames(taxa_table) <- taxa_table$ID; taxa_table$ID <- NULL
metadata <- read.csv('data/HMP2_metadata.tsv', sep = '\t')
rownames(metadata) <- metadata$ID; metadata$ID <- NULL

param_list <- list(input_data = taxa_table, input_metadata = metadata, min_abundance = 0, min_prevalence = 0, output = 'tmp/', 
                   min_variance = 0, normalization = 'TSS', transform = 'LOG', analysis_method = 'LM', 
                   formula = '~ diagnosis + dysbiosisUC + dysbiosisCD + antibiotics + age + (1 | subject)', 
                   save_scatter = FALSE, save_models = FALSE, plot_heatmap = T, plot_scatter = F, 
                   max_significance = 0.1, augment = TRUE, iterative_mode = TRUE, cores=1)
fit_out <- Maaslin3(param_list)

write.table(rbind(fit_out$fit_data_non_zero$results, 
                  fit_out$fit_data_binary$results), 
            'results.tsv', sep = '\t', row.names = F)
```

To run MaAsLin 3 with the inferred abundance options, run the following from the MaAsLin 3 directory:
```
abundance <- read.csv('data/synthetic_abundance.tsv', sep = '\t')
metadata <- read.csv('data/synthetic_metadata.tsv', sep = '\t')
scaling_factors <- read.csv('data/scaling_factors.tsv', sep = '\t')

param_list <- list(input_data = abundance, input_metadata = metadata, min_abundance = 0, min_prevalence = 0.0, output = tmp_fit_out, 
                    min_variance = 0, normalization = 'TSS', transform = 'LOG', analysis_method = 'LM', 
                    fixed_effects = colnames(metadata)[colnames(metadata) != "ID"], save_scatter = FALSE, 
                    save_models = F, plot_heatmap = F, plot_scatter = F, max_significance = 0.1, augment = T, iterative_mode = F,
                    unscaled_abundance = scaling_factors)
fit_out <- Maaslin3(param_list)
write.table(rbind(fit_out$fit_data_non_zero$results, 
                  fit_out$fit_data_binary$results), 
            'results_inferred_absolute.tsv', sep = '\t', row.names = F)
```

The file `scaling_factors.tsv` gives the scaling factors for normalization (how much of the spiked feature there is on the absolute scale).
