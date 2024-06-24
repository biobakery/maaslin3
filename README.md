# MaAsLin 3 #

[MaAsLin3](http://huttenhower.sph.harvard.edu/maaslin3)  is the next generation of MaAsLin (Microbiome Multivariable Association with Linear Models). This repository contains the MaAsLin 3 code as an early shell of a Bioconductor package.

If you use the MaAsLin2 software, please cite our manuscript: 

William A. Nickols, Jacob T. Nearing, Kelsey N. Thompson, Curtis Huttenhower MaAsLin 3: Refining and extending generalized multivariate linear models for meta-omic association discovery. (In progress)

If you have questions, please direct it to :   
[MaAsLin3 Forum](https://forum.biobakery.org/c/Downstream-analysis-and-statistics/MaAsLin2)    

--------------------------------------------

## Contents ##
* [Description](#description)
* [Requirements](#requirements)
* [Installation](#installation)
* [How to Run](#how-to-run)
    * [Input Files](#input-files)
    * [Output Files](#output-files)
    * [Run a Demo](#run-a-demo)
    * [Options](#options)
* [Visualization](#visualization)
* [Troubleshooting](#troubleshooting)

## Description ##
MaAsLin2 finds associations between microbiome multi-omics features and complex metadata in population-scale epidemiological studies. The software includes multiple analysis methods (including support for multiple covariates and repeated measures), filtering, normalization, and transform options to customize analysis for your specific study. 

## Requirements ##
MaAsLin2 is an R package that can be run on the command line or as an R function.The following packages are required dependencies:
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
## Installation ##

#### Run MaAsLin 3 on HMP2 example data:
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

#### Run MaAsLin 3 with the inferred abundance options:
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

##### Session Info #####

Session info from running the demo in R can be displayed with the following command.

```{r}
sessionInfo()
```