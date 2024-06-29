# MaAsLin 3 #

[MaAsLin3](http://huttenhower.sph.harvard.edu/maaslin3)  is the next generation of MaAsLin (Microbiome Multivariable Association with Linear Models). This repository contains the MaAsLin 3 code as an early shell of a Bioconductor package.

If you use the MaAsLin2 software, please cite our manuscript: 

William A. Nickols, Jacob T. Nearing, Kelsey N. Thompson, Curtis Huttenhower MaAsLin 3: Refining and extending generalized multivariate linear models for meta-omic association discovery. (In progress)

If you have questions, please direct it to:  [MaAsLin3 Forum](https://forum.biobakery.org/c/Downstream-analysis-and-statistics/MaAsLin2)    

--------------------------------------------

## Contents ##
- [MaAsLin 3](#maaslin-3)
  - [Contents](#contents)
  - [Description](#description)
  - [Requirements](#requirements)
  - [Installation](#installation)
      - [Install using Github and devtools](#install-using-github-and-devtools)
  - [How to Run](#how-to-run)
      - [Run MaAsLin 3 on HMP2 example data:](#run-maaslin-3-on-hmp2-example-data)
      - [Run MaAsLin 3 with the inferred abundance options:](#run-maaslin-3-with-the-inferred-abundance-options)
        - [Session Info](#session-info)

## Description ##
MaAsLin3 finds associations between microbiome multi-omics features and complex metadata in population-scale epidemiological studies. The software includes multiple analysis methods (including support for multiple covariates and repeated measures), filtering, normalization, and transform options to customize analysis for your specific study. 

## Requirements ##
MaAsLin3 is an R package that can be run on the command line or as an R function.The following packages are required dependencies:
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
multcomp
```

## Installation ##

#### Install using Github and devtools
```
library("devtools")
install_github("biobakery/MaAsLin3")
library("maaslin3")
```


## How to Run ##
#### Run MaAsLin 3 on HMP2 example data:

Note that we include `reads_filtered` as a fixed effect since variable read depth over the samples is likely to create prevalence effects. Because these data are compositional, setting `median_comparison_abundance = TRUE` is recommended so that the abundance coefficients are tested against the median coefficient. By contrast, we can set `median_comparison_prevalence = FALSE` since we do not the typical bug to have no prevalence association with the included variables.

```
#Read features table 
taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package = "maaslin3")
taxa_table <- read.csv(taxa_table_name, sep = '\t')
rownames(taxa_table) <- taxa_table$ID; taxa_table$ID <- NULL

#Read metadata table
metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")
metadata <- read.csv(metadata_name, sep = '\t')
rownames(metadata) <- metadata$ID; metadata$ID <- NULL

#Prepare parameter lists 
param_list <- list(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = 'output', 
                   normalization = 'TSS', 
                   transform = 'LOG', 
                   formula = '~ diagnosis + dysbiosisUC + dysbiosisCD + antibiotics + age + reads_filtered + (1 | subject)', 
                   save_scatter = FALSE, 
                   save_models = FALSE, 
                   plot_heatmap = FALSE, 
                   plot_scatter = FALSE, 
                   max_significance = 0.1, 
                   augment = TRUE, 
                   median_comparison_abundance = TRUE, 
                   median_comparison_prevalence = FALSE, 
                   cores=1)

#Run MaAsLin3
fit_out <- maaslin3::maaslin3(param_list)
```

The outputs can now be found in the `fit_out` object or in the directory `output`. Particularly for prevalence associations, make sure to check a feature's number of non-zeros since rare microbes can occasionally produce large effect sizes and small p-values in complex models.

#### Run MaAsLin 3 with the inferred abundance options:
```
abundance_file <- system.file("extdata", "synthetic_abundance.tsv", package = "maaslin3")
abundance <- read.csv(abundance_file, sep = '\t')

metadata_file <- system.file("extdata", "synthetic_metadata.tsv", package = "maaslin3")
metadata <- read.csv(metadata_file, sep = '\t')

scaling_factors_file <- system.file("extdata", "scaling_factors.tsv", package = "maaslin3")
scaling_factors <- read.csv(scaling_factors_file, sep = '\t')

#Prepare parameter lists 
param_list <- list(input_data = abundance, 
                   input_metadata = metadata, 
                   output = 'output_absolute/', 
                   normalization = 'TSS', 
                   transform = 'LOG', 
                   fixed_effects = colnames(metadata)[colnames(metadata) != "ID"], 
                   save_scatter = FALSE, 
                   save_models = F, 
                   plot_heatmap = F, 
                   plot_scatter = F, 
                   max_significance = 0.1, 
                   augment = T, 
                   median_comparison_abundance = TRUE, 
                   median_comparison_prevalence = FALSE, 
                   unscaled_abundance = scaling_factors)

#Run MaAsLin3
fit_out <- maaslin3::maaslin3(param_list)
```

The file `scaling_factors.tsv` gives the scaling factors for normalization (how much of the spiked feature there is on the absolute scale).

##### Session Info #####

Session info from running the demo in R can be displayed with the following command.

```{r}
sessionInfo()
```