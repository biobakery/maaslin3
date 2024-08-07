\name{maaslin_contrast_test}
\alias{maaslin_contrast_test}
\title{
Perform a contrast from a fit MaAsLin 3 model.
}
\description{
Perform a contrast test (\code{lmerTest::contest} for mixed effects linear; 
\code{multcomp::glht} for all others) using a named contrast matrix and 
right hand side. One contrast test is applied per row of the matrix.
}
\usage{
maaslin_contrast_test(fits,
                    contrast_mat,
                    rhs)
}
\arguments{
    \item{fits}{The list of fits from \code{\link{maaslin_fit}} with 
    \code{save_models = TRUE}.}
    \item{contrast_mat}{A matrix with one row per contrast test to run. The
    columns will be matched to the coefficients of the model by name. 
    Contrast vector coefficients need not be specified if they would be zero. 
    If row names are provided, they will be used to label the test in the
    results.}
    \item{rhs}{The right hand size of the contrast test. The length should
    be the same as the number of rows in the \code{contrast_mat}.}
}
\value{
A dataframe with the following columns:
\describe{
    \item{(1)}{\code{feature}: The feature involved in the test.}
    \item{(2)}{\code{test}: The name of the \code{contrast_mat} row for the
    test or the row number of the \code{contrast_mat} if no row names are
    specified.}
    \item{(3)}{\code{coef}: The coefficient used in the hypothesis test: this
    is the dot product of the contrast vector and the model coefficients.}
    \item{(4)}{\code{rhs}: The right hand side against which the coefficient
    is test.}
    \item{(5)}{\code{stderr}: The standard error of the coefficient.}
    \item{(7)}{\code{pval_individual}: The (uncorrected) p-value of the
    test.}
    \item{(8)}{\code{error}: Any error produced during the testing.
    NA otherwise.}
    \item{(10)}{\code{model}: \code{abundance} for the abundance models and
    \code{prevalence} for the prevalence models.}
    }
}
\author{
    William Nickols<willnickols@g.harvard.edu>,\cr
    Jacob Nearing<nearing@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{
    # Read features table
    taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package =
    "maaslin3")
    taxa_table <- read.csv(taxa_table_name, sep = '\t')
    
    # Read metadata table
    metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package =
    "maaslin3")
    metadata <- read.csv(metadata_name, sep = '\t')
    
    metadata$diagnosis <-
        factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
    metadata$dysbiosis_state <-
        factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC',
        'dysbiosis_CD'))
    metadata$antibiotics <-
        factor(metadata$antibiotics, levels = c('No', 'Yes'))
    
    #Run MaAsLin3
    fit_out <- maaslin3::maaslin3(input_data = taxa_table,
                                input_metadata = metadata,
                                output = 'output',
                                formula = '~ diagnosis + dysbiosis_state +
                                antibiotics + age + reads',
                                plot_summary_plot = FALSE,
                                plot_associations = FALSE,
                                save_models = TRUE)
                                
    contrast_mat <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 
        ncol = 4, nrow = 2, byrow = TRUE)
        
    colnames(contrast_mat) <- c("diagnosisUC",
                                "diagnosisCD",
                                "dysbiosis_statedysbiosis_UC",
                                "dysbiosis_statedysbiosis_CD")
                                
    rownames(contrast_mat) <- c("diagnosis_test", "dysbiosis_test")
                                
    maaslin_contrast_test(fits = fit_out$fit_data_abundance$fits,
                            contrast_mat = contrast_mat,
                            rhs = c(0, 0))
                                
    unlink('output', recursive=TRUE)
    logging::logReset()
}
