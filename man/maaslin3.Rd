\name{maaslin3}
\alias{maaslin3}
\title{
MaAsLin 3 is the next generation of MaAsLin, a multivariable
    statistical framework for finding abundance and prevalence 
    associations between clinical metadata and high-dimensional 
    microbial multi-omics data.
}
\description{
This wrapper for all MaAsLin 3 steps finds abundance and prevalence 
    associations between microbiome meta-omics features and 
    complex metadata in population-scale epidemiological studies. 
    The software includes multiple analysis methods (including 
    support for multiple covariates, repeated measures, and ordered predictors),
    filtering, normalization, and transform options to customize analysis for your specific study.
}
\usage{
maaslin3(param_list)
}
\arguments{
    \item{param_list}{
    The parameter list containing the following required named items:
    \describe{
      \item{(1)}{\code{input_data}: A data frame of feature abundances or read counts or a filepath to a tab-delimited file with abundances. It should be formatted with features as columns and samples as rows (or the transpose). The column and row names should be the feature names and sample names respectively.}
      \item{(2)}{\code{input_metadata}: A data frame of per-sample metadata or a filepath to a tab-delimited file with metadata. It should be formatted with variables as columns and samples as rows (or the transpose). The column and row names should be the variable names and sample names respectively.}
      \item{(3)}{\code{output}: The output folder to write results.}
    }
}
}
\value{
    List containing the results from applying the model.
}
\author{
    Himel Mallick<himel.stat.iitk@gmail.com>,\cr
    Ali Rahnavard<gholamali.rahnavard@gmail.com>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{

}
