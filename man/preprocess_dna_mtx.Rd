\name{preprocess_dna_mtx}
\alias{preprocess_dna_mtx}
\title{
Pre-process the DNA covariates for metatranscriptomics
}
\description{
Pre-process the DNA covariates for metatranscriptomics by total-sum-scaling DNA abundances per sample and then, for each sample in each feature:
\enumerate{
      \item{Log 2 transforming the DNA abundance if the DNA abundance is >=0}
      \item{Setting the DNA abundance to \code{log2([minimum non-zero relative abundance in the dataset] / 2)} if the corresponding RNA abundance is non-zero but the DNA abundance is zero}
      \item{Setting the DNA abundance to NA if both are zero}
}
When the DNA is present, the RNA data can be modeled as usual in MaAsLin 3 with \code{log2(DNA)} as a covariate. When the DNA is not present, if the RNA is present, we assume the DNA was missed due to finite read depth, so the DNA abundance is imputed with a small pseudo-count. When neither the DNA nor RNA is present, we assume the gene/microbe was not in the sample and therefore no information about the transcription level can be obtained. Setting the DNA covariate to NA has the effect of dropping the sample when fitting the relevant feature's model in MaAsLin 3. Unlike most MaAsLin functions that will infer the samples from the row names and column names, \bold{the \code{rna_table} must be formated as samples (rows) by features (columns)}.
}
\usage{
preprocess_dna_mtx(dna_table, rna_table)
}
\arguments{
\item{dna_table}{The samples (rows) by features (columns) data frame of DNA or taxon abundances to preprocess. These can be relative abundances or counts.}
\item{rna_table}{The samples (rows) by features (columns) data frame of RNA to preprocess. These can be relative abundances or counts.}
}
\value{
A list containing the following named items:
    \enumerate{
      \item{\code{dna_table}: The table of log2 transformed DNA relative abundances with NAs for any feature-sample pairs for which both the DNA and RNA abundances were 0.}
      \item{\code{rna_table}: The table of total sum scaled RNA abundances. These are not log2 transformed.}
    }
}
\author{
    William Nickols<willnickols@g.harvard.edu>,\cr
    Jacob Nearing<nearing@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
