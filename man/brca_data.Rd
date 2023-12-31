% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{brca_data}
\alias{brca_data}
\title{BRCA gene expression dataset}
\format{
Matrix with 40 columns and 18638 rows. Row names are Entrez IDs
and column names are the  TCGA identifyers of the samples.
}
\source{
\url{https://cancergenome.nih.gov/}
}
\usage{
data(brca_data)
}
\value{
Matrix with 40 columns and 18638 rows. Row names are Entrez IDs
and column names are the  TCGA identifyers of the samples.
}
\description{
Gene expression of 40 samples from the BRCA-US project from
The Cancer Genome Atlas (TCGA).
}
\details{
Gene expression matrix with 40 samples taken from the BRCA-US project from
The Cancer Genome Atlas (TCGA). The data has been log-transformed and
normalized with TMM.
}
\keyword{datasets}
