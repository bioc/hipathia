% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chart.R
\name{multiple_pca_plot}
\alias{multiple_pca_plot}
\title{Plots multiple components of a PCA}
\usage{
multiple_pca_plot(
  fit,
  group = NULL,
  sample_colors = NULL,
  comps = seq_len(3),
  plot_variance = FALSE,
  legend = TRUE,
  cex = 2,
  pch = 20,
  main = "Multiple PCA plot",
  save_png = NULL
)
}
\arguments{
\item{fit}{princomp object as returned by \code{do_pca}}

\item{group}{Vector with the group to which each sample belongs.
The samples must be ordered as in \code{path_vals}.
By default, all samples will be assigned to the same class.}

\item{sample_colors}{Named character vector of colors. The names of the
colors must be the classes in \code{group}. Each sample will be
assigned the color corresponding to its class, taken from the
\code{group} vector. By default a color will be assigned
automatically to each class.}

\item{comps}{Vector with the components to be plot}

\item{plot_variance}{Logical, whether to plot the cumulative variance.}

\item{legend}{Boolean, whether to plot a legend in the plot.
Default is TRUE.}

\item{cex}{Graphical parameter from \code{par()} function.}

\item{pch}{Graphical parameter from \code{par()} function.}

\item{main}{Main title of the image}

\item{save_png}{Path to the file where the image as PNG will be saved.
By default, the image is not saved.}
}
\value{
Plots multiple components of a PCA
}
\description{
Plots multiple components of a PCA analysis computed with \code{do_pca}
}
\examples{
data(path_vals)
sample_group <- brca_design[colnames(path_vals),"group"]
pca_model <- do_pca(path_vals[seq_len(ncol(path_vals)),])
multiple_pca_plot(pca_model, sample_group, cex = 3, plot_variance = TRUE)

}
