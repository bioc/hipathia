##
## data.R
## Data documentation script
##
## Written by Marta R. Hidalgo
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##



#' BRCA gene expression dataset as SummarizedExperiment
#'
#' A dataset containing a matrix with the Gene expression of 40 samples from the
#' BRCA-US project from The Cancer Genome Atlas (TCGA), and their experimental
#' design, containing 20 "Tumor" samples 20 "Normal" samples.
#'
#' The gene expression matrix includes 40 samples. The data has been
#' log-transformed and normalized with TMM.
#'
#' @format SummarizedExperiment. The assay is a matrix with 40 columns and
#' 18638 rows. Row names are Entrez IDs and column names are the TCGA
#' identifyers of the samples. The colData() is a data.frame with 1 column and
#' 40 rows, including the experimental design of the 40 samples from the BRCA-US
#' project from TCGA. Field \code{group} is the type of sample, either "Tumor"
#' or "Normal".
#'
#' @return SummarizedExperiment including a matrix with 40 columns and 18638
#' rows. Row names are Entrez IDs and column names are the TCGA identifyers of
#' the samples.
#'
#' @usage data(brca)
#'
#' @source \url{https://cancergenome.nih.gov/}
#'
"brca"



#' Results object
#'
#' Results object returned by \code{hipathia::hipathia} function, after calling
#' \code{hidata <- hipathia(brca, pathways, verbose=TRUE, uni.terms = TRUE,
#' GO.terms = TRUE)}
#'
#' @format MultiAssayExperiment object of 4 listed experiments, with the
#' activity values of nodes, paths and functional annotations for each sample:
#' Nodes includes a matrix with 6826 rows
#' Paths includes a matrix with 1876 rows
#' Uni.terms includes a matrix with 142 rows
#' GO.terms includes a matrix with 1654 rows
#'
#' @return Object of results, including nodes, pathways and functional
#' information.
#'
#' @usage data(hidata)
#'
"hidata"


#' Wilcoxon and limma comparison object for nodes, pathways and functional
#' annotations
#'
#' Comparison object returned by \code{hipathia::DAcomp} function, after
#' calling
#' \code{DAdata <- DAcomp(hidata, "group", g1 = "Tumor", g2 = "Normal")}
#'
#' @format List object with 4 entries:
#' Nodes includes a matrix with 6826 rows and 8 columns
#' Paths includes a matrix with 1876 rows and 13 columns
#' Uni.terms includes a matrix with 142 rows and 6 columns
#' GO.terms includes a matrix with 1654 rows and 6 columns
#'
#' @return List of tibbles with the comparison results
#'
#' @usage data(DAdata)
#'
"DAdata"

#' Pathways object including pathways has03320 and hsa04012.
#'
#' Pathways object returned by \code{hipathia::load_pathways} function, after
#' calling
#' \code{pathways <- load_pathways(species = "hsa",
#' pathways_list = c("hsa03320", "hsa04012"))}
#'
#' @format Pathways object
#'
#' @return Pathways object including pathways has03320 and hsa04012.
#'
#' @usage data(pathways)
#'
"pathways"

#' BRCA gene expression dataset
#'
#' Gene expression of 40 samples from the BRCA-US project from
#' The Cancer Genome Atlas (TCGA).
#'
#' Gene expression matrix with 40 samples taken from the BRCA-US project from
#' The Cancer Genome Atlas (TCGA). The data has been log-transformed and
#' normalized with TMM.
#'
#' @format Matrix with 40 columns and 18638 rows. Row names are Entrez IDs
#' and column names are the  TCGA identifyers of the samples.
#'
#' @return Matrix with 40 columns and 18638 rows. Row names are Entrez IDs
#' and column names are the  TCGA identifyers of the samples.
#'
#' @usage data(brca_data)
#'
#' @source \url{https://cancergenome.nih.gov/}
#'
"brca_data"


#' BRCA experimental design
#'
#' Experimental design of the gene expression matrix \code{brca_data} with
#' 40 samples taken from the BRCA-US project from The Cancer Genome Atlas
#' (TCGA). 20 samples are "Tumor" samples and 20 samples are "Normal" samples.
#'
#' @format Dataframe with 1 column and 40 rows, including the experimental
#' design of the 40 samples from the BRCA-US project from TCGA. Field
#' \code{group} is the type of sample, either "Tumor" or "Normal".
#'
#' @return Dataframe with 1 column and 40 rows, including the experimental
#' design of the 40 samples from the BRCA-US project from TCGA. Field
#' \code{group} is the type of sample, either "Tumor" or "Normal".
#'
#' @usage data(brca_design)
#'
#' @source \url{https://cancergenome.nih.gov/}
#'
"brca_design"




#' Normalized BRCA gene expression dataset
#'
#' Experimental design matrix once expression matrix \code{brca_data} has been
#' translated to Entrez geens with \code{translate_matrix} and normalized using
#' \code{normalize_data}.
#'
#' To create the data, the following functions have been called:
#' \code{trans_data <- translate_matrix(brca_data, "hsa")}
#' \code{exp_data <- normalize_data(trans_data)}
#'
#' @format Matrix with 40 columns and 3184 rows. Row names are Entrez IDs
#' and column names are the  TCGA identifyers of the samples.
#'
#' @return Matrix with 40 columns and 3184 rows. Row names are Entrez IDs
#' and column names are the  TCGA identifyers of the samples.
#'
#' @usage data(exp_data)
"exp_data"




#' Results object
#'
#' Results object returned by \code{hipathia::hipathia} function, after calling
#' \code{results <- hipathia(exp_data, pathways, verbose=TRUE)}
#'
#' @format Object of results, including pathways information.
#'
#' @return Object of results, including pathways information.
#'
#' @usage data(results)
#'
"results"


#' Wilcoxon comparison of pathways object
#'
#' Comparison object returned by \code{hipathia::do_wilcoxon} function, after
#' calling
#' \code{comp <- do_wilcoxon(path_vals, sample_group, g1 = "Tumor", g2 =
#' "Normal")}
#' \code{path_names <- get_path_names(pathways, rownames(comp))}
#' \code{comp <- cbind(path_names, comp)}
#'
#' @format Table with 1868 rows and 5 columns
#'
#' @return Pathway comparison result
#'
#' @usage data(comp)
#'
"comp"

#' Pathways matrix of the BRCA gene expression dataset
#'
#' Matrix of pathway activation values for the BRCA dataset. This matrix is
#' extracted from the Results object returned by the \code{hipathia} function
#' by means of the \code{get_paths_matrix} function.
#'
#' \code{path_vals <- get_paths_matrix(results)}
#'
#' @format Matrix with 40 columns and 1868 rows. Row names are Pathway IDs
#' and column names are the TCGA identifyers of the samples.
#'
#' @return Matrix with 40 columns and 1868 rows. Row names are Pathway IDs
#' and column names are the TCGA identifyers of the samples.
#'
#' @usage data(path_vals)
#'
"path_vals"


#' Gene Ontology matrix of the BRCA gene expression dataset
#'
#' Matrix of Gene Ontology terms activation values for the BRCA dataset.
#' This matrix is computed from the Results object returned by the
#' \code{hipathia} function by means of the \code{quantify_terms} function.
#'
#' \code{go_vals <- quantify_terms(results, pathways, "GO")}
#'
#' @format Matrix with 40 columns and 1654 rows. Row names are Gene Ontology
#' terms and column names are the TCGA identifyers of the samples.
#'
#' @return  Matrix with 40 columns and 1654 rows. Row names are Gene Ontology
#' terms and column names are the TCGA identifyers of the samples.
#'
#' @usage data(go_vals)
#'
"go_vals"



##' Uniprot matrix of the BRCA gene expression dataset
##'
##' Matrix of Uniprot functions activation values for the BRCA dataset.
##' This matrix is computed from the Results object returned by the
##' \code{hipathia} function by means of the \code{quantify_terms} function.
##'
##' \code{uni_vals <- quantify_terms(results, pathways, "uniprot")}
##'
##' @format Matrix with 40 columns and 142 rows. Row names are Uniprot functions
##' and column names are the TCGA identifyers of the samples.
##'
##' @return  Matrix with 40 columns and 142 rows. Row names are Uniprot
##' functions and column names are the TCGA identifyers of the samples.
##'
##' @usage data(uni_vals)
##'
#"uni_vals"


