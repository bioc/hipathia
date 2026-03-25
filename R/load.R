##
## load.R
## Load objects from Zenodo
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


get_hpannot_version <- function(){
    return("v3")
}

get_package_folder <- function(){
  return(find.package("hipathia"))
}

#' Loads a RData object and downloads it from Zenodo if necessary
#'
#' @param file File name of the object to load
#'
#' #@examples
#' #get_file("xref_hsa_v3.rda")
#' #get_file("meta_graph_info_hsa_v3.rda")
#'
#' @return Object in the file
#' @import zen4R
#'
get_file <- function(file){
  pfold <- get_package_folder()
  v <- get_hpannot_version()
  filepath <- file.path(pfold, "data", v)
  if(!file.exists(filepath))
    dir.create(filepath)
  files <- list.files(filepath)
  if(!file %in% files){
    download_zenodo("10.5281/zenodo.18268423", path = filepath,
      files = list(file))
  }
  filename <- load(file.path(filepath, file))
  x <- get(filename)
  return(x)
}


#' Loads annotations object
#'
#' @param db Database to be used. Either "GO" or "uniprot".
#' @param species Species of the samples.
#'
#' #@examples
#' #load_annofuns("GO", "hsa")
#' #load_annofuns("uniprot", "hsa")
#'
#' @return Annotations object
#' @import zen4R
#'
load_annofuns <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    v <- get_hpannot_version()
    file <- paste0("annofuns_", db, "_", species, "_", v, ".rda")
    annofuns <- get_file(file)
    return(annofuns)
}


#' Loads object with graph information
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_mgi("hsa")
#'
#' @return Graph information object
#' @import zen4R
#'
load_mgi <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("meta_graph_info_", species, "_", v, ".rda")
    mgi <- get_file(file)
    return(mgi)
}


#' Loads object with pseudo graph information
#'
#' @param species Species of the samples.
#' @param group_by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include.
#'
#' #@examples
#' #load_pseudo_mgi("hsa", "uniprot")
#'
#' @return Pseudo graph information object
#' @import zen4R
#'
load_pseudo_mgi <- function(species, group_by){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_grouping(group_by))
        stop("Grouping not accepted")
    v <- get_hpannot_version()
    file <- paste0("pmgi_", species, "_", group_by, "_", v, ".rda")
    pmgi <- get_file(file)
    return(pmgi)
}


#' Loads table of references
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_xref("hsa")
#'
#' @return Table of references
#' @import zen4R
#'
load_xref <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("xref_", species, "_", v, ".rda")
    xref <- get_file(file)
    return(xref)
}


#' Loads table of translation from HGNC to Entrez
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_entrez_hgnc("hsa")
#'
#' @return Table of translation from HGNC to Entrez
#' @import zen4R
#'
load_entrez_hgnc <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("entrez_hgnc_", species, "_", v, ".rda")
    entrez_hgnc <- get_file(file)
    return(entrez_hgnc)
}


#' Loads functional annotations to genes
#'
#' Loads functional annotations from HGNC to the selected database.
#'
#' @param db Database to be used. Either "GO" or "uniprot".
#' @param species Species of the samples.
#'
#' #@examples
#' #load_annots("GO", "hsa")
#'
#' @return Functional annotations from HGNC to the selected database.
#' @import zen4R
#'
load_annots <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    v <- get_hpannot_version()
    file <- paste0("annot_", db, "_", species, "_", v, ".rda")
    annot <- get_file(file)
    return(annot)
}


#' Loads GO graph information
#'
#' #@examples
#' #load_gobp_frame()
#'
#' @return GO graph information
#' @import zen4R
#'
load_gobp_frame <- function(){
  v <- get_hpannot_version()
    file <- paste0("go_bp_frame_", v, ".rda")
    gbf <- get_file(file)
    return(gbf)
}


#' Loads GO graph
#'
#' #@examples
#' #load_gobp_net()
#'
#' @return GO graph
#' @import zen4R
#'
load_gobp_net <- function(){
    v <- get_hpannot_version()
    file <- paste0("go_bp_net_", v, ".rda")
    gbn <- get_file(file)
    return(gbn)
}

