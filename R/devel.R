DApathways <- function(hidata, groups, expdes, g2 = NULL, 
                    path.method = "wilcoxon", node.method = "limma", 
                    fun.method = "wilcoxon",
                    order = FALSE, paired = FALSE, adjust = TRUE,
                    conf.level = 0.05, sel_assay = 1){

  if(is.null(g2) & (any(c(path.method, node.method) == "wilcoxon") | 
                    (any(c("uni.terms", "GO.terms") %in% names(hidata)) & 
                     fun.method == "wilcoxon")))
    stop("Wilcoxon comparison method needs two groups to compare, 
         introduced in arguments expdes and g2 (ex. expdes = 'case', g2 = 'control').
         Please provide both arguments or change comparison method to 'limma'.")
  
  # Node comparison
  if(node.method == "wilcoxon"){
    node.comp <- do_wilcoxon(data = hidata[["nodes"]], 
                             group = groups, 
                             g1 = expdes, 
                             g2 = g2, 
                             paired = paired, 
                             adjust = adjust, 
                             sel_assay = sel_assay, 
                             order = order)
  }else if(node.method == "limma"){
    node.comp <- do_limma(data = hidata[["nodes"]], 
                          groups = groups, 
                          expdes = expdes, 
                          g2 = g2, 
                          sel_assay = sel_assay, 
                          order = order)
  }
  
  # Pathways comparison
  if(path.method == "wilcoxon"){
    path.comp <- do_wilcoxon(data = hidata[["paths"]], 
                             group = groups, 
                             g1 = expdes, 
                             g2 = g2, 
                             paired = paired, 
                             adjust = adjust, 
                             sel_assay = sel_assay, 
                             order = order)
  }else if(path.method == "limma"){
    path.comp <- do_limma(data = hidata[["paths"]], 
                          groups = groups, 
                          expdes = expdes, 
                          g2 = g2, 
                          sel_assay = sel_assay, 
                          order = order)
  }
  path.comp$measured.nodes <- get_measured_nodes(hidata, rownames(path.comp))
  
  DAdata <- list(nodes = node.comp, paths = path.comp)
  
  # Uniprot comparison
  if("uni.terms" %in% names(hidata)){
    if(fun.method == "wilcoxon"){
      uni.comp <- do_wilcoxon(data = hidata[["uni.terms"]], 
                               group = groups, 
                               g1 = expdes, 
                               g2 = g2, 
                               paired = paired, 
                               adjust = adjust, 
                               sel_assay = sel_assay, 
                               order = order)
    }else if(fun.method == "limma"){
      uni.comp <- do_limma(data = hidata[["uni.terms"]], 
                            groups = groups, 
                            expdes = expdes, 
                            g2 = g2, 
                            sel_assay = sel_assay, 
                            order = order)
    }
    DAdata$uni.terms <- uni.comp
  }
  
  # GO comparison
  if("GO.terms" %in% names(hidata)){
    if(fun.method == "wilcoxon"){
      GO.comp <- do_wilcoxon(data = hidata[["GO.terms"]], 
                              group = groups, 
                              g1 = expdes, 
                              g2 = g2, 
                              paired = paired, 
                              adjust = adjust, 
                              sel_assay = sel_assay, 
                              order = order)
    }else if(fun.method == "limma"){
      GO.comp <- do_limma(data = hidata[["GO.terms"]], 
                           groups = groups, 
                           expdes = expdes, 
                           g2 = g2, 
                           sel_assay = sel_assay, 
                           order = order)
    }
    DAdata$GO.terms <- GO.comp
  }
  
  return(DAdata)
}
