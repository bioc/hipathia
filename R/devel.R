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
    node.comp <- data.frame(ID = rowData(hidata[["nodes"]])$name,
                            name = rowData(hidata[["nodes"]])$label,
                            node.comp)

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
    mesdf <- get_measured_nodes(hidata)
    path.comp <- data.frame(ID = rowData(hidata[["paths"]])$path.ID,
                            name = rowData(hidata[["paths"]])$path.name,
                            path.comp,
                            num.nodes = mesdf[rownames(path.comp), "num.nodes"],
                            num.gene.nodes = mesdf[rownames(path.comp), "num.gene.nodes"],
                            num.measured.nodes = mesdf[rownames(path.comp), "num.measured.nodes"],
                            ratio.measured.gene.nodes = mesdf[rownames(path.comp), "ratio.measured.gene.nodes"],
                            nodes = rowData(hidata[["paths"]])$path.nodes)

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
      uni.comp <- data.frame(ID = rownames(assay(hidata[["uni.terms"]])),
                              name = rownames(assay(hidata[["uni.terms"]])),
                             uni.comp)
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
      GO.comp <- data.frame(ID = rownames(rowData(hidata[["GO.terms"]])),
                            GO.comp)
      DAdata$GO.terms <- GO.comp
  }

  return(DAdata)
}


DAsummary <- function(DAdata, conf.level = 0.05){
    # Summary
    summ <- lapply(names(DAdata), function(feat){
        data <- DAdata[[feat]]
        summdf <- data.frame(total = nrow(data),
                             sigs = sum(data$FDRp.value < conf.level),
                             UPs = sum(data$FDRp.value < conf.level & data$statistic > 0),
                             DOWNs = sum(data$FDRp.value < conf.level & data$statistic < 0))
    })
    names(summ) <- names(DAdata)
    summ <- do.call(rbind, summ)
    return(summ)
}

pathway_summary <- function(DAdata, conf = 0.05){
    require(dplyr)
    comp <- DAdata$paths
    comp$pathway.ID <- sapply(strsplit(comp$ID, split = "-"), "[[", 2)
    comp$pathway.name <- sapply(strsplit(comp$name, split = ":"), "[[", 1)
    allpathways <- unique(comp[,c("pathway.ID", "pathway.name")])
    summ <- lapply(allpathways$pathway.ID, function(pathway){
        mini <- comp[comp$pathway.ID == pathway,]
        data.frame(n.sigs = sum(mini$FDRp.value < conf),
                   n.ups = sum(mini$FDRp.value < conf & mini$UP.DOWN == "UP"),
                   n.downs = sum(mini$FDRp.value < conf & mini$UP.DOWN == "DOWN"),
                   n.paths = nrow(mini),
                   p.sigs = sum(mini$FDRp.value < conf)/nrow(mini),
                   p.ups = sum(mini$FDRp.value < conf &
                                         mini$UP.DOWN == "UP")/nrow(mini),
                   p.downs = sum(mini$FDRp.value < conf &
                                           mini$UP.DOWN == "DOWN")/nrow(mini))
    })
    summ <- do.call("rbind", summ)
    summ <- tibble(ID = allpathways$pathway.ID, name = allpathways$pathway.name, summ)
    # rownames(summ) <- summ$ID
    summ <- summ[order(-summ$p.sigs, -summ$n.paths),]
    return(summ)
}
