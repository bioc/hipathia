DApathways <- function(hidata, groups, expdes, g2 = NULL,
                    path.method = "wilcoxon", node.method = "limma",
                    fun.method = "wilcoxon",
                    order = FALSE, paired = FALSE, adjust = TRUE,
                    conf.level = 0.05, sel_assay = 1){

    require(dplyr)
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
    node.comp <- tibble(ID = rowData(hidata[["nodes"]])$name,
                        name = rowData(hidata[["nodes"]])$label,
                        node.comp,
                        type = rowData(hidata[["nodes"]])$node.type)

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
    mesdf <- get_measured_nodes(hidata)[rownames(path.comp),]
    alt <- get_altered_nodes(hidata, node.comp, conf.level)[rownames(path.comp),]
    path.comp <- tibble(ID = rowData(hidata[["paths"]])$path.ID,
                        name = rowData(hidata[["paths"]])$path.name,
                        path.comp,
                        N.nodes = mesdf$num.nodes,
                        N.gene.nodes = mesdf$num.gene.nodes,
                        N.measured.nodes = mesdf$num.measured.nodes,
                        ratio.measured.gene.nodes = mesdf$ratio.measured.gene.nodes,
                        nodes = rowData(hidata[["paths"]])$path.nodes,
                        N.DA.nodes = alt$N.DA.nodes,
                        DA.nodes = alt$DA.nodes)

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
      uni.comp <- tibble(ID = rownames(assay(hidata[["uni.terms"]])),
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
      GO.comp <- tibble(ID = rownames(rowData(hidata[["GO.terms"]])),
                            GO.comp)
      DAdata$GO.terms <- GO.comp
  }

  return(DAdata)
}


DAsummary <- function(DAdata, plot = FALSE, conf.level = 0.05,
                      n.paths = 10, ratio = F){
    # Summary
    summ <- lapply(names(DAdata), function(feat){
        data <- DAdata[[feat]]
        summdf <- data.frame(feature = feat,
                             total = nrow(data),
                             sigs = sum(data$FDRp.value < conf.level),
                             UPs = sum(data$FDRp.value < conf.level & data$statistic > 0),
                             DOWNs = sum(data$FDRp.value < conf.level & data$statistic < 0))
    })
    summ <- tibble(do.call(rbind, summ))
    Psumm <- pathway_summary(DAdata, conf.level)
    DAsumm <- list(all = summ, by.pathway = Psumm)
    p <- nsig_plot(summ)
    g <- summary_plot(Psumm, n.paths = n.paths, ratio = ratio)
    if(plot == TRUE){
        DAsumm$nsig_plot <- p
        DAsumm$summary_plot <- g
    }
    return(DAsumm)
}

pathway_summary <- function(DAdata, conf = 0.05){
    require(dplyr)
    # PATHS
    comp <- DAdata$paths
    comp$pathway.ID <- sapply(strsplit(comp$ID, split = "-"), "[[", 2)
    comp$pathway.name <- sapply(strsplit(comp$name, split = ":"), "[[", 1)
    allpathways <- unique(comp[,c("pathway.ID", "pathway.name")])
    summp <- lapply(allpathways$pathway.ID, function(pathway){
        mini <- comp[comp$pathway.ID == pathway,]
        pdf <- data.frame(sigs = sum(mini$FDRp.value < conf),
                   UPs = sum(mini$FDRp.value < conf & mini$statistic > 0),
                   DOWNs = sum(mini$FDRp.value < conf & mini$statistic < 0),
                   total = nrow(mini),
                   ratio.sigs = sum(mini$FDRp.value < conf)/nrow(mini),
                   ratio.UPs = sum(mini$FDRp.value < conf &
                                   mini$statistic > 0)/nrow(mini),
                   ratio.DOWNs = sum(mini$FDRp.value < conf &
                                     mini$statistic < 0)/nrow(mini))
    })
    summp <- do.call("rbind", summp)
    # NODES
    ndata <- DAdata$nodes
    ndata$pathway.ID <- sapply(strsplit(ndata$ID, split = "-"), "[[", 2)
    summn <- lapply(allpathways$pathway.ID, function(pathway){
        mini <- ndata[ndata$pathway.ID == pathway,]
        ndf <- data.frame(sig.nodes = sum(mini$FDRp.value < conf),
                          UP.nodes = sum(mini$FDRp.value < conf & mini$statistic > 0),
                          DOWN.nodes = sum(mini$FDRp.value < conf & mini$statistic < 0),
                          gene.nodes = sum(mini$type == "gene"),
                          total.nodes = nrow(mini))
    })
    summn <- do.call("rbind", summn)
    # TOGETHER
    summ <- tibble(ID = allpathways$pathway.ID,
                   name = allpathways$pathway.name,
                   summp,
                   summn)
    # rownames(summ) <- summ$ID
    summ <- summ[order(-summ$ratio.sigs, -summ$total),]
    return(summ)
}


summary_plot <- function(Psumm, n.paths = 10, ratio = F){
    require(ggplot2)
    require(reshape2)
    require(ggpubr)
    require(MetBrewer)
    require(rcartocolor)
    require(RColorBrewer)

    pdata <- Psumm[1:n.paths,]
    pdata$name <- factor(pdata$name, levels = pdata$name[n.paths:1])

    palette <- c("#da1f1f", "#1f9cda") # Clásico hipathia
    palette <- c("#ff8a43", "#089099") # Colores hipathia
    palette <- c("#c7e9b4", "#41b6c4") # Frío bonito
    palette <- met.brewer(name = "Signac", 2, direction = 1) # Cálido bonito

    d1 <- mutate(pdata, Not = total - UPs - DOWNs) %>%
        mutate(UP = UPs) %>%
        mutate(DOWN = DOWNs) %>%
        select(c(name, UP, DOWN, Not))
    data1 <- melt(d1)
    data1$variable <- factor(data1$variable, levels = unique(data1$variable)[c(3,1,2)])
    g1 <- ggplot(data1, aes(x = name, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(name = "Status", values = c("#dfe0df", palette)) +
        # scale_fill_met_d("Hiroshige", direction = 1) +
        ylab("Total significant paths") +
        xlab("Pathway") +
        ggtitle("Top altered pathways") +
        theme_minimal() +
        theme(legend.position = "left") +
        coord_flip()

    d3 <- select(pdata, c(name, UP.nodes, DOWN.nodes))
    colnames(d3) <- c("name", "UP", "DOWN")
    data3 <- melt(d3, "name")
    data3$nodes <- pdata$gene.nodes
    g3 <- ggplot(data3, aes(x = name, y = variable)) +
        geom_point(aes(color = variable, size = nodes)) +
        geom_point(aes(size = nodes - 5), color = "white") +
        geom_point(aes(color = variable, size = value)) +
        scale_color_manual(name = "Status", values = palette) +
        ylab("DE nodes") +
        ggtitle("") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()) +
        coord_flip()
    g1 + geom_point(aes(x = "UP", size = select(d3, UP)))

    if(ratio == TRUE){
        d2 <- select(pdata, c(name, ratio.sigs, ratio.UPs, ratio.DOWNs))
        colnames(d2) <- c("name", "Sig", "UP", "DOWN")
        data2 <- melt(d2)
        g2 <- ggplot(data2, aes(x = variable, y = name, fill = value)) +
            geom_tile() +
            # scale_fill_distiller(palette = "YlOrBr", direction = 1) +
            scale_fill_met_c("Hokusai3", direction = 1) +
            theme_minimal() +
            ggtitle("") +
            theme(legend.position = "top") +
            xlab("Ratio SP") +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank())
        all <- ggarrange(g1, g2, g3, widths = c(0.7, 0.15, 0.15), ncol = 3,
                         common.legend = TRUE, legend = "right")
    }else{
        all <- ggarrange(g1, g3, widths = c(0.8, 0.2), ncol = 2,
                         common.legend = TRUE, legend = "right")
    }
    print(all)
}


nsig_plot <- function(summ){
    require(ggplot2)
    require(reshape2)

    palette <- c("#089099", "#ff8a43", "#5bc6cf", "#befcff") # Colores hipathia

    colnames(summ) <- c("feature", "Total", "Sig", "UP", "DOWN")
    data <- melt(summ)
    g <- ggplot(data, aes(x = feature, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(name = "Type", values = palette) +
        ylab("") +
        xlab("Feature") +
        ggtitle("Results summary") +
        theme_minimal() +
        theme(legend.position = "left")

    print(g)
    return(g)
}

define_colors <- function(colors, no.col = NULL){
    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#50b7ae", "white", "#f16a34")
        }else if(colors == "classic"){
            colors <- c("#1f9cda","white","#da1f1f")
        }else if(colors == "okee"){
            colors <- met.brewer("OKeeffe1", 7, direction = -1)[c(2,4,6)]
        }else if(colors == "hiro"){
            colors <- met.brewer("Hiroshige", 9, direction = -1)[c(3,7,9)]
        }else if(colors == "new"){
            colors <- c("#089099", "#eee8a9", "#ff8a43")
        }else if(colors == "edges"){
            colors <- c("#447fdd", "gray", "#da6c42")
        }else if(colors == "vg"){
            colors <- c("#60a8ff", "#ff9368")
        }
    }
    down_col <- colors[1]
    no_col <- ifelse(is.null(no.col), colors[2], no.col)
    up_col <- ifelse(length(colors) == 3, colors[3], colors[2])
    both <- rgb(colorRamp(c(up_col, down_col))(0.5)/256)
    return(list(down = down_col, no = no_col, up = up_col, both = both))
}

get_edges_df <- function(g){
    require(dplyr)
    tibble(from = get.edgelist(g)[,1],
               to = get.edgelist(g)[,2]) %>%
        mutate(name = paste(from, to, sep = "_"))
}

get_edges_status <- function(pg, edgename, DApaths, adjust = TRUE){

    # Compute matrix of edge inclusion per subpath
    edgesinsub <- sapply(pg$effector.subgraphs, function(es){
        esell <- get_edges_df(es)
        return(edgename %in% esell$name)
    })
    rownames(edgesinsub) <- edgename

    # Filter DA paths in this pathway & stablish their UP/DOWN status
    name <- pg$path.id
    isname <- sapply(strsplit(DApaths$ID, "-"), "[[", 2) == name
    comp <- dplyr::filter(DApaths, isname)
    if(adjust == TRUE){pv <- comp$FDRp.value}else{pv <- comp$p.value}
    pathsig <- ifelse(pv < 0.05, comp$`UP/DOWN`, "NOT")
    names(pathsig) <- comp$ID

    # Summarize edge status across all including paths
    edgestatus <- apply(edgesinsub, 1, function(v){
        status <- pathsig[names(v)[v == TRUE]]
        if("UP" %in% status & !"DOWN" %in% status)
            return("UP")
        if("DOWN" %in% status & !"UP" %in% status)
            return("DOWN")
        if("UP" %in% status & "DOWN" %in% status)
            return("Both")
        if(!"UP" %in% status & !"DOWN" %in% status)
            return("None")
    })
    return(edgestatus)
}

prepare_DAedges <- function(DApaths, name, pathways, cols, conf = 0.05, adjust = TRUE){
    require(dplyr)
    pg <- pathways$pathigraphs[[name]]

    # Define colors
    color.edge.type <- c(cols$up, cols$down, cols$both, "lightgray", "gainsboro") # c(met.brewer("Egypt", 4), "gainsboro") # c("#0571b0", "green", "#ca0020", "#ffc868", "gainsboro")
    names(color.edge.type) <- c("UP", "DOWN", "Both", "None", "function")

    # Create edges tibble
    edges <- get_edges_df(pg$graph) %>%
        mutate(status = get_edges_status(pg, name, DApaths, adjust),
               functional = grepl("_func", to),
               type = ifelse(functional, "function", status),
               color = color.edge.type[type],
               width = ifelse(functional, 1, 10),
               arrows = ifelse(functional, "none", "to"),
               # hidden = functional,
               dashed = E(pg$graph)$relation == -1)
    return(edges)
}

prepare_edges <- function(name, pathways,conf = 0.05, adjust = TRUE){
    require(dplyr)
    pg <- pathways$pathigraphs[[name]]

    # Create edges tibble
    edges <- get_edges_df(pg$graph) %>%
        mutate(functional = grepl("_func", to),
               width = ifelse(functional, 1, 10),
               color = "lightgray",
               arrows = ifelse(functional, "none", "to"),
               dashed = E(pg$graph)$relation == -1)
    return(edges)
}

prepare_DAnodes <- function(DAdata, name, pathways, cols,
                            conf = 0.05, adjust = TRUE, no.col = NULL){

    DAnodes <- DAdata[["nodes"]]
    DApaths <- DAdata[["paths"]]
    g <- pathways$pathigraphs[[name]]$graph

    isname <- sapply(strsplit(DAnodes$ID, "-"), "[[", 2) == name
    DAnodes <- dplyr::filter(DAnodes, isname)
    isname <- sapply(strsplit(DApaths$ID, "-"), "[[", 2) == name
    DApaths <- dplyr::filter(DApaths, isname)
    # if(adjust == TRUE){
    #     pathsUP <- DApaths %>% filter(FDRp.value < conf, statistic > 0) %>% select(ID)
    #     pathsDOWN <- DApaths %>% filter(FDRp.value < conf, statistic < 0) %>% select(ID)
    # }else{
    #     pathsUP <- DApaths %>% filter(FDRp.value < conf, statistic > 0) %>% select(ID)
    #     pathsDOWN <- DApaths %>% filter(FDRp.value < conf, statistic < 0) %>% select(ID)
    # }
    # pathsUP <- gsub("P", "N", unlist(pathsUP))
    # pathsDOWN <- gsub("P", "N", unlist(pathsDOWN))
    effectors <- gsub("P", "N", unlist(DApaths$ID))

    if(adjust == TRUE){pv <- DAnodes$FDRp.value}else{pv <- DAnodes$p.value}
    color <- get_colors_from_pval(tolower(DAnodes$`UP/DOWN`),
                                  pv,
                                  up_col = cols$up,
                                  down_col = cols$down,
                                  no_col = cols$no,
                                  conf = conf)
    names(color) <- DAnodes$ID
    toadd <- V(g)$name[!V(g)$name %in% names(color)]
    coltoadd <- rep("white", length(toadd))
    names(coltoadd) <- toadd
    color <- c(color, coltoadd)
    color <- color[V(g)$name]

    # color.border <- color
    # color.border[pathsUP] <- "#812804"
    # color.border[pathsDOWN] <- "#005dcf"
    #
    # border.width <- rep(1, length(V(g)))
    # names(border.width) <- V(g)$name
    # border.width[pathsUP] <- 5
    # border.width[pathsDOWN] <- 5

    isUP <- sapply(V(g)$name, function(n){
        if(n %in% DAnodes$ID){
            if(adjust == TRUE){
                filter(DAnodes, ID == n) %>% select(FDRp.value) < conf &
                    filter(DAnodes, ID == n) %>% select(statistic) > 0
            }else{
                filter(DAnodes, ID == n) %>% select(p.value) < conf &
                    filter(DAnodes, ID == n) %>% select(statistic) > 0
            }
        }else{FALSE}
    })
    isDOWN <- sapply(V(g)$name, function(n){
        if(n %in% DAnodes$ID){
            if(adjust == TRUE){
                filter(DAnodes, ID == n) %>% select(FDRp.value) < conf &
                    filter(DAnodes, ID == n) %>% select(statistic) < 0
            }else{
                filter(DAnodes, ID == n) %>% select(p.value) < conf &
                    filter(DAnodes, ID == n) %>% select(statistic) < 0
            }
        }else{FALSE}
    })
    group <- rep("gene", length(V(g)))
    names(group) <- V(g)$name
    group[isUP] <- "gene UP"
    group[isDOWN] <- "gene DOWN"
    group[effectors] <- "effector"
    group[names(group) %in% effectors & isUP] <- "effector UP"
    group[names(group) %in% effectors & isDOWN] <- "effector DOWN"
    group[V(g)$shape == "circle"] <- "metabolite"
    group[grepl("_func", V(g)$name)] <- "function"

    nodes <- tibble(id = V(g)$name,
                    group = group,
                    label = V(g)$label,
                    title = V(g)$tooltip,
                    size = 10,
                    color = color,
                    # color.border = color.border,
                    # borderWidth = border.width,
                    font.size = 35,
                    x = V(g)$x,
                    y = V(g)$y)
    return(nodes)
}

prepare_nodes <- function(name, pathways, conf = 0.05, adjust = TRUE,
                          no.col = "BlanchedAlmond"){

    g <- pathways$pathigraphs[[name]]$graph

    group <- rep("gene", length(V(g)))
    group[V(g)$shape == "circle"] <- "metabolite"
    group[grepl("_func", V(g)$name)] <- "function"

    nodes <- tibble(id = V(g)$name,
                    group = group,
                    label = V(g)$label,
                    title = V(g)$tooltip,
                    color = no.col,
                    color.border = no.col,
                    size = 10,
                    font.size = 35,
                    x = V(g)$x,
                    y = V(g)$y)
    return(nodes)
}


plotVG <- function(name, pathways, DAdata = NULL, colors = "vg",
                   conf = 0.05, adjust = TRUE, main = "Pathway",
                   submain = "", no.col = "BlanchedAlmond",
                   height = "800px"){

    cols <- define_colors(colors, no.col)

    if(is.null(DAdata)){
        nodes <- prepare_nodes(name, pathways, conf, adjust, no.col)
        edges <- prepare_edges(name, pathways, conf, adjust)
        submain <- "KEGG database"
        ledges <- data.frame(label = c("Relation", "function"),
                             color = c("lightgray", "gainsboro"),
                             width = c(10, 1))
    }else{
        nodes <- prepare_DAnodes(DAdata, name, pathways, cols, conf, adjust, no.col)
        edges <- prepare_DAedges(DAdata[["paths"]], name, pathways, cols, conf, adjust)
        submain <- "Differential activation plot"
        ledges <- data.frame(label = c("UP", "DOWN", "Both", "None", "function"),
                             color = c(cols$up, cols$down, cols$both, "lightgray", "gainsboro"),
                             width = c(10, 10, 10, 10, 1))
    }

    pname <- pathways$pathigraphs[[name]]$path.name

    plotVisGraphDE(nodes, edges, ledges, main = pname, submain = submain,
                   cols = cols, height = height)
}


plotVisGraphDE <- function(nodes, edges, ledges, main = "Pathway",
                           submain = "Differential activation plot",
                           cols = list(no = "BlanchedAlmond", up = "red", down = "blue"),
                           height = "800px"){
    require(visNetwork, quietly = TRUE)

    coords <- matrix(c(nodes$x, -nodes$y), ncol = 2)

    visNetwork(nodes, edges, height = height, width = "100%",
               main = main, submain = submain) %>%
        visGroups(groupname = "gene",
                  shape = "box",
                  color = list(background = cols$no,
                               border = cols$no,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "#6d7698"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "effector",
                  shape = "elipse",
                  color = list(background = cols$no,
                               border = cols$no,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "#6d7698"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "gene UP",
                  shape = "box",
                  color = list(background = cols$up,
                               border = cols$up,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "white"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "effector UP",
                  shape = "elipse",
                  color = list(background = cols$up,
                               border = cols$up,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "white"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "gene DOWN",
                  shape = "box",
                  color = list(background = cols$down,
                               border = cols$down,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "white"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "effector DOWN",
                  shape = "elipse",
                  color = list(background = cols$down,
                               border = cols$down,
                               highlight = list(background = "#ff9368",
                                                border = "#e4882e")),
                  font = list(color = "white"),
                  labelHighlightBold = FALSE,
                  shadow = list(enabled = TRUE,
                                size = 5)) %>%
        visGroups(groupname = "metabolite",
                  shape = "square",
                  color = cols$no,
                  shadow = FALSE) %>%
        visGroups(groupname = "function",
                  shape = "text",
                  font = list(color = "darkgray", align = "right"),
                  shadow = FALSE) %>%
        visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
        visOptions(highlightNearest = list(enabled = TRUE,
                                           degree = 100,
                                           algorithm = "hierarchical",
                                           hover = F,
                                           labelOnly = F),
                   # nodesIdSelection = list(enabled = TRUE,
                   #                         main = "Select by gene",
                   #                         values = nodes$label[groups == "gene"]),
                   # selectedBy = list(variable = "label",
                   #                   main = "Select by function",
                   #                   values = unique(nodes$label[groups == "function"]),
                   #                   multiple = TRUE)
        ) %>%
        visLegend(position = "right", useGroups = T, main = "Legend",
                  addEdges = ledges)

}


DAreport <- function(DAdata, pathways, conf.level = 0.05, adjust = TRUE,
                     group_by = "pathway", colors = "classic",
                     output_folder = NULL, path = NULL, verbose = TRUE){

    nodecomp <- as.data.frame(DAdata[["nodes"]])
    rownames(nodecomp) <- nodecomp$ID
    colors_de <- node_color(nodecomp, pathways, group_by = group_by,
                            colors = colors, conf = conf.level,
                            adjust = adjust)

    comp <- as.data.frame(select(DAdata[["paths"]], ID:FDRp.value))
    rownames(comp) <- comp$ID
    path <- create_report(comp, pathways, node_colors = colors_de,
                          output_folder = output_folder, path = path,
                          verbose = verbose)
    return(path)
}