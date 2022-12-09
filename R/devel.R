##
## devel.R
## Devel functions of package Hipathia
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

#' Compares the gene expression, pathway activation level and the function
#' activation level of the
#'
#' @param hidata Either a SummarizedExperiment object or a matrix, returned by
#' function \code{hipathia}.
#' @param groups Either a character indicating the name of the column in colData
#' including the classes to compare, or a character vector with the class to
#' which each sample belongs.
#' Samples must be ordered as in \code{hidata}.
#' @param expdes String, either an equation expression to pas to \code{limma},
#' or the label of the first group to be compared
#' @param g2 String, label of the second group to be compared, if not specified
#' in \code{expdes}.
#' @param path.method String, method to be used when comparing pathways.
#' Options include \code{wilcoxon} (default, performs a Wilcoxon test comparing
#' conditions \code{expdes} and \code{g2} - in this case, mandatory parameter)
#' and \code{limma} (performs a limma DE analysis using functions \code{lmFit},
#' \code{contrasts.fit} and \code{eBayes} using the formula in \code{expdes} or
#' comparing conditions \code{expdes} and \code{g2}.
#' @param node.method String, method to be used when comparing nodes.
#' Options include \code{wilcoxon} (performs a Wilcoxon test comparing
#' conditions \code{expdes} and \code{g2} - in this case, mandatory parameter)
#' and \code{limma} (default, performs a limma DE analysis using functions
#' \code{lmFit}, \code{contrasts.fit} and \code{eBayes} using the formula in
#' \code{expdes} or comparing conditions \code{expdes} and \code{g2}.
#' @param fun.method String, method to be used when comparing functions.
#' Options include \code{wilcoxon} (default, performs a Wilcoxon test comparing
#' conditions \code{expdes} and \code{g2} - in this case, mandatory parameter)
#' and \code{limma} (performs a limma DE analysis using functions \code{lmFit},
#' \code{contrasts.fit} and \code{eBayes} using the formula in \code{expdes} or
#' comparing conditions \code{expdes} and \code{g2}.
#' @param order Boolean, whether to order the results table by the
#' \code{FDRp.value} column. Default is FALSE.
#' @param paired Boolean, whether the samples to be compared are paired.
#' If TRUE, function \code{wilcoxsign_test} from package \code{coin} is
#' used. If FALSE, function \code{wilcox.test} from package \code{stats}
#' is used.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param conf.level Numeric, cut off for significance. Default is 0.05.
#' @param sel_assay Character or integer, indicating the assay to be normalized
#' in the SummarizedExperiment. Default is 1.
#'
#' @return List including comparison results for nodes, pathways and functions,
#' if present.
#'
#' @examples
#' data(hidata)
#' comp <- DAcomp(hidata, groups = "group", expdes = "Tumor", g2 = "Normal")
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom tibble tibble
#' @importFrom methods is
#'
DAcomp <- function(hidata, groups, expdes, g2 = NULL,
                    path.method = "wilcoxon", node.method = "limma",
                    fun.method = "wilcoxon",
                    order = FALSE, paired = FALSE, adjust = TRUE,
                    conf.level = 0.05, sel_assay = 1){

    # require(dplyr)
  if(is.null(g2) & (any(c(path.method, node.method) == "wilcoxon") |
                    (any(c("uni.terms", "GO.terms") %in% names(hidata)) &
                     fun.method == "wilcoxon")))
    stop("Wilcoxon comparison method needs two groups to compare,
         introduced in arguments expdes and g2 (ex. expdes = 'case', g2 =
         'control'). Please provide both arguments or change comparison method
         to 'limma'.")

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
                        name = rowData(hidata[["nodes"]])$node.name,
                        label = rowData(hidata[["nodes"]])$label,
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
    alt <- get_altered_nodes(hidata,
                             node.comp, conf.level)[rownames(path.comp),]
    path.comp <- tibble(ID = rowData(hidata[["paths"]])$path.ID,
                        name = rowData(hidata[["paths"]])$path.name,
                        path.comp,
                        N.nodes = mesdf$num.nodes,
                        N.gene.nodes = mesdf$num.gene.nodes,
                        N.measured.nodes = mesdf$num.measured.nodes,
                        ratio.measured.gene.nodes =
                            mesdf$ratio.measured.gene.nodes,
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


#' Lists and plots the top \code{n} altered nodes, paths and functions (Uniprot
#' keywords and/or GO terms, if present).
#'
#' @param DAdata List of comparison results, returned by function \code{DAcomp}.
#' @param n Number of top features to show.
#' @param conf.level Numeric, cut off for significance. Default is 0.05.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param colors String with the color scheme or vector of colors to be used.
#' See  \code{define_colors} for available options. Default is "hiro".
#'
#' @return Plot and list of tables including top \code{n} altered features for
#' nodes, paths and functions if present.
#'
#' @examples
#' data(DAdata)
#' DAtop(DAdata)
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr recode_factor
#' @importFrom dplyr mutate
#'
DAtop <- function(DAdata, n = 10, conf.level = 0.05, adjust = TRUE,
                  colors = "hiro"){
    colors <- define_colors(colors)
    toplist <- lapply(names(DAdata), function(feat){
        DA <- DAdata[[feat]]
        if(feat == "nodes") DA$name <- paste(DA$name, "(node)")
        if(adjust == TRUE){
            newn <- min(n, sum(DA$FDRp.value < conf.level))
            DA[order(DA$p.value, decreasing = FALSE),][seq_len(newn),] %>%
                mutate(logPV = abs(log10(FDRp.value)) * sign(statistic),
                       feature = feat)
        }else{
            newn <- min(n, sum(DA$p.value < conf.level))
            DA[order(DA$p.value, decreasing = FALSE),][seq_len(newn),] %>%
                mutate(logPV = abs(log10(p.value)) * sign(statistic),
                       feature = feat)
        }
    })
    names(toplist) <- names(DAdata)
    top <- do.call(rbind, lapply(toplist, function(tl) select(tl, c(name, logPV,
                                                                    feature))))
    top$name <- factor(top$name, levels = top$name[nrow(top):1])
    top$feature <- factor(top$feature,
                          levels = c("nodes", "paths",
                                     names(DAdata)[!names(DAdata) %in%
                                                       c("nodes", "paths")]))
    top$feature <- recode_factor(top$feature, nodes = "Nodes", paths = "Paths",
                                 uni.terms = "Uniprot", GO.terms = "GO terms")
    dir <- c("UP", "DOWN")
    names(dir) <- c("1", "-1")
    top$direction <- dir[paste(sign(top$logPV))]

    print(ggplot(top, aes(x = name, y = logPV, color = direction)) +
              geom_point(stat = "identity") +
              scale_color_manual(name = "Status", values = c(colors$down,
                                                             colors$up)) +
              # scale_fill_met_d("Hiroshige", direction = 1) +
              ylab("abs(Log10 of Adjusted P-value) * direction") +
              xlab("") +
              facet_wrap(~ feature, scales = "free_y") +
              ggtitle(paste("Top", n, "altered features")) +
              theme_minimal() +
              theme(legend.position = "bottom") +
              coord_flip())

    return(toplist)
}

#' Lists and plots the top \code{n} altered pathways, taking into account the
#' number of altered .
#'
#' @param DAdata List of comparison results, returned by function \code{DAcomp}.
#' @param n Number of top features to show.
#' @param conf.level Numeric, cut off for significance. Default is 0.05.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param ratio Boolean, whether to plot the ratio of significant paths with
#' respect to the total paths in the pathway. Default is FALSE.
#' @param colors String with the color scheme or vector of colors to be used.
#' See  \code{define_colors} for available options. Default is "hiro".
#' @param order.by String, how to order table of results. Available options
#' include \code{ratio} (default, uses the ratio of significant paths with
#' respect to the total paths in the pathway) and \code{number} (uses the number
#' of significant paths in the pathway).
#'
#' @return Plot and tibble including top \code{n} altered pathways.
#'
#' @export
#' @examples
#' data(DAdata)
#' DAsummary(DAdata)
#'
DAsummary <- function(DAdata, n = 10, conf.level = 0.05, adjust = TRUE,
                      ratio = FALSE, colors = "hiro", order.by = "number"){
    # Summary
    Psumm <- pathway_summary(DAdata, conf.level, adjust = adjust,
                             order.by = order.by)
    g <- summary_plot(Psumm, n.paths = n, ratio = ratio, colors = colors)
    return(Psumm)
}

#' Table and plot of total number of altered and not altered nodes, paths and
#' functions (Uniprot keywords and/or GO terms, if present).
#'
#' @param DAdata List of comparison results, returned by function \code{DAcomp}.
#' @param conf.level Numeric, cut off for significance. Default is 0.05.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param colors String with the color scheme or vector of colors to be used.
#' See  \code{define_colors} for available options. Default is "hiro".
#'
#' @return Plot and tibble including the number of total, altered, UP- and
#' DOWN-regulated features for nodes, paths and functions if present.
#'
#' @examples
#' data(DAdata)
#' DAoverview(DAdata)
#'
#' @export
#' @importFrom tibble tibble
#'
DAoverview <- function(DAdata, conf.level = 0.05, adjust = TRUE,
                       colors = "hiro"){
    # Summary
    summ <- lapply(names(DAdata), function(feat){
        data <- DAdata[[feat]]
        if(adjust == TRUE){
            summdf <- data.frame(feature = feat,
                                 total = nrow(data),
                                 sigs = sum(data$FDRp.value < conf.level),
                                 UPs = sum(data$FDRp.value < conf.level &
                                               data$statistic > 0),
                                 DOWNs = sum(data$FDRp.value < conf.level &
                                                 data$statistic < 0))
        }else{
            summdf <- data.frame(feature = feat,
                                 total = nrow(data),
                                 sigs = sum(data$p.value < conf.level),
                                 UPs = sum(data$p.value < conf.level &
                                               data$statistic > 0),
                                 DOWNs = sum(data$p.value < conf.level &
                                                 data$statistic < 0))
        }
    })
    summ <- tibble(do.call(rbind, summ))
    p <- nsig_plot(summ, colors)
    return(summ)
}


#' @importFrom tibble tibble
pathway_summary <- function(DAdata, conf = 0.05, adjust = TRUE,
                            order.by = "ratio"){

    # PATHS
    comp <- DAdata$paths
    comp$pathway.ID <- sapply(strsplit(comp$ID, split = "-"), "[[", 2)
    comp$pathway.name <- sapply(strsplit(comp$name, split = ":"), "[[", 1)
    allpathways <- unique(comp[,c("pathway.ID", "pathway.name")])
    summp <- lapply(allpathways$pathway.ID, function(pathway){
        mini <- comp[comp$pathway.ID == pathway,]
        if(adjust == TRUE){
            pdf <- data.frame(sigs = sum(mini$FDRp.value < conf),
                              UPs = sum(mini$FDRp.value < conf &
                                            mini$statistic > 0),
                              DOWNs = sum(mini$FDRp.value < conf &
                                              mini$statistic < 0),
                              total = nrow(mini),
                              ratio.sigs = sum(mini$FDRp.value < conf)/nrow(mini),
                              ratio.UPs = sum(mini$FDRp.value < conf &
                                                  mini$statistic > 0)/nrow(mini),
                              ratio.DOWNs = sum(mini$FDRp.value < conf &
                                                    mini$statistic < 0)/nrow(mini))
        }else{
            pdf <- data.frame(sigs = sum(mini$p.value < conf),
                              UPs = sum(mini$p.value < conf &
                                            mini$statistic > 0),
                              DOWNs = sum(mini$p.value < conf &
                                              mini$statistic < 0),
                              total = nrow(mini),
                              ratio.sigs = sum(mini$p.value < conf)/nrow(mini),
                              ratio.UPs = sum(mini$p.value < conf &
                                                  mini$statistic > 0)/nrow(mini),
                              ratio.DOWNs = sum(mini$p.value < conf &
                                                    mini$statistic < 0)/nrow(mini))
        }
    })
    summp <- do.call("rbind", summp)
    # NODES
    ndata <- DAdata$nodes
    ndata$pathway.ID <- sapply(strsplit(ndata$ID, split = "-"), "[[", 2)
    summn <- lapply(allpathways$pathway.ID, function(pathway){
        mini <- ndata[ndata$pathway.ID == pathway,]
        if(adjust == TRUE){
            ndf <- data.frame(sig.nodes = sum(mini$FDRp.value < conf),
                              UP.nodes = sum(mini$FDRp.value < conf &
                                                 mini$statistic > 0),
                              DOWN.nodes = sum(mini$FDRp.value < conf &
                                                   mini$statistic < 0),
                              gene.nodes = sum(mini$type == "gene"),
                              total.nodes = nrow(mini))
        }else{
            ndf <- data.frame(sig.nodes = sum(mini$p.value < conf),
                              UP.nodes = sum(mini$p.value < conf &
                                                 mini$statistic > 0),
                              DOWN.nodes = sum(mini$p.value < conf &
                                                  mini$statistic < 0),
                              gene.nodes = sum(mini$type == "gene"),
                              total.nodes = nrow(mini))
        }
    })
    summn <- do.call("rbind", summn)
    # TOGETHER
    summ <- tibble(ID = allpathways$pathway.ID,
                   name = allpathways$pathway.name,
                   summp,
                   summn)
    # rownames(summ) <- summ$ID
    if(order.by == "ratio"){
        summ <- summ[order(-summ$ratio.sigs, -summ$sigs),]
    }else if(order.by == "number"){
        summ <- summ[order(-summ$sigs, -summ$ratio.sigs),]
    }else{
        stop("Not supported order.by parameter")
    }
    return(summ)
}


#' @import ggplot2
#' @importFrom MetBrewer scale_fill_met_c
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
summary_plot <- function(Psumm, n.paths = 10, ratio = FALSE, colors = "vg"){

    pdata <- Psumm[seq_len(n.paths),]
    pdata$name <- factor(pdata$name, levels = pdata$name[n.paths:1])

    palette <- define_colors(colors)

    d1 <- mutate(pdata, Not = total - UPs - DOWNs) %>%
        mutate(UP = UPs) %>%
        mutate(DOWN = DOWNs) %>%
        select(c(name, UP, DOWN, Not))
    data1 <- melt(d1, "name")
    data1$variable <- factor(data1$variable,
                             levels = unique(data1$variable)[c(3,1,2)])
    g1 <- ggplot(data1, aes(x = name, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(name = "Status",
                          values = c("#dfe0df", palette$up, palette$down)) +
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
        scale_color_manual(name = "Status",
                           values = c(palette$up, palette$down)) +
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
        data2 <- melt(d2, "name")
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


#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr recode_factor
nsig_plot <- function(summ, colors = "vg"){

    palette <- c("#089099", "#ff8a43", "#5bc6cf", "#befcff") # Colores hipathia
    palette <- define_colors(colors)
    d1 <- mutate(summ, Not = total - UPs - DOWNs) %>%
        mutate(UP = UPs) %>%
        mutate(DOWN = DOWNs) %>%
        select(c(feature, UP, DOWN, Not))
    d1$feature <- factor(d1$feature,
                         levels = c("nodes", "paths",
                                    d1$feature[!d1$feature %in%
                                                   c("nodes", "paths")]))
    d1$feature <- recode_factor(d1$feature, nodes = "Nodes", paths = "Paths",
                  uni.terms = "Uniprot", GO.terms = "GO terms")

    data1 <- melt(d1, "feature")
    # data1$feature <- factor(data1$feature, levels = levels(data1$feature)[length(levels(data1$feature)):1])
    data1$variable <- factor(data1$variable,
                             levels = unique(data1$variable)[c(3,1,2)])
    g <- ggplot(data1, aes(x = feature, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(name = "Status",
                          values = c("#dfe0df", palette$up, palette$down)) +
        ylab("") +
        xlab("Feature") +
        ggtitle("Results overview") +
        theme_minimal() +
        theme(legend.position = "right") #+
        # coord_flip()

    print(g)
    return(g)
}

#' Color palettes to be used in plots.
#'
#' @param colors String with the color scheme or vector of colors to be used.
#' Available predefined options include: \code{hipathia}, \code{classic},
#' \code{soft}, \code{okee}, \code{hiro}, \code{new}, \code{vg}, \code{orchid}.
#'
#' @param no.col String with the color given to non-significant nodes, if not
#' given in parameter \code{colors}.
#'
#' @return Plot and list of tables including top \code{n} altered features for
#' nodes, paths and functions if present.
#'
#' @examples
#' define_colors("hiro")
#'
#' @export
#' @importFrom grDevices rgb
#' @importFrom grDevices colorRamp
#'
define_colors <- function(colors, no.col = NULL){
    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#50b7ae", "white", "#f16a34")
        }else if(colors == "classic"){
            colors <- c("#1f9cda","white","#da1f1f")
        }else if(colors == "soft"){
            colors <- c("#50B7AE", "#ffa17a")
        }else if(colors == "okee"){
            colors <- c("#447fdd", "#da6c42")
        }else if(colors == "hiro"){
            colors <- c("#72bcd5", "#f7aa58", "#e76254")
        }else if(colors == "new"){
            colors <- c("#089099", "#eee8a9", "#ff8a43")
        }else if(colors == "vg"){
            colors <- c("#60a8ff", "#ff9368")
        }else if(colors == "orchid"){
            colors <- c("#9da6ef", "#d8443c")
        }
    }
    down_col <- colors[1]
    no_col <- ifelse(is.null(no.col), colors[2], no.col)
    up_col <- ifelse(length(colors) == 3, colors[3], colors[2])
    both <- rgb(colorRamp(c(up_col, down_col))(0.5)/256)
    return(list(down = down_col, no = no_col, up = up_col, both = both))
}

#' @importFrom tibble tibble
get_edges_df <- function(g){
    tibble(from = get.edgelist(g)[,1],
               to = get.edgelist(g)[,2]) %>%
        mutate(name = paste(from, to, sep = "_"))
}

#' @importFrom dplyr filter
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
    comp <- filter(DApaths, isname)
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

#' @importFrom dplyr mutate
prepare_DAedges <- function(DApaths, name, pathways, cols, conf = 0.05,
                            adjust = TRUE){
    # require(dplyr)
    pg <- pathways$pathigraphs[[name]]

    # Define colors
    color.edge.type <- c(cols$up, cols$down, cols$both, "lightgray",
                         "gainsboro") # c(met.brewer("Egypt", 4), "gainsboro") # c("#0571b0", "green", "#ca0020", "#ffc868", "gainsboro")
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

#' @importFrom dplyr mutate
prepare_edges <- function(name, pathways,conf = 0.05, adjust = TRUE){
    # require(dplyr)
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

#' @importFrom tibble tibble
#' @importFrom dplyr filter
prepare_DAnodes <- function(DAdata, name, pathways, cols,
                            conf = 0.05, adjust = TRUE, no.col = NULL){

    DAnodes <- DAdata[["nodes"]]
    DApaths <- DAdata[["paths"]]
    g <- pathways$pathigraphs[[name]]$graph

    isname <- sapply(strsplit(DAnodes$ID, "-"), "[[", 2) == name
    DAnodes <- filter(DAnodes, isname)
    isname <- sapply(strsplit(DApaths$ID, "-"), "[[", 2) == name
    DApaths <- filter(DApaths, isname)
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

#' @importFrom tibble tibble
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


#' Plots a pathway with or without the comparison information, using the
#' visNetwork library.
#'
#' @param name KEGG ID of the pathway to plot.
#' @param pathways Pathways object.
#' @param DAdata List of comparison results, returned by function \code{DAcomp}.
#' @param colors String with the color scheme or vector of colors to be used.
#' See  \code{define_colors} for available options. Default is "hiro".
#' @param conf Numeric, cut off for significance. Default is 0.05.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param main Title of the plot.
#' @param submain Subtitle of the plot.
#' @param no.col String with the color given to non-significant nodes.
#' @param height Height of the plot. Default is "800px".
#'
#' @return Plot of the pathway.
#'
#' @examples
#' data(pathways)
#' plotVG("hsa03320", pathways)
#'
#' data(DAdata)
#' plotVG("hsa04012", pathways, DAdata)
#'
#' @import visNetwork
#' @export
#'
plotVG <- function(name, pathways, DAdata = NULL, colors = "hiro",
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
        nodes <- prepare_DAnodes(DAdata, name, pathways, cols, conf, adjust,
                                 no.col)
        edges <- prepare_DAedges(DAdata[["paths"]], name, pathways, cols, conf,
                                 adjust)
        submain <- "Differential activation plot"
        ledges <- data.frame(label = c("UP", "DOWN", "Both", "None", "function"),
                             color = c(cols$up, cols$down, cols$both,
                                       "lightgray", "gainsboro"),
                             width = c(10, 10, 10, 10, 1))
    }

    pname <- pathways$pathigraphs[[name]]$path.name

    plotVisGraphDE(nodes, edges, ledges, main = pname, submain = submain,
                   cols = cols, height = height)
}


#' @import visNetwork
plotVisGraphDE <- function(nodes, edges, ledges, main = "Pathway",
                           submain = "Differential activation plot",
                           cols = list(no = "BlanchedAlmond", up = "red",
                                       down = "blue"),
                           height = "800px"){
    # require(visNetwork, quietly = TRUE)

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
                                           hover = FALSE,
                                           labelOnly = FALSE),
                   # nodesIdSelection = list(enabled = TRUE,
                   #                         main = "Select by gene",
                   #                         values = nodes$label[groups == "gene"]),
                   # selectedBy = list(variable = "label",
                   #                   main = "Select by function",
                   #                   values = unique(nodes$label[groups == "function"]),
                   #                   multiple = TRUE)
        ) %>%
        visLegend(position = "right", useGroups = TRUE, main = "Legend",
                  addEdges = ledges)

}


#' Create visualization HTML
#'
#' Saves the results of a DAdata comparison for the Hipathia pathway values
#' into a folder, and creates a HTML from which to visualize the results on
#' top of the pathways. The results are stored into the specified folder.
#' If this folder does not exist, it will be created. The parent folder must
#' exist.
#'
#' @examples
#' data(DAdata)
#' data(pathways)
#' DAreport(DAdata, pathways)
#'
#' @param DAdata List of comparison results, returned by function \code{DAcomp}.
#' @param pathways Pathways object as returned by the \code{load_pathways}
#' function
#' @param conf.level Level of significance. By default 0.05.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method. Default is TRUE.
#' @param group_by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include. Default is set to "pathway".
#' @param colors String with the color scheme or vector of colors to be used.
#' See  \code{define_colors} for available options. Default is "hiro".
#' @param output_folder Name of the folder in which the report will be stored.
#' @param path Absolute path to the parent directory in which `output_folder`
#' will be saved. If it is not provided, it will be created in a temp folder.
#' @param verbose Boolean, whether to show details about the results of the
#' execution
#'
#' @return Saves the results and creates a report to visualize them through
#' a server in the specified \code{output_folder}. Returns the folder where
#' the report has been stored.
#'
#' @export
#'
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
