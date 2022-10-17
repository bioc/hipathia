DApathways <- function(data, groups, expdes, g2 = NULL, 
                       
                    path.method = "wilcoxon", node.method = "limma", 
                    order = FALSE, paired = FALSE, adjust = TRUE,
                    conf.level = 0.05, sel_assay = 1){
  
  # Pathways comparison
  if(path.method == "wilcoxon"){
    path.comp <- do_wilcoxon(data = data, 
                             group = groups, 
                             g1 = expdes, 
                             g2 = g2, 
                             paired = paired, 
                             adjust = adjust, 
                             sel_assay = sel_assay, 
                             order = order)
  }else if(path.method == "limma"){
    path.comp <- do_limma(data = data, 
                          groups = groups, 
                          expdes = expdes, 
                          g2 = g2, 
                          sel_assay = sel_assay, 
                          order = order)
  }
   
  # Node comparison
  
  
}
