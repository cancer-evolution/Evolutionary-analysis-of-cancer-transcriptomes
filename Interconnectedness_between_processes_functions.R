
find_pathway_network <- function(pathway_name){
  local_pathway <- pathways[[pathway_name]]
  pathway_network <- network[intersect(which(network$Gene1 %in% local_pathway), which(network$Gene2 %in% local_pathway)),]
  return(pathway_network)
}

find_minimum_distance <- function(local_pathway2, genes_in_local_pathway1, shortest_paths_matrix){
  genes_in_local_pathway2 <- pathways[[local_pathway2]]
  return(min(shortest_paths_matrix[which(rownames(shortest_paths_matrix) %in% genes_in_local_pathway2),
                                   which(colnames(shortest_paths_matrix) %in% genes_in_local_pathway1)]))
}

identify_overlap <- function(pathway_name, genes_in_local){
  genes_in_pathway <- pathways[[pathway_name]]
  return(sum(genes_in_local %in% genes_in_pathway))
}

find_mode_of_distance_1 <- function(local_pathway2, genes_in_local_pathway1, shortest_paths_matrix){
  genes_in_local_pathway2 <- pathways[[local_pathway2]]
  local_matrix <- shortest_paths_matrix[which(rownames(shortest_paths_matrix) %in% genes_in_local_pathway2),
                                        which(colnames(shortest_paths_matrix) %in% genes_in_local_pathway1)]
  
  return(sum(local_matrix == 1))
}

find_mode_of_distance_1v2 <- function(local_pathway2, genes_in_local_pathway1, shortest_paths_matrix){
  genes_in_local_pathway2 <- pathways[[local_pathway2]]
  local_matrix <- shortest_paths_matrix[which(rownames(shortest_paths_matrix) %in% genes_in_local_pathway2),
                                        which(colnames(shortest_paths_matrix) %in% genes_in_local_pathway1)]
  if(class(local_matrix) == "matrix" && nrow(local_matrix) != 0 && ncol(local_matrix) != 0){
    coordinates_distance_1 <- which(local_matrix==1, arr.ind=T)
    edges_distance_1 <- cbind(rownames(local_matrix)[coordinates_distance_1[,1]],
                              colnames(local_matrix)[coordinates_distance_1[,2]])
    
    return(edges_distance_1) 
  } else {
    return(cbind(rownames(shortest_paths_matrix)[which(rownames(shortest_paths_matrix) %in% genes_in_local_pathway2)],
                 colnames(shortest_paths_matrix)[which(colnames(shortest_paths_matrix) %in% genes_in_local_pathway1)])) 
  }
}

find_mode_of_distance_1v3 <- function(local_pathway2, genes_in_local_pathway1, network_v2){
  genes_in_local_pathway2 <- pathways[[local_pathway2]]
  
  local_network_1 <- network_v2[intersect(which(network_v2[,1] %in% genes_in_local_pathway1),
                                        which(network_v2[,2] %in% genes_in_local_pathway2)),]
  
  local_network_2 <- network_v2[intersect(which(network_v2[,2] %in% genes_in_local_pathway1),
                                          which(network_v2[,1] %in% genes_in_local_pathway2)),]
  
  edges_distance_1 <- rbind(local_network_1, local_network_2)
  
  return(edges_distance_1) 
 
}


get_pathway_components_summary <- function (pathways){
  pathway_components <- vector()
  for (local_pathway_name in names(pathways)){ 
    local_pathway <- pathways[[local_pathway_name]]
    pathway_network <- network[intersect(which(network$Gene1 %in% local_pathway), which(network$Gene2 %in% local_pathway)),]
    genes_in_network <- unique(c(as.character(network$Gene1[which(network$Gene1 %in% local_pathway)]), as.character(network$Gene2[which(network$Gene2 %in% local_pathway)])))
    percentage_genes_in_network <- sum(genes_in_network %in% local_pathway)/length(local_pathway)*100
    if (nrow(pathway_network) == 0){
      pathway_components <- rbind(pathway_components, c(local_pathway_name, 0, 0, percentage_genes_in_network))
    }else{
      vertices <- sort(unique(c(as.character(pathway_network$Gene1), as.character(pathway_network$Gene2))))
      pathway_network_NEL <- graphNEL(nodes=vertices)
      pathway_network_NEL <- addEdge(as.character(pathway_network$Gene1), as.character(pathway_network$Gene2), pathway_network_NEL, 1)
      pathway_components <- rbind(pathway_components, c(local_pathway_name, length(connComp(pathway_network_NEL)), paste(sapply(connComp(pathway_network_NEL), length), collapse=','), percentage_genes_in_network))
    }
  }
  colnames(pathway_components) <- c("Pathway","Number_connected_components", "Size_of_connected_components","Percentage_pathway_genes_in_network")
  pathway_components <- as.data.frame(pathway_components)
  return(pathway_components)
}

obtain_overlap_between_pathways <- function(){
  number_overlap_between_pathways <- sapply(names(pathways), function(local){
    genes_in_local <- pathways[[local]]
    overlap_local <- sapply(names(pathways), identify_overlap, genes_in_local)
  })
  
  number_overlap_between_pathways <- as.data.frame(number_overlap_between_pathways)
  
  rowCol <- expand.grid(rownames(number_overlap_between_pathways),
                        colnames(number_overlap_between_pathways)) 
  labs <- rowCol[as.vector(upper.tri(number_overlap_between_pathways,diag=F)),]
  overlap_between_pathways_df <- cbind(labs, number_overlap_between_pathways[upper.tri(number_overlap_between_pathways,diag=F)])
  
  colnames(overlap_between_pathways_df) <- c("Pathway1", "Pathway2", "Number_genes_overlap")
 
  return(overlap_between_pathways_df)
}

load_ssGSEA_object <- function(pathway_database, paired){
  if (pathway_database %in% c("KEGG", "Reactome")){
    signature_ssGSEA <- 
    ssGSEA_results <- list()
    
    if(paired==TRUE){
      ssGSEA_results #List, with names as tumour types, and items as the results matrix of ssGSEA calculated with GSVA per sample (tumour and normal paired samples)
    }else{
      ssGSEA_results #List, with names as tumour types, and items as the results matrix of ssGSEA calculated with GSVA per sample (tumour and normal)
    }
      ssGSEA_results_temp <- ssGSEA_results
      ssGSEA_results <- list()
    for (tissue_type in tumours){
      ssGSEA_results[["normal"]][[tissue_type]] <- ssGSEA_results_temp[[tissue_type]][,which(substr(colnames(ssGSEA_results_temp[[tissue_type]]), 6,6) == "N")]
      ssGSEA_results[["tumour"]][[tissue_type]] <- ssGSEA_results_temp[[tissue_type]][,which(substr(colnames(ssGSEA_results_temp[[tissue_type]]), 6,6) == "T")]
    }
    return(ssGSEA_results)
  }
}

calculate_correlation_of_ssGSEA_scores <- function(ssGSEA_results, tumours){
  cor_ssGSEA <- list()
  for(type in c("normal", "tumour")){
    for(tumour in tumours){
      cor_ssGSEA[[tumour]][[type]] <- as.data.frame(cor(t(ssGSEA_results[[type]][[tumour]]), method="sp"))
      cor_ssGSEA[[tumour]][[type]] <- as.data.frame(cor(t(ssGSEA_results[[type]][[tumour]]), method="sp"))
      
      rowCol <- expand.grid(rownames(cor_ssGSEA[[tumour]][[type]]),
                            colnames(cor_ssGSEA[[tumour]][[type]])) 
      labs <- rowCol[as.vector(upper.tri(cor_ssGSEA[[tumour]][[type]],diag=F)),]
      cor_ssGSEA[[tumour]][[type]] <- cbind(labs, cor_ssGSEA[[tumour]][[type]][upper.tri(cor_ssGSEA[[tumour]][[type]],diag=F)])
      
      colnames(cor_ssGSEA[[tumour]][[type]]) <- c("Pathway1", "Pathway2", "Correlation")
      t <- cor_ssGSEA[[tumour]][[type]][,c(2,1,3)]
      colnames(t) <- c("Pathway1", "Pathway2", "Correlation")
      cor_ssGSEA[[tumour]][[type]] <- rbind(cor_ssGSEA[[tumour]][[type]],
                                            t)
      
      cor_ssGSEA[[tumour]][[type]] <- as.data.frame(cor_ssGSEA[[tumour]][[type]])      
    }
  }

  return(cor_ssGSEA)
}

calculate_correlation_of_ssGSEA_scores_sample_subset <- function(ssGSEA_samples){
      correlations <- as.data.frame(cor(t(ssGSEA_samples), method="sp"))
      rowCol <- expand.grid(rownames(correlations),
                            colnames(correlations)) 
      labs <- rowCol[as.vector(upper.tri(correlations,diag=F)),]
      correlations <- cbind(labs, correlations[upper.tri(correlations,diag=F)])
      
      colnames(correlations) <- c("Pathway1", "Pathway2", "Correlation")
      t <- correlations[,c(2,1,3)]
      colnames(t) <- c("Pathway1", "Pathway2", "Correlation")
      correlations <- rbind(correlations,t)
      
      correlations <- as.data.frame(correlations)      
  
  return(correlations)
}


#excluding from the calculation the edges that appear within any of the pathways
find_number_edges_connect_pathways_no_overlap <- function(pathway_database, network_database, network_pathways, network_v2){
  load(paste(network_database, "_shortest_paths_matrix_entrez.Rdata", sep=""))  #matrix of shortest paths obtained form igraph
  
  connection_distance_1_between_pathways <- lapply(names(pathways), function(local_pathway){
    genes_in_local_pathway1 <- pathways[[local_pathway]]
    edges_connecting_pathways <- lapply(names(pathways), function(x) {find_mode_of_distance_1v2(x,genes_in_local_pathway1, shortest_paths_matrix)})
 
    print(local_pathway)
    return(edges_connecting_pathways)
  })
  
  names(connection_distance_1_between_pathways) <- names(pathways)
  connection_distance_1_between_pathways <- lapply(names(connection_distance_1_between_pathways), function(p1){
    names(connection_distance_1_between_pathways[[p1]]) <- names(pathways)
    return(connection_distance_1_between_pathways[[p1]])
  })
  names(connection_distance_1_between_pathways) <- names(pathways)
  
  #Find edges that join pathways and that are within pathways
  number_edges <- do.call('rbind', lapply(names(connection_distance_1_between_pathways), function(p1_name){
    print(p1_name)
    # p1_edges <- connection_distance_1_between_pathways[[p1_name]]
    # do.call('rbind', lapply(names(p1_edges), function(p2_name){
    #   if(ncol(p1_edges[[p2_name]]) == 2 && nrow(p1_edges[[p2_name]])>0){
    #     p2_edges <- graph.data.frame(p1_edges[[p2_name]],directed=FALSE)  
    #     return(c(p1_name, p2_name, gsize(p2_edges), gsize(p2_edges %s% network_pathways[[p1_name]]),
    #              gsize(p2_edges %s% network_pathways[[p2_name]])))        
    #   } else{
    #     return(c(p1_name, p2_name, 0, 0, 0))
    #   }
    
    edges_between_p1 <- connection_distance_1_between_pathways[[p1_name]] #same if p2_name
   do.call('rbind', lapply(names(edges_between_p1), function(p2_name){
      if(ncol(edges_between_p1[[p2_name]]) == 2 && nrow(edges_between_p1[[p2_name]])>0){
        edges_between_graph <- graph.data.frame(edges_between_p1[[p2_name]], directed=FALSE)
        edges_between_graph <- simplify(edges_between_graph, remove.loops=TRUE, remove.multiple=TRUE)
        
        p1_edges <- network_pathways[[p1_name]]
        p2_edges <- network_pathways[[p2_name]]
        #Total number of edges between processes, internal edges of p1 minus the internal edges it shares with p2, internal edges of p2 (including those shared with p1)
        return(c(p1_name, p2_name, gsize(edges_between_graph), gsize(p1_edges %s% edges_between_graph) - gsize(p1_edges %s% p2_edges),
                 gsize(p2_edges %s% edges_between_graph)))        
      } else{
        return(c(p1_name, p2_name, 0, 0, 0))
      }
#      if(gsize(p2_edges) != 0){
#         return(c(p1_name, p2_name, gsize(p2_edges), gsize(p2_edges %s% network_pathways[[p1_name]]),
#                  gsize(p2_edges %s% network_pathways[[p2_name]])))        
#       } else{
#         return(c(p1_name, p2_name, 0, 0, 0))
#       }

    }))
  }))
  colnames(number_edges) <- c("Pathway1", "Pathway2", "Number_edges_distance_1", "Number_overlaping_edges_pathway_1", "Number_overlaping_edges_pathway_2")
  
  return(number_edges)
}

load_clinical_data <- function(tumour){
  clinical_data_TCGA <- read.delim(paste("nationwidechildrens.org_clinical_patient_", tumour, ".txt", sep=""))  #Metadata from TCGA's website
  
  colnames(clinical_data_TCGA) <- unlist(clinical_data_TCGA[1,])
  clinical_data_TCGA <- clinical_data_TCGA[-c(1,2),]
  patients_clinical <- data.frame(do.call('rbind', strsplit(as.character(clinical_data_TCGA$bcr_patient_barcode),'-',fixed=TRUE)))
  patients_clinical <- paste(patients_clinical$X3, "_T", sep="")
  
  clinical_data <- data.frame(pT = clinical_data_TCGA$pathologic_T, pN = clinical_data_TCGA$pathologic_N, 
                              pM = clinical_data_TCGA$pathologic_M, pS = clinical_data_TCGA$pathologic_stage)
  
  if ("gleason_score" %in% colnames(clinical_data_TCGA)){
    clinical_data <- data.frame(clinical_data, Gleason=clinical_data_TCGA$gleason_score)
  }
  if ("neoplasm_histologic_grade" %in%  colnames(clinical_data_TCGA)){
    clinical_data <- data.frame(clinical_data, Grade=clinical_data_TCGA$neoplasm_histologic_grade)
  }
  if("lab_proc_her2_neu_immunohistochemistry_receptor_status" %in% colnames(clinical_data_TCGA)){
    clinical_data <- data.frame(clinical_data, her2_status=clinical_data_TCGA$lab_proc_her2_neu_immunohistochemistry_receptor_status)    
  }
  
  clinical_data <- droplevels(clinical_data)
  clinical_data$Patient <- patients_clinical
  
  levels(clinical_data$pT) <- unique(c(levels(clinical_data$pT), "T1", "T2", "T3", "T4"))
  clinical_data$pT <- replace(clinical_data$pT, which(clinical_data$pT %in% c("T1", "T1a", "T1b", "T1c")), "T1")
  clinical_data$pT <- replace(clinical_data$pT, which(clinical_data$pT %in% c("T2", "T2a", "T2b", "T2c")), "T2")
  clinical_data$pT <- replace(clinical_data$pT, which(clinical_data$pT %in% c("T3", "T3a", "T3b", "T3c")), "T3")
  clinical_data$pT <- replace(clinical_data$pT, which(clinical_data$pT %in% c("T4", "T4a", "T4b", "T4d")), "T4")
  
  levels(clinical_data$pN) <- unique(c(levels(clinical_data$pN), "N1", "N2", "N3", "N0"))
  clinical_data$pN <- replace(clinical_data$pN, which(clinical_data$pN %in% c("N0", "N0 (i-)", "N0 (i+)", "N0 (mol+)")), "N0")
  clinical_data$pN <- replace(clinical_data$pN, which(clinical_data$pN %in% c("N1", "N1a", "N1b", "N1c", "N1mi")), "N1")
  clinical_data$pN <- replace(clinical_data$pN, which(clinical_data$pN %in% c("N2", "N2a", "N2b", "N2c")), "N2")
  clinical_data$pN <- replace(clinical_data$pN, which(clinical_data$pN %in% c("N3", "N3a", "N3b", "N3c")), "N3")
  
  levels(clinical_data$pM) <- unique(c(levels(clinical_data$pM), "M1", "M0"))
  clinical_data$pM <- replace(clinical_data$pM, which(clinical_data$pM %in% c("M0", "cM0 (i+)")), "M0")
  clinical_data$pM <- replace(clinical_data$pM, which(clinical_data$pM %in% c("M1", "M1a", "M1b", "M1c")), "M1")
  
  levels(clinical_data$pS) <- unique(c(levels(clinical_data$pS), "StageI", "StageII", "StageIII", "StageIV"))
  clinical_data$pS <- replace(clinical_data$pS, which(clinical_data$pS %in% c("Stage I", "Stage IA", "Stage IB")), "StageI")
  clinical_data$pS <- replace(clinical_data$pS, which(clinical_data$pS %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")), "StageII")
  clinical_data$pS <- replace(clinical_data$pS, which(clinical_data$pS %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")), "StageIII")
  clinical_data$pS <- replace(clinical_data$pS, which(clinical_data$pS %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC")), "StageIV")
  
  if (tumour == "PRAD"){
    clinical_data$Gleason <- as.character(clinical_data$Gleason)
    clinical_data$Gleason <- sapply(clinical_data$Gleason, function(gle){
      if (as.numeric(gle) < 10){
        return(paste(0, gle, sep=""))
      }else
        return(gle)
    })
    clinical_data$Gleason <- paste("Gleason_", clinical_data$Gleason, sep="")

    clinical_data$Gleason <- as.factor(clinical_data$Gleason)
  }

  
  clinical_data <- droplevels(clinical_data)
  return(clinical_data)
}

calculate_correlation_of_ssGSEA_scores_v2 <- function(local_ssGSEA_results){
  result <- as.data.frame(cor(t(local_ssGSEA_results), method="sp"))
  rowCol <- expand.grid(rownames(result),
                        colnames(result)) 
  labs <- rowCol[as.vector(upper.tri(result,diag=F)),]
  result <- cbind(labs, result[upper.tri(result,diag=F)])
  
  colnames(result) <- c("Pathway1", "Pathway2", "Correlation")
  
  result <- as.data.frame(result)      
  return(result)
}

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

plot_circle_of_correlations <- function(edge_table, tumour, type, all_nodes, cutoff, expression, network_type="multiple", layout="circular", center=NULL){
  #  edge_table <- edge_table[which(edge_table[,3]>0.01)]
  if (nrow(edge_table) > 0){
    #all_nodes_age=as.character(pathway_ages[match(all_nodes, pathway_ages[,1]),2])
    
    #edge_table <- edge_table[order(all_nodes_age),]
    
    weight_scaled <- abs(as.numeric(edge_table[,4]))
    edge_table <- as.data.frame(edge_table)   
    edge_table$color <- NA
    if(network_type == "multiple"){
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="MC_MC", rgb(red=0, green=0, blue=1, alpha=weight_scaled[edge_table[,5]=="MC_MC"]))
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_UC", rgb(red=1, green=0, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_UC"]))
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_MC", rgb(red=0, green=1, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_MC"]))
      
    }else{
      edge_table$color <- grey.colors(n=length(weight_scaled), alpha=weight_scaled)
    }
    edge_table$lty <- sapply(edge_table[,4], function(x){ 
      if(x>0){
        return(1) 
      } else{
        return(2)
      } })
    
    graph <- graph.data.frame(edge_table, directed=FALSE)
    graph <- simplify(graph, remove.loops=TRUE, remove.multiple=TRUE, edge.attr.com="first")
    graph <- add.vertices(graph, nv = length(all_nodes[!(all_nodes %in% V(graph)$name)]), name=all_nodes[!(all_nodes %in% V(graph)$name)])
    #graph <- subgraph.edges(graph, eids=which(E(graph)$"Degree_of_interaction" > cutoff), delete.vertices=FALSE)
    
    nodes <- V(graph)$name
    
    nodes_ages <- as.character(pathway_ages[match(nodes, pathway_ages[,1]),2])
    node_shapes <- replace(nodes_ages, nodes_ages=="UC", "circle")
    node_shapes <- replace(node_shapes, node_shapes=="MC", "square")
    node_shapes <- replace(node_shapes, is.na(node_shapes), "rectangle")
    if(expression == "qusage"){
      if(type %in% c("normal")){
        nodes_fill <- rep("white", length(nodes))
      } else{
        qusage_values <- load_qusage_values() 
        if (tumour == "tumour-median"){
          qusage_values <- apply(qusage_values,2,median, na.rm=T)
        }else{
          qusage_values <- qusage_values[which(rownames(qusage_values) == tumour),]
        }
        qusage_values <- qusage_values[names(qusage_values) %in% nodes]
        #number_vector <- round(seq(-1.61, 0.8, 0.01),2)
        
        #     color_vector_1 <- colorpanel(n=sum(number_vector<0), low="blue", high="white")
        #     color_vector_2 <- colorpanel(n=sum(number_vector>0), low="white", high="red")
        #     color_vector <- c(color_vector_1, "white", color_vector_2)
        #     nodes_fill <- round(unname(qusage_values[match(nodes, names(qusage_values))]),2)
        #     nodes_fill <- color_vector[match(nodes_fill,number_vector)]
        #     nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
        nodes_fill <- qusage_values[names(qusage_values) %in% nodes]
        nodes_fill <- ifelse(nodes_fill>0, "red", "blue")
        nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
      }
      
    } else if (expression=="ssGSEA"){
      #Load ssGSEA_results. List with tumours as names, and matrix of ssGSEA scores obtained with GSVA as items
      if (tumour == "tumour-median"){
        ssGSEA_results_df <- vector()
        for (tumour_type in tumours){
          ssGSEA_results_df <- cbind(ssGSEA_results_df,
                                     ssGSEA_results[[tumour_type]])
        }
        if (type == "normal"){
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="N")]
                } else{
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="T")]
        }
        ssGSEA_results_median <- apply(ssGSEA_results_df2, 1, median)
        nodes_fill <- ssGSEA_results_median[names(ssGSEA_results_median) %in% nodes]

        color_vector <- colorpanel(n=length(nodes_fill), low="white", high="purple")
        nodes_fill_names <- names(sort(nodes_fill))
        nodes_fill <- color_vector
        names(nodes_fill) <- nodes_fill_names
        nodes_fill <- nodes_fill[match(nodes, names(nodes_fill))]        
      }
    }
    
    E(graph)$weight <- 2.5
    
    V(graph)$color <- unname(nodes_fill)[match(V(graph)$name, names(nodes_fill))]
    V(graph)$vertex.shape <- node_shapes
    
    nodes_with_edges <- unique(as.vector(get.edgelist(graph)))
    V(graph)$size <- ifelse(V(graph)$name %in% nodes_with_edges, 7, 3)
    
    if(layout == "circular"){
      coords <- layout_in_circle(graph)
      lab.locs <- radian.rescale(x=1:length(V(graph)$name), direction=-1, start=0)
    } else if (layout =="star"){
      coords <- layout.star(graph, center=center)
    }
    
    V(graph)$name <- sapply(strsplit(V(graph)$name, " "), function(x){
      if(length(x) != 1){
        paste(x, "", collapse="\n")        
      } else{
        return(paste(x, "", collapse=" "))
      }
    })
   if(network_type == "multiple" && layout=="circular"){
     
      par(bg = "black")
      V(graph)$label.cex = .55  #.55
      p <- plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="white", 
                  vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                  main=paste(tumour, network_database, sep="\n"),
                  vertex.shape= V(graph)$vertex.shape,
                  vertex.frame.color="white",
                  vertex.label.dist=.45,
                  vertex.label.degree=lab.locs)   
      #title(paste(tumour, network_database, type, sep="\n"), col.main = "white")
    } else if (network_type == "individual" && layout=="circle"){
      par(bg="white")
      V(graph)$label.cex = 1
      plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="black", 
                  vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                  vertex.shape= V(graph)$vertex.shape,
                  vertex.frame.color="white", 
                  vertex.label.dist=.45,
                  vertex.label.degree=lab.locs,
                  vertex.size=20)     
      title(tumour, col.main = "black")
    } else if (network_type == "individual" && layout=="star"){
      par(bg="white")
      V(graph)$label.cex =1
      plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="black", 
                  vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                  vertex.shape= V(graph)$vertex.shape,
                  vertex.frame.color="white", 
                  vertex.label.dist=.45,
                  vertex.size=20)     
      title(tumour, col.main = "black")
      
    }
  } else{
    par(bg = "black")
    g <- barabasi.game(10)
    plot(g, vertex.label.color="black", vertex.color="black", edge.color="black")
  }
}

plot_circle_of_differences <- function(edge_table, tumour, type, all_nodes, cutoff, expression){
  #  edge_table <- edge_table[which(edge_table[,3]>0.01)]
  if (nrow(edge_table) > 0){
    weight_scaled <- abs(as.numeric(edge_table[,4]))
    edge_table <- as.data.frame(edge_table)   
    edge_table$color <- NA
    edge_table$color <- replace(edge_table$color, edge_table[,5]=="MC_MC", rgb(red=0, green=0, blue=1, alpha=weight_scaled[edge_table[,5]=="MC_MC"]))
    edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_UC", rgb(red=1, green=0, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_UC"]))
    edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_MC", rgb(red=0, green=1, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_MC"]))
    edge_table$lty <- sapply(edge_table$Correlation, function(x){ 
      if(x>0){
        return(1) 
      } else{
        return(2)
      } })
    
    graph <- graph.data.frame(edge_table, directed=FALSE)
    graph <- simplify(graph, remove.loops=TRUE, remove.multiple=TRUE, edge.attr.com="first")
    graph <- add.vertices(graph, nv = length(all_nodes[!(all_nodes %in% V(graph)$name)]), name=all_nodes[!(all_nodes %in% V(graph)$name)])
    #graph <- subgraph.edges(graph, eids=which(E(graph)$"Degree_of_interaction" > cutoff), delete.vertices=FALSE)
    
    nodes <- V(graph)$name
    
    nodes_ages <- as.character(pathway_ages[match(nodes, pathway_ages[,1]),2])
    node_shapes <- replace(nodes_ages, nodes_ages=="UC", "circle")
    node_shapes <- replace(node_shapes, node_shapes=="MC", "square")
    node_shapes <- replace(node_shapes, is.na(node_shapes), "rectangle")
    if(expression == "qusage"){
      if(type == "normal"){
        nodes_fill <- rep("white", length(nodes))
      } else{
        qusage_values <- load_qusage_values() 
        if (tumour == "tumour-median"){
          qusage_values <- apply(qusage_values,2,median, na.rm=T)
        }else{
          qusage_values <- qusage_values[which(rownames(qusage_values) == tumour),]
        }
        qusage_values <- qusage_values[names(qusage_values) %in% nodes]
        #number_vector <- round(seq(-1.61, 0.8, 0.01),2)
        
        #     color_vector_1 <- colorpanel(n=sum(number_vector<0), low="blue", high="white")
        #     color_vector_2 <- colorpanel(n=sum(number_vector>0), low="white", high="red")
        #     color_vector <- c(color_vector_1, "white", color_vector_2)
        #     nodes_fill <- round(unname(qusage_values[match(nodes, names(qusage_values))]),2)
        #     nodes_fill <- color_vector[match(nodes_fill,number_vector)]
        #     nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
        nodes_fill <- qusage_values[names(qusage_values) %in% nodes]
        nodes_fill <- ifelse(nodes_fill>0, "red", "blue")
        nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
      }
      
    } else if (expression=="ssGSEA"){
      #Load ssGSEA_results. List with tumours as names, and matrix of ssGSEA scores obtained with GSVA as items
      
      if (tumour == "tumour-median"){
        ssGSEA_results_df <- vector()
        for (tumour_type in tumours){
          ssGSEA_results_df <- cbind(ssGSEA_results_df,
                                     ssGSEA_results[[tumour_type]])
        }
        if (type == "normal"){
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="N")]
         } else{
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="T")]
         }
        ssGSEA_results_median <- apply(ssGSEA_results_df2, 1, median)
         nodes_fill <- ssGSEA_results_median[names(ssGSEA_results_median) %in% nodes]

        color_vector <- colorpanel(n=length(nodes_fill), low="white", high="purple")
        nodes_fill_names <- names(sort(nodes_fill))
        nodes_fill <- color_vector
        names(nodes_fill) <- nodes_fill_names
        nodes_fill <- nodes_fill[match(nodes, names(nodes_fill))]        
      }
    }
    
    E(graph)$weight <- 2.5
    
    V(graph)$color <- unname(nodes_fill)[match(V(graph)$name, names(nodes_fill))]
    V(graph)$vertex.shape <- node_shapes
    
    nodes_with_edges <- unique(as.vector(get.edgelist(graph)))
    V(graph)$size <- ifelse(V(graph)$name %in% nodes_with_edges, 7, 3)
    V(graph)$label.cex =.6
    par(bg = "black")
    
    coords <- layout_in_circle(graph, order = order(nodes_ages, V(graph)$name))
    
    V(graph)$name <- sapply(strsplit(V(graph)$name, " "), function(x){
      if(length(x) <= 2){
        paste(x, "", collapse="\n")        
      } else{
        return(paste(x, "", collapse=" "))
      }
    })
    
    plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="white", 
                vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                main=paste(tumour, network_database, sep="\n"),
                vertex.shape= V(graph)$vertex.shape,
                vertex.frame.color="white")
    title(paste(tumour, network_database, type, sep="\n"), col.main = "white")
    
  } else{
    par(bg = "black")
    g <- barabasi.game(10)
    plot(g, vertex.label.color="black", vertex.color="black", edge.color="black")
  }
}




polar.layout <- function(radii, angles) {
  cbind(radii*cos(angles), -radii*sin(angles))        
}

load_qusage_values <- function(){
  tumours <- c("LUAD", "LUSC", "BRCA","PRAD","LIHC", "COAD","STAD")
  
  qusage_results <- list()
  qusage_results_matrix <- vector()
  for (type in tumours){
    qusage_results[[type]] <- #Load matrix with qusage results. Columns: pathway.name,log.fold.change,p.Value,FDR
    qusage_results_names <- qusage_results[[type]][,1]
    qusage_results[[type]] <- t(apply(qusage_results[[type]],1,function(local){
      if (is.na(local[4]) || as.numeric(local[4]) >= 0.05){
        return(rep(NA,4))
      } else{
        return(local)
      }
    }))
    qusage_results_matrix <- cbind(qusage_results_matrix, qusage_results[[type]][,2])
  }
  qusage_results_matrix <- apply(qusage_results_matrix,1,as.numeric)
  gonames <- Term(GOTERM)
  colnames(qusage_results_matrix) <- qusage_results_names
  rownames(qusage_results_matrix) <- tumours
  return(qusage_results_matrix)
}


calculate_degree_interaction <- function(network){
  genes_in_network <- unique(c(network[,1], network[,2]))
  network_graph <- graph.data.frame(network, directed=FALSE)
  network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
  
  network_v2 <- get.data.frame(network_graph)
  colnames(network_v2) <- c("Gene1", "Gene2")
  #Create networks for each pathway
  network_pathways <- lapply(names(pathways), function(name_local_pathway){
    genes_local_pathway <- pathways[[name_local_pathway]]
    network_local <- network_v2[intersect(which(network_v2[,1] %in% genes_local_pathway),
                                          which(network_v2[,2] %in% genes_local_pathway)),]
    #return(graph.data.frame(network_local, directed=FALSE))
    return(induced.subgraph(network_graph, which(V(network_graph)$name %in% genes_local_pathway)))
  }
  )
  names(network_pathways) <- names(pathways)
  
  connection_distance_1_between_pathways2 <- find_number_edges_connect_pathways_no_overlap(pathway_database, network_database, network_pathways, network_v2)
  
  #Counting only genes that are not shared between the groups but that are in the network
  connection_distance_1_between_pathways2 <- cbind(connection_distance_1_between_pathways2, t(apply(connection_distance_1_between_pathways2, 1, function(row){
    path1 <- as.character(row[1])
    path2 <- as.character(row[2])
    path1_v <- pathways[[path1]]  #V(network_pathways[[path1]])$name
    path2_v <- pathways[[path2]]  #V(network_pathways[[path2]])$name
    path1_v <- path1_v[path1_v %in% genes_in_network]
    path2_v <- path2_v[path2_v %in% genes_in_network]
    return(c(length(path1_v) - sum(path1_v %in% path2_v),
            length(path2_v) - sum(path2_v %in% path1_v)))
    })))
    
  #Using the total number of nodes of each group
  #connection_distance_1_between_pathways2 <- cbind(connection_distance_1_between_pathways2, 
  #                                                 as.numeric(as.character(sapply(as.character(connection_distance_1_between_pathways2[,1]), function(path){
  #                                                   vcount(network_pathways[[path]])
  #                                                 }))),
  #                                                 as.numeric(as.character(sapply(as.character(connection_distance_1_between_pathways2[,2]), function(path){
  #                                                   vcount(network_pathways[[path]])
  #                                                 }))))
  colnames(connection_distance_1_between_pathways2)[6:7] <- c("Number_genes_pathway_1", "Number_genes_pathway_2")
  
  number_of_valid_connections <- as.numeric(connection_distance_1_between_pathways2[,3])-as.numeric(connection_distance_1_between_pathways2[,4])-as.numeric(connection_distance_1_between_pathways2[,5])
  total_number_possible_connections <- as.numeric(connection_distance_1_between_pathways2[,6])*as.numeric(connection_distance_1_between_pathways2[,7])
  connection_distance_1_between_pathways2 <- cbind(connection_distance_1_between_pathways2,
                                                   number_of_valid_connections/total_number_possible_connections)
  colnames(connection_distance_1_between_pathways2)[8] <- "Degree_of_interaction"
  
  connection_distance_1_between_pathways2 <- do.call('rbind', apply(as.data.frame(connection_distance_1_between_pathways2), 1, function(row){
    if(row[1] != row[2]){
      return(row)
    }
  }))
  
#  connection_distance_1_between_pathways2 <- connection_distance_1_between_pathways2[match(interaction(cor_ssGSEA[[tumour]][[type]][,1], cor_ssGSEA[[tumour]][[type]][,2]),
#                                                                                           interaction(connection_distance_1_between_pathways2[,1], connection_distance_1_between_pathways2[,2])),]
  return(connection_distance_1_between_pathways2)
}

add_ages <- function(dataframe){
  dataframe <- as.data.frame(dataframe)
  pathway_ages <- #Data frame of GOlims and Ages
  
  dataframe <- data.frame(dataframe, Ages1=pathway_ages[match(dataframe$Pathway1, pathway_ages[,1]),2],
                     Ages2=pathway_ages[match(dataframe$Pathway2, pathway_ages[,1]),2])
  dataframe <- data.frame(dataframe, Interaction_type=paste(dataframe$Ages1, dataframe$Ages2, sep="_"))
  dataframe$Interaction_type <- replace(dataframe$Interaction_type, which(dataframe$Interaction_type == "MC_UC"), "UC_MC")
  dataframe$Interaction_type <- factor(dataframe$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))
  return(dataframe)
}

intersect_edge_type_with_cutoff <- function(df, edge_type){
  if (edge_type != "all"){
    return(intersect(which(df$Interaction_type == edge_type),
                     intersect(which(df$Degree_of_interaction >= cutoff[[edge_type]]),
                               which(!is.infinite(df$Degree_of_interaction)))))
  } else{
    a <- vector()
    for (edge_type in c("UC_UC", "UC_MC", "MC_MC")){
      a <- c(a, intersect(which(df$Interaction_type == edge_type),
                          intersect(which(df$Degree_of_interaction >= cutoff[[edge_type]]),
                           which(!is.infinite(df$Degree_of_interaction)))))
    }
    return(sort(a))
  }
}

build_cor_number_edges_matrix <- function(cor_number_edges, correlation_type){
  cor_number_edges_matrix <- vector()
  for(tumour in tumours){
    cor_number_edges_matrix <- cbind(cor_number_edges_matrix,
                                       cor_number_edges[[tumour]][[correlation_type]])
  }
  rownames(cor_number_edges_matrix) <- as.character(interaction(cor_number_edges[[tumour]][,1],
                                                              cor_number_edges[[tumour]][,2]))
  colnames(cor_number_edges_matrix) <- tumours
  return(cor_number_edges_matrix)
}

plot_heatmap_of_interconnectedness <- function(t){
  t$Degree_of_interaction <- as.numeric(as.character(t$Degree_of_interaction))
  temp <- t[,c(2,1,3:ncol(t))]
  colnames(temp)[1:2] <- c("Pathway1", "Pathway2")
  t2 <- rbind(t, temp)
  
  t3 <- dcast(t2[,c("Pathway1", "Pathway2", "Degree_of_interaction")], Pathway1 ~ Pathway2)
  rownames(t3) <- t3[,1]
  t3 <- t3[,-1]
  
  #Replace Inf with 0
  t32 <-apply(t3,2,function(x){
    replace(x, is.infinite(x), 0)
  })
  
  # t32 <-apply(t3,2,function(x){
  #   replace(x, is.nan(x), 0)
  # })
  
  h32 <- hclust(dist(t32))
  #plot(h32, cex=.7)
  
  ordering <- h32$labels[h32$order]
  t2$Pathway1 <- factor(t2$Pathway1, levels=ordering)
  t2$Pathway2 <- factor(t2$Pathway2, levels=ordering)
  # t2$Color1 <- ifelse(t2$Ages1 == "UC", "red",
  #                     ifelse(t2$Ages1 == "MC", "blue", "black"))
  # t2$Color2 <- ifelse(t2$Ages2 == "UC", "red",
  #                     ifelse(t2$Ages2 == "MC", "blue", "black"))
  # t2$Color1[is.na(t2$Color1)] <- "black"
  # t2$Color2[is.na(t2$Color2)] <- "black"
  
  ordering <- ordering[ordering != "transposition"]
  ages_ordering <- pathway_ages[match(ordering, pathway_ages[,1]), 2]
  ages_ordering_color <- ifelse(ages_ordering == "UC", "red",
                                ifelse(ages_ordering == "MC", "blue", "black"))
  ages_ordering_color[is.na(ages_ordering_color)] <- "black"
  
  t2 <- t2[t2[,1] != "transposition",]
  t2 <- t2[t2[,2] != "transposition",]
  
  t2$Degree_of_interaction <- as.numeric(as.character(t2$Degree_of_interaction))
  t2$Degree_of_interaction[is.infinite(t2$Degree_of_interaction)] <- 0
  g <- ggplot(t2, aes(Pathway1, Pathway2))+
    geom_tile(aes(fill = Degree_of_interaction),colour = "white") + 
    scale_fill_gradient(low = "white",high = "blue", guide_legend(title="Interconnectedness"))+
    xlab("")+
    ylab("")+
    theme(axis.text.x = element_text(angle = 90, hjust=1, colour=ages_ordering_color, size =7),
          axis.text.y = element_text(colour=ages_ordering_color, size = 7))
  print(g)
}


plot_circle_of_correlations_PNAS <- function(edge_table, tumour, type, all_nodes, cutoff, expression, network_type="multiple", layout="circular", center=NULL){
  #  edge_table <- edge_table[which(edge_table[,3]>0.01)]
  if (nrow(edge_table) > 0){
    #all_nodes_age=as.character(pathway_ages[match(all_nodes, pathway_ages[,1]),2])
    
    #edge_table <- edge_table[order(all_nodes_age),]
    
    weight_scaled <- abs(as.numeric(edge_table[,4]))
    edge_table <- as.data.frame(edge_table)   
    edge_table$color <- NA
    if(network_type == "multiple"){
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="MC_MC", rgb(red=0, green=0, blue=1, alpha=weight_scaled[edge_table[,5]=="MC_MC"]))
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_UC", rgb(red=1, green=0, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_UC"]))
      edge_table$color <- replace(edge_table$color, edge_table[,5]=="UC_MC", rgb(red=0, green=1, blue=0, alpha=weight_scaled[edge_table[,5]=="UC_MC"]))
      
    }else{
      edge_table$color <- grey.colors(n=length(weight_scaled), alpha=weight_scaled)
    }
    edge_table$lty <- sapply(edge_table[,4], function(x){ 
      if(x>0){
        return(1) 
      } else{
        return(2)
      } })
    
    graph <- graph.data.frame(edge_table, directed=FALSE)
    graph <- simplify(graph, remove.loops=TRUE, remove.multiple=TRUE, edge.attr.com="first")
    graph <- add.vertices(graph, nv = length(all_nodes[!(all_nodes %in% V(graph)$name)]), name=all_nodes[!(all_nodes %in% V(graph)$name)])
    #graph <- subgraph.edges(graph, eids=which(E(graph)$"Degree_of_interaction" > cutoff), delete.vertices=FALSE)
    
    nodes <- V(graph)$name
    
    nodes_ages <- as.character(pathway_ages[match(nodes, pathway_ages[,1]),2])
    node_shapes <- replace(nodes_ages, nodes_ages=="UC", "circle")
    node_shapes <- replace(node_shapes, node_shapes=="MC", "square")
    node_shapes <- replace(node_shapes, is.na(node_shapes), "rectangle")
    if(expression == "qusage"){
      if(type %in% c("normal")){
        nodes_fill <- rep("white", length(nodes))
      } else{
        qusage_values <- load_qusage_values() 
        if (tumour == "tumour-median"){
          qusage_values <- apply(qusage_values,2,median, na.rm=T)
        }else{
          qusage_values <- qusage_values[which(rownames(qusage_values) == tumour),]
        }
        qusage_values <- qusage_values[names(qusage_values) %in% nodes]
        #number_vector <- round(seq(-1.61, 0.8, 0.01),2)
        
        #     color_vector_1 <- colorpanel(n=sum(number_vector<0), low="blue", high="white")
        #     color_vector_2 <- colorpanel(n=sum(number_vector>0), low="white", high="red")
        #     color_vector <- c(color_vector_1, "white", color_vector_2)
        #     nodes_fill <- round(unname(qusage_values[match(nodes, names(qusage_values))]),2)
        #     nodes_fill <- color_vector[match(nodes_fill,number_vector)]
        #     nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
        nodes_fill <- qusage_values[names(qusage_values) %in% nodes]
        nodes_fill <- ifelse(nodes_fill>0, "red", "blue")
        nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")
      }
      
    } else if (expression=="ssGSEA"){
      #Load ssGSEA_results
      if (tumour == "tumour-median"){
        ssGSEA_results_df <- vector()
        for (tumour_type in tumours){
          ssGSEA_results_df <- cbind(ssGSEA_results_df,
                                     ssGSEA_results[[tumour_type]])
        }
        if (type == "normal"){
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="N")]
          #other <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="T")]
        } else{
          ssGSEA_results_df2 <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="T")]
          #other <- ssGSEA_results_df[,which(substr(colnames(ssGSEA_results_df), 6,6)=="N")]
        }
        ssGSEA_results_median <- apply(ssGSEA_results_df2, 1, median)
        #other_median <- apply(other, 1, median)
        #ssGSEA_results_median <- ssGSEA_results_median/other_median
        nodes_fill <- ssGSEA_results_median[names(ssGSEA_results_median) %in% nodes]
        #nodes_fill <- ifelse(nodes_fill>1, "red", "blue")
        #nodes_fill <- replace(nodes_fill, is.na(nodes_fill), "white")        
        
        color_vector <- colorpanel(n=length(nodes_fill), low="white", high="purple")
        nodes_fill_names <- names(sort(nodes_fill))
        nodes_fill <- color_vector
        names(nodes_fill) <- nodes_fill_names
        nodes_fill <- nodes_fill[match(nodes, names(nodes_fill))]        
      }
    }
    
    E(graph)$weight <- 2.5
    
    V(graph)$color <- unname(nodes_fill)[match(V(graph)$name, names(nodes_fill))]
    V(graph)$vertex.shape <- node_shapes
    
    nodes_with_edges <- unique(as.vector(get.edgelist(graph)))
    V(graph)$size <- ifelse(V(graph)$name %in% nodes_with_edges, 7, 3)
    
    
    #  labs <- radian.rescale(x=1:length(V(graph_temp)$name), direction = -1, start=0)
    
    if(layout == "circular"){
      #coords <- layout_in_circle(graph, order = order(nodes_ages, V(graph)$name))
      coords <- layout_in_circle(graph)
      #V(graph)$name[seq(2,length(V(graph)$name),2)] <- paste("", V(graph)$name[seq(2,length(V(graph)$name),2)], sep="\n")
      lab.locs <- radian.rescale(x=1:length(V(graph)$name), direction=-1, start=0)
    } else if (layout =="star"){
      coords <- layout.star(graph, center=center)
    }
    
    V(graph)$name <- sapply(strsplit(V(graph)$name, " "), function(x){
      if(length(x) != 1){
        paste(x, "", collapse="\n")        
      } else{
        return(paste(x, "", collapse=" "))
      }
    })
    if(network_type == "multiple" && layout=="circular"){
      
      par(bg = "black")
      V(graph)$label.cex = .55  #.55
      p <- plot.igraph(graph, layout = coords, vertex.label=NA, vertex.label.color="white", 
                       vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                       main=paste(tumour, network_database, sep="\n"),
                       vertex.shape= V(graph)$vertex.shape,
                       vertex.frame.color="white",
                       vertex.label.dist=.45,
                       vertex.label.degree=lab.locs)   
      #title(paste(tumour, network_database, type, sep="\n"), col.main = "white")
    } else if (network_type == "individual" && layout=="circle"){
      par(bg="white")
      V(graph)$label.cex = 1
      plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="black", 
                  vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                  vertex.shape= V(graph)$vertex.shape,
                  vertex.frame.color="white", 
                  vertex.label.dist=.45,
                  vertex.label.degree=lab.locs,
                  vertex.size=20)     
      title(tumour, col.main = "black")
    } else if (network_type == "individual" && layout=="star"){
      par(bg="white")
      V(graph)$label.cex =1
      plot.igraph(graph, layout = coords, vertex.label=V(graph)$name, vertex.label.color="black", 
                  vertex.color=V(graph)$color, edge.width=E(graph)$weight, 
                  vertex.shape= V(graph)$vertex.shape,
                  vertex.frame.color="white", 
                  vertex.label.dist=.45,
                  vertex.size=20)     
      title(tumour, col.main = "black")
      
    }
  } else{
    par(bg = "black")
    g <- barabasi.game(10)
    plot(g, vertex.label.color="black", vertex.color="black", edge.color="black")
  }
}