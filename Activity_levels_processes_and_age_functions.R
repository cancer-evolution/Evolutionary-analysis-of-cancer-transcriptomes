prepare_QSarray_object <- function(expression_data, labels, contrast, pathways)
{
  #Step-by-step
  comparison <- makeComparison(expression_data, labels, contrast, var.equal=FALSE)
  gene_set_aggregation <- aggregateGeneSet(comparison, pathways)
  results_with_VIF <- calcVIF(expression_data, gene_set_aggregation, useCAMERA=FALSE)
  return(results_with_VIF)
}

calculate_p_values <- function(results_with_VIF, test_type)
{
  
  log.fold.change = results_with_VIF$path.mean
  if ("selfcontained" %in% test_type){
    p.Value = pdf.pVal(results_with_VIF, alternative=c("two.sided"),
                       direction=FALSE, addVIF = TRUE, selfContained=TRUE)   
  } else if(test_type == "competitive"){
    p.Value = pdf.pVal(results_with_VIF, alternative=c("two.sided"),
                       direction=FALSE, addVIF = TRUE, selfContained=FALSE)
  }
  FDR = p.adjust(p.Value, method="fdr")
  pathway.name = names(log.fold.change)
  
  results = data.frame(pathway.name,log.fold.change,p.Value,FDR)
  rownames(results) = NULL
  return(results)
}

qusage_script <- function(tumour_type, pathway_type)
{
  #Prepare expression data following limma pipeline (removing lowly expressed genes)
  #and use voom object for qusage
  
  expression_normal <- #concatenated gene expression from TCGA across all patients.
  expression_tumour <- #concatenated gene expression from TCGA across all patients.
  
  samples_low_correlation <- #Samples to be removed (those with low correlation with other samples) 
  
  if (substr(samples_low_correlation[1],9,12) != "")  #if there are samples to be removed
  {
    expression_normal <- expression_normal[,!(substr(colnames(expression_normal),9,12) %in% substr(samples_low_correlation,9,12))]
    expression_tumour <- expression_tumour[,!(substr(colnames(expression_tumour),9,12) %in% substr(samples_low_correlation,9,12))]
  }
  
  expression_data <- cbind(expression_normal,expression_tumour[,2:ncol(expression_tumour)])
  rownames(expression_data) <- sapply(as.character(expression_data[,1]), function(x) {
    unname(unlist(strsplit(x,"|", fixed=TRUE)))[2] })
  
  expression_type <- c(rep(0,ncol(expression_normal)-1),rep(1,ncol(expression_tumour)-1))
  expression_header <- c("geneID",rep(0,ncol(expression_normal)-1),rep(1,ncol(expression_tumour)-1))
  colnames(expression_data) <- expression_header
  
  expression_y <- DGEList(counts=expression_data[,2:ncol(expression_data)], genes=expression_data[,1])
  high_genes_expression_normal <- which(rowSums(cpm(expression_y[,1:(ncol(expression_normal)-1)])>1)==(ncol(expression_normal)-1))
  high_genes_expression_tumour <- which(rowSums(cpm(expression_y[,1:(ncol(expression_tumour)-1)])>1)==(ncol(expression_tumour)-1))
  
  #Genes that have a cpm > 1 in all expression_normal and expression_tumour samples
  expression_all_high <- which((rowSums(cpm(expression_y)>1) == ncol(expression_y)))
  
  #Genes to be kept
  expression_genes_to_keep <- c(high_genes_expression_normal, high_genes_expression_tumour, expression_all_high)
  
  expression_y <- expression_y[sort(unique(expression_genes_to_keep)),]
  length(sort(unique(expression_genes_to_keep)))
  
  expression_y$samples$lib.size <- colSums(expression_y$counts)
  
  expression_y <- calcNormFactors(expression_y, method="upperquartile")
  
  expression_design <- model.matrix(~expression_type)
  
  expression_v <- voom(expression_y,expression_design,plot=FALSE)
  
  
  ######qusage
  
  print (paste("Starting qusage of ", tumour_type, "...", sep=""))
  
  expression_data_qusage <- expression_v$E
  expression_data <- as.matrix(expression_data_qusage)
  
  labels <- as.factor(c(rep("N",ncol(expression_normal)-1),rep("T",ncol(expression_tumour)-1)))
  contrast <- "T-N"
  

  slim_terms <- getOBOCollection("goslim_generic.obo") # Downloaded from GO website
  go_bp_children <- GOBPOFFSPRING$"GO:0008150"
    
  slim_terms <- ids(slim_terms)[which(ids(slim_terms) %in% go_bp_children)]
    
  gonames <- Term(GOTERM)
  slim_terms <- cbind(slim_terms, unname(gonames[match(slim_terms, names(gonames))])) 
  load(pathways) #list with GOID as names, entrez as items
  pathways <- pathways[slim_terms[,1]]
  pathways <- pathways[!is.na(names(pathways))]
  names(pathways) <- unname(gonames[match(names(pathways), names(gonames))])

  #save(pathways, file="GOslims_parents.Rdata") #entrez id

  QSarray_object <- prepare_QSarray_object(expression_data, labels, contrast, pathways)
  print (paste("Finished preparing QSarray object of ", tumour_type, "!", sep=""))
  
  qusage_results_sc <- calculate_p_values(QSarray_object, "selfcontained")
  filename_sc = paste(tumour_type, "_", pathway_type, "_selfcontained_qusage_results_parents.txt", sep="")
  write.table(qusage_results_sc, file = filename_sc, sep="\t", row.names=FALSE, quote=FALSE)
  print (paste("Finished selfcontained test of ", tumour_type, "!", sep=""))
  
  qusage_results_comp <- calculate_p_values(QSarray_object, "competitive")
  filename_comp = paste(tumour_type, "_", pathway_type, "_competitive_qusage_results_parents.txt", sep="")
  write.table(qusage_results_comp, file = filename_comp, sep="\t", row.names=FALSE, quote=FALSE)
  print (paste("Finished competitive test of ", tumour_type, "!", sep=""))
  
}

ssgsea_script <- function(type, pathway_type, paired){
  expression_data <- prepare_expression(type, paired)
  if (pathway_type == "GOslims"){
    load(pathways) #list with GOID as names, entrez as items 
  }
  return(gsva(as.matrix(expression_data), pathways, method="ssgsea", verbose=FALSE, rnaseq=TRUE, ssgsea.norm=TRUE))
}

prepare_expression <- function(type, paired){
  expression_normal <- #concatenated gene expression from TCGA across all patients.
  expression_tumour <- #concatenated gene expression from TCGA across all patients.
  samples_low_correlation <- #samples to be removed
  
  rownames(expression_normal) <- expression_normal[,1]
  rownames(expression_tumour) <- expression_tumour[,1]
  
  if (substr(samples_low_correlation[1],9,12) != "")  #if there are samples to be removed
  {
    expression_normal <- expression_normal[,!(substr(colnames(expression_normal),9,12) %in% substr(samples_low_correlation,9,12))]
    expression_tumour <- expression_tumour[,!(substr(colnames(expression_tumour),9,12) %in% substr(samples_low_correlation,9,12))]
  }
  
  if(paired==TRUE){
    normal_samples <- substr(colnames(expression_normal), 9,12)
    tumour_samples <- substr(colnames(expression_tumour), 9,12)
    expression_normal <- expression_normal[,which(normal_samples %in% tumour_samples)]
    expression_tumour <- expression_tumour[,which(tumour_samples %in% normal_samples)]
    expression_normal <- expression_normal[,match(substr(colnames(expression_tumour), 9,12),
                                                 substr(colnames(expression_normal), 9,12))]
  }
  
  expression_data <- cbind(expression_normal,expression_tumour[,2:ncol(expression_tumour)])
  expression_data_header <- c("gene_id", paste(substr(colnames(expression_normal)[2:ncol(expression_normal)], 9,12), "_N", sep=""),
                              paste(substr(colnames(expression_tumour)[2:ncol(expression_tumour)], 9,12), "_T", sep=""))
  
  colnames(expression_data) <- expression_data_header
  expression_data[,1] <- sapply(as.character(expression_data[,1]), function(x) {
    unname(unlist(strsplit(x,"|", fixed=TRUE)))[2] })
  expression_data <- expression_data[!(expression_data[,1] == "?"),]
  expression_data <- expression_data[!duplicated(expression_data[,1]),]
  rownames(expression_data) <- expression_data[,1]
  expression_data <- expression_data[,-1]
  
  return(expression_data)
}

load_TAI <- function(){
  TAI_l <- list()
  tumours <- c("LUAD", "LUSC", "BRCA","PRAD","LIHC", "COAD", "STAD") 
  for (i in tumours)
  {
    TAI_normal <- #TAI of normal samples. Colums: Patient ID, raw_TAI, corrected_TAI, Phylostrata proportions 1:16
    TAI_normal <- cbind(paste(as.character(TAI_normal[,1]), "_N", sep=""), TAI_normal[,3])
    
    TAI_tumour <- #TAI of tumour samples. Colums: Patient ID, raw_TAI, corrected_TAI, Phylostrata proportions 1:16
    TAI_tumour <- cbind(paste(as.character(TAI_tumour[,1]), "_T", sep=""), TAI_tumour[,3])
    
    TAI <- rbind(TAI_normal, TAI_tumour)
    colnames(TAI) <- c("Patient", "TAI")
    TAI <- TAI[order(TAI[,2], decreasing=TRUE),]
    TAI_l[[i]] <- TAI
  }
  return(TAI_l)
}

load_clinical_data <- function(tumour){
  clinical_data_TCGA <- read.delim(paste("nationwidechildrens.org_clinical_patient_", tumour, ".txt", sep=""))  #Downloaded from TCGA's website
  
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
  clinical_data <- droplevels(clinical_data)
  return(clinical_data)
}

qusage_script_explicit <- function(pathways, expression_data, expression_design){
  
  expression_design <- unname(expression_design[,2])
  labels <- as.factor(c(rep("N",sum(expression_design == 0)),rep("T",sum(expression_design == 1))))
  contrast <- "T-N"
  
  for( x in (names(pathways))){
    if(length(pathways[[x]]) == 0){
      pathways[[x]] <- NULL
    } else if (sum(pathways[[x]] %in% rownames(expression_data)) <= 1){
      pathways[[x]] <- NULL
    }
  }
  if(class(pathways) == "matrix"){
    t <- colnames(pathways)
    colnames(pathways) <- NULL
    pathways <- list(as.vector(pathways))
    names(pathways) <- t
  }
  
  if(length(pathways) != 0){
    QSarray_object <- prepare_QSarray_object(expression_data, labels, contrast, pathways)
    
    qusage_results_sc <- calculate_p_values(QSarray_object, "selfcontained")
    return(qusage_results_sc)
  }
}

calculate_difference_between_UC_MC_components <- function(start_matrix, absolute){
  temp_matrix <- start_matrix[duplicated(interaction(start_matrix$GOslims, start_matrix$Tumour))|
                                duplicated(interaction(start_matrix$GOslims, start_matrix$Tumour), fromLast=TRUE) ,]
  if (absolute==TRUE){
    temp_matrix$FC <- abs(temp_matrix$FC)
  }
  
  return_matrix <- data.frame(temp_matrix[c(TRUE,FALSE),c(1,3)],
                              difference_log_FC = diff(temp_matrix$FC)[c(TRUE,FALSE)],
                              row.names=NULL)
  
  return_matrix <- data.frame(return_matrix,
                              Number_genes = unname(apply(return_matrix,1,function(row){
                                slim <- as.character(row[1])
                                return(c(length(slim_genes[[slim]])))})),
                              SlimAge= GO_ages[match(return_matrix[,1],
                                                         GO_ages[,1]),9])
  return_matrix$SlimAge <- factor(return_matrix$SlimAge, levels=c("UC", "MC"))
  return_matrix$difference_log_FC <- return_matrix$difference_log_FC*(-1)
  
  qusage_slims <- data.frame()
  for(tumour_type in tumours){
    temp <- read.delim(paste(tumour_type,"_goslims_selfcontained_qusage_results.txt", sep=""))
    column_names <- colnames(temp)
    temp <- t(apply(temp, 1, as.character))
    qusage_slims <- rbind(qusage_slims, data.frame(temp, Tumour=tumour_type))
  }
  colnames(qusage_slims)[1:4] <- column_names
  
  return_matrix <- data.frame(return_matrix,
                              qusage_slims[match(
                                interaction( return_matrix$GOslims, return_matrix$Tumour),
                                interaction( qusage_slims$pathway.name, qusage_slims$Tumour)
                              ),2])
  colnames(return_matrix)[6] <- "log_FC_of_slim"
  return_matrix[,6] <- as.numeric(as.character(return_matrix[,6]))
  
  return_matrix$SlimAge <- factor(return_matrix$SlimAge, levels=c(levels(return_matrix$SlimAge), "NA"))
  return_matrix$SlimAge <- replace(return_matrix$SlimAge,
                                   is.na(return_matrix$SlimAge),
                                   "NA")
  
  return(return_matrix)
}

create_df_of_error_bars <- function(df){
  error_bars_difference_df <- aggregate(difference_log_FC~GOslims, df, range)
  error_bars_FC_df <- aggregate(log_FC_of_slim~GOslims, df, range)
  FC_median_df <- aggregate(log_FC_of_slim~GOslims, df, median)
  difference_median_df <- aggregate(difference_log_FC~GOslims, df, median)
  
  error_bars_df <- data.frame(FC_median_df[,1], FC_median_df[,2], error_bars_FC_df[,2],
                              difference_median_df[,2], error_bars_difference_df[,2])
  colnames(error_bars_df) <- c("GOSlim", "Median_FC","Min_logFC","Max_logFC",
                               "Median_logFC_difference","Min_logFC_difference", "Max_logFC_difference")
  error_bars_df <- data.frame(error_bars_df, SlimAge=GO_ages[match(error_bars_df$GOSlim, 
                                                                       GO_ages$GOSlims), "Age"])
  error_bars_df$SlimAge <- factor(error_bars_df$SlimAge, levels=c("UC", "MC"))
  return(error_bars_df)
}

plot_error_FC_diff_error_bars <- function(df){
  ggplot(df, aes(x=Median_FC, y=Median_logFC_difference))+
    geom_errorbar(aes(ymin=Min_logFC_difference, ymax=Max_logFC_difference, colour=Direction), width=0, size=.25)+
    geom_errorbarh(aes(xmin=Min_logFC, xmax=Max_logFC, colour=Direction), height =0, size=.25)+
    scale_shape_manual(values=c(21,24), labels=c("UC slim", "MC slim"))+
    geom_point(aes(fill=Direction, colour=Direction, shape=Driven),size=4)+
    geom_point(aes(fill=Direction, shape=Driven), colour="black", size=4, show.legend=FALSE)+
    ##geom_point(aes(colour=Direction, fill=Driven, shape=SlimAge, size=Degree_ratio))+
    ##geom_point(aes(fill=Direction, colour=Driven, shape=SlimAge, size=log_degree_ratio), show.legend=FALSE)+
    #scale_colour_manual(values = c("red","blue"))+
    #scale_fill_manual(name = "Driven by",
    #                  values = c("green","purple"),
    #                  labels=c("Driven by UC", "Driven by MC"))+
    #scale_size_continuous(range = c(0.1,10))+
    geom_vline(xintercept = 0, colour="black", linetype = "longdash", size=0.5)+
    geom_hline(yintercept = 0, colour="black", linetype = "longdash", size=0.5)+
    #geom_text_repel(data=subset(df, (Median_logFC_difference>0 & Median_FC < 0) |
    #                              (Median_logFC_difference<0 & Median_FC > 0)), aes(label=GOSlim), size=3.5)+
    #geom_text_repel(data=subset(df, Degree_ratio < 1), aes(label=GOSlim), size=3.5)+
    #geom_text_repel(data=df, aes(label=GOSlim), size=3.5)+
    scale_y_continuous(breaks=seq(-1.25, 0.5, 0.25))+
    scale_x_continuous(breaks=seq(-2, 0.75, 0.25))+
    xlab("Median logFC of cellular process")+
    ylab("Median difference of the logFC of\nthe UC and MC components")+
    theme_bw()+
    guides(colour=NULL)+
    theme(legend.title=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size=18),axis.text.x= element_text(size=19, angle = 45, hjust = 1),
          axis.text.y= element_text(size=19),
          axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
          legend.text=element_text(size=17), legend.title=element_text(size=17))
}

print_trees_attributes_pathway_of_interest <- function(pathway_of_interest){
  children_of_pathway_of_interest <- GOBOFFSPRING_is_a[[pathway_of_interest]]
  id_interest <- names(gonames)[match(pathway_of_interest, unname(gonames))]
  children_of_pathway_of_interest <- c(children_of_pathway_of_interest, id_interest)
  number_genes <- go_genes[children_of_pathway_of_interest]
  number_genes <- number_genes[!is.na(names(number_genes))]
  number_genes <- sapply(number_genes, length)
  
  pathway_of_interest_df <- cbind(unname(gonames)[match(children_of_pathway_of_interest, names(gonames))],
                                  children_of_pathway_of_interest, 
                                  number_genes[match(children_of_pathway_of_interest, names(number_genes))])
  rownames(pathway_of_interest_df) <- NULL
  
  
  for(type in tumours){
    pathway_of_interest_df <- cbind(pathway_of_interest_df,
                                    qusage_slim_children[[type]][[pathway_of_interest]][match(pathway_of_interest_df[,1], qusage_slim_children[[type]][[pathway_of_interest]][,1]),2])
  }
  
  colnames(pathway_of_interest_df) <- c("GOterm", "GOID", "Number of genes", tumours)
  pathway_of_interest_df <- as.data.frame(pathway_of_interest_df)
  pathway_of_interest_df[4:ncol(pathway_of_interest_df)] <- apply(pathway_of_interest_df[4:ncol(pathway_of_interest_df)],2,as.numeric) 
  pathway_of_interest_df <- cbind(pathway_of_interest_df, 
                                  rowMeans(pathway_of_interest_df[,4:ncol(pathway_of_interest_df)], na.rm=TRUE))
  
  pathway_of_interest_df <- pathway_of_interest_df[!is.na(pathway_of_interest_df[ncol(pathway_of_interest_df)]),]
  colnames(pathway_of_interest_df)[ncol(pathway_of_interest_df)] <- "Mean_logFC"
  
  write.table(pathway_of_interest_df, file=paste(pathway_of_interest, "_qusageFC2.txt", sep=""), quote=FALSE,
              sep="\t", row.names=FALSE)
  
  #Trim tree to only include nodes with logFC and their parents if necessary
  
  all_parents <- unique(unname(unlist(GOBPANCESTOR_is_a[pathway_of_interest_df[,1]])))
  all_parents <- unname(gonames)[match(all_parents, names(gonames))]
  all_parents <- unique(c(all_parents, as.character(pathway_of_interest_df[,1])))
  
  all_children_pathway_interest <- unname(gonames)[match(unlist(GOBOFFSPRING_is_a[pathway_of_interest]), names(gonames))]
  all_parents <- all_parents[all_parents %in% all_children_pathway_interest]
  all_parents <- c(all_parents, pathway_of_interest)
  go_interactions_interest <- go_interactions_is_a[intersect(which(go_interactions_is_a[,1] %in% all_parents), which(go_interactions_is_a[,3] %in% all_parents)),]
  
  write.table(go_interactions_interest, file = paste(pathway_of_interest, "_trimmed_tree.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
}

extract_network_components <- function(slim, network, slim_genes){
  local_genes <- slim_genes[[slim]]
  local_network <- network[intersect(which(network[,"Gene1"] %in% local_genes),
                                     which(network[,"Gene2"] %in% local_genes)),]
  return(local_network)
}