#Functions for figure 1

read_expression <- function(tumour_type, sample_type){
   exp <-   #concatenated gene expression from TCGA across all patients. First column geneID, second column, gene phylostrata. After that, expression values
  if (sample_type == "normal"){
    colnames(exp)[3:ncol(exp)] <- paste(substr(colnames(exp)[3:ncol(exp)],9,12), "_N", sep="")
  } else{
    colnames(exp)[3:ncol(exp)] <- paste(substr(colnames(exp)[3:ncol(exp)],9,12), "_T", sep="")
  }
  samples_low_correlation <- readLines(paste(tumour_type, "_samples_low_correlation.txt", sep=""))  #Samples to be removed (those with low correlation with other samples) 
  if (substr(samples_low_correlation[1],9,12) != "")  #if there are samples to be removed
  {
    exp <- exp[,!(substr(colnames(exp),1,4) %in% substr(samples_low_correlation,9,12))]
  }
  return(exp)
}

read_expression_for_Domazet <- function(tumour_type, sample_type){
  exp <- read.delim(paste(tumour_type, "_", sample_type, "_normalized_expression.txt", sep=""))
  if (sample_type == "normal"){
    colnames(exp)[2:ncol(exp)] <- paste(substr(colnames(exp)[2:ncol(exp)],9,12), "_N", sep="")
  } else{
    colnames(exp)[2:ncol(exp)] <- paste(substr(colnames(exp)[2:ncol(exp)],9,12), "_T", sep="")
  }
  samples_low_correlation <- readLines(paste(tumour_type, "_samples_low_correlation.txt", sep=""))
  if (substr(samples_low_correlation[1],9,12) != "")  #if there are samples to be removed
  {
    exp <- exp[,!(substr(colnames(exp),1,4) %in% substr(samples_low_correlation,9,12))]
  }
  gene_ids <- unlist(strsplit(as.character(exp$gene_id), "\\|"))
  entrez <- gene_ids[seq(2,length(gene_ids), 2)]
  exp$gene_id <- entrez
  
  gene_ages <- read.delim("Domazet_gene_ages.txt") #Gene ages downloaded from Domazet paper
  gene_ages_entrez <- gene_annotations[match(gene_ages$Gene_ID, gene_annotations$Ensembl.Gene.ID), "EntrezGene.ID"]
  exp_ages <- cbind(gene_ages_entrez, gene_ages[,"Phylostratum"])
  
  exp <- cbind(exp[,1], exp_ages[match(as.character(exp[,1]), as.character(exp_ages[,1])), 2], exp[,2:ncol(exp)])
  colnames(exp)[1:2] <- c("GeneID", "Phylostratum")
  exp <- exp[!is.na(exp$Phylostratum),]
  return(exp)
}

find_genes_with_go <- function(high_level){
  if (substr(high_level,1,3) != "GO:"){
    GO_id <- names(gonames)[which(unname(gonames) == high_level)]
  } else{
    GO_id <- high_level
  }
  GO_id_off <- as.list(GOBPOFFSPRING)[[GO_id]]
  return(unique(gene_association[which(gene_association[,5] %in% GO_id_off),3]))
}

prepare_TAI_ESTIMATE <- function(TAI_all_raw, tumour_type){
  ESTIMATE_scores_normal <- #Estimate scores for normal samples in .gct
  colnames(ESTIMATE_scores_normal) <- as.character(unlist(unname(ESTIMATE_scores_normal[2,])))
  ESTIMATE_scores_normal <- ESTIMATE_scores_normal[-c(1,2,3,4),3:ncol(ESTIMATE_scores_normal)]
  colnames(ESTIMATE_scores_normal) <- paste(substr(colnames(ESTIMATE_scores_normal),9,12), "_N", sep="")
  
  ESTIMATE_scores_tumour <- #Estimate scores for normal samples in .gct
  colnames(ESTIMATE_scores_tumour) <- as.character(unlist(unname(ESTIMATE_scores_tumour[2,])))
  ESTIMATE_scores_tumour <- ESTIMATE_scores_tumour[-c(1,2,3,4),3:ncol(ESTIMATE_scores_tumour)]
  colnames(ESTIMATE_scores_tumour) <- paste(substr(colnames(ESTIMATE_scores_tumour),9,12), "_T", sep="")
  
  ESTIMATE_scores_all <- as.numeric(c(as.matrix(ESTIMATE_scores_normal[1,]), as.matrix(ESTIMATE_scores_tumour[1,])))
  names(ESTIMATE_scores_all) <-c(colnames(ESTIMATE_scores_normal), colnames(ESTIMATE_scores_tumour))
  
  TAI_ESTIMATE <- cbind(TAI_all_raw, ESTIMATE_scores_all[match(rownames(TAI_all_raw), names(ESTIMATE_scores_all))])
  colnames(TAI_ESTIMATE)[ncol(TAI_ESTIMATE)] <- c("Estimate_score")
  TAI_ESTIMATE <- cbind(TAI_ESTIMATE, rownames(TAI_ESTIMATE))
  colnames(TAI_ESTIMATE)[ncol(TAI_ESTIMATE)] <- c("Patient")
  return(TAI_ESTIMATE)
}

prepare_TAI_ESTIMATE_new_prostate <- function(TAI_all_raw, tumour_type){
  ESTIMATE_scores <- #Read in ESTIMATE results in gct
  colnames(ESTIMATE_scores) <- as.character(unlist(unname(ESTIMATE_scores[2,])))
  ESTIMATE_scores <- ESTIMATE_scores[-c(1,2,3,4),3:ncol(ESTIMATE_scores)]
  ESTIMATE_scores <- unlist(ESTIMATE_scores[1,])
  ESTIMATE_scores_names <- names(ESTIMATE_scores)
  ESTIMATE_scores <- as.numeric(as.character(ESTIMATE_scores))
  names(ESTIMATE_scores) <- ESTIMATE_scores_names
  
  TAI_ESTIMATE <- cbind(TAI_all_raw, ESTIMATE_scores[match(rownames(TAI_all_raw), names(ESTIMATE_scores))])
  colnames(TAI_ESTIMATE)[ncol(TAI_ESTIMATE)] <- c("Estimate_score")
  TAI_ESTIMATE <- cbind(TAI_ESTIMATE, rownames(TAI_ESTIMATE))
  colnames(TAI_ESTIMATE)[ncol(TAI_ESTIMATE)] <- c("Patient")
  return(TAI_ESTIMATE)
}
step_TAI_calculation <- function(expression){
  step1 <- apply(expression[,3:ncol(expression)],2,function(sample){
    return(sample*expression[,2])
  })
  step2 <- apply(step1,2,sum)
  step3 <- step2/apply(expression[,3:ncol(expression)],2,sum)
  names(step3) <- colnames(expression)[3:length(colnames(expression))]
  
  #phylostrata proportions
  result <- step3
  total_library_size <- apply(expression[,3:ncol(expression)], 2, sum)
  for(i in 1:16){
    temp <- expression[which(expression[,2] == i),]
    temp2 <- apply(temp[,3:ncol(temp)], 2, sum)
    temp3 <- temp2/total_library_size
    result <- cbind(result, temp3)
  }
  colnames(result) <- c("TAI", paste("Phy", 1:16, sep="_"))
  return(result)
}

calculate_TAI <- function(expression, tumour_type){
  TAI <- step_TAI_calculation(expression)
  TAI_ESTIMATE <- prepare_TAI_ESTIMATE(TAI, tumour_type)
  
  #rownames(TAI_normal_raw) <- paste(substr(rownames(TAI_normal_raw),9,12), "_N", sep="")
  #rownames(TAI_tumour_raw) <- paste(substr(rownames(TAI_tumour_raw),9,12), "_T", sep="")
  
  return(TAI_ESTIMATE)
}

calculate_TAI_new_prostate <- function(expression, tumour_type){
  TAI <- step_TAI_calculation(expression)
  TAI_ESTIMATE <- prepare_TAI_ESTIMATE_new_prostate(TAI, tumour_type)
  
  #rownames(TAI_normal_raw) <- paste(substr(rownames(TAI_normal_raw),9,12), "_N", sep="")
  #rownames(TAI_tumour_raw) <- paste(substr(rownames(TAI_tumour_raw),9,12), "_T", sep="")
  
  return(TAI_ESTIMATE)
}

plotPhyDist <- function (tumour_type, difference_phy_means, sum_phy_variances, ylim){
  plot(0, 0, type="n", ann=FALSE, axes=FALSE, ylim=ylim)
  rect(-1.5,-0.075,17,0,col="aliceblue",border='black')
  rect(-1.5,0,17,0.075,col="lavenderblush",border='black')
  par(new = TRUE)
  plot(difference_phy_means, type="l", ylim=ylim, lwd=2, col="black", xaxt="n", xlab="Phylostrata", ylab="Difference of expression proportions", main=tumour_type)
  axis(1, at=1:16, labels=1:16)
  arrows(1:16, difference_phy_means-sum_phy_variances, 1:16, difference_phy_means+sum_phy_variances, code=3, angle=90, length=0.1, col='black')
}

loess_correction_plus_median <- function(ESTIMATE_TAI_type){
  ESTIMATE_TAI_all <- rbind(ESTIMATE_TAI_type$normal, ESTIMATE_TAI_type$tumour)
  ESTIMATE_TAI_all <- as.data.frame(ESTIMATE_TAI_all)
  ESTIMATE_TAI_all$TAI <- as.numeric(as.character(ESTIMATE_TAI_all$TAI))
  ESTIMATE_TAI_all$Estimate_score <- as.numeric(as.character(ESTIMATE_TAI_all$Estimate_score))
  loess_model <- loess(ESTIMATE_TAI_all$TAI~ESTIMATE_TAI_all$Estimate_score)
  TAI_predictions <- predict(loess_model, data.frame(x=ESTIMATE_TAI_all$Estimate_score))
  TAI_corrected_all <- ESTIMATE_TAI_all$TAI-TAI_predictions
  
  ESTIMATE_TAI_type$normal <- as.data.frame(ESTIMATE_TAI_type$normal)
  ESTIMATE_TAI_type$tumour <- as.data.frame(ESTIMATE_TAI_type$tumour)
  ESTIMATE_TAI_type$normal$TAI <- as.numeric(as.character(ESTIMATE_TAI_type$normal$TAI))
  ESTIMATE_TAI_type$tumour[,"TAI"] <- as.numeric(as.character(ESTIMATE_TAI_type$tumour$TAI))
  
  median_all_raw_TAI <- median(c(ESTIMATE_TAI_type$normal[,"TAI"], ESTIMATE_TAI_type$tumour[,"TAI"]))
  TAI_corrected_all_median <- TAI_corrected_all+median_all_raw_TAI
  
  TAI_corrected_normal <- as.numeric(TAI_corrected_all_median[1:nrow(ESTIMATE_TAI_type$normal)])
  TAI_corrected_tumour <- as.numeric(TAI_corrected_all_median[(nrow(ESTIMATE_TAI_type$normal)+1):length(TAI_corrected_all_median)])
  
  tumour_TAI <- as.numeric(as.character(ESTIMATE_TAI_type$tumour[,"TAI"]))
  tumour_ESTIMATE <- as.numeric(as.character(ESTIMATE_TAI_type$tumour[,"Estimate_score"]))
  tumour_phy <- ESTIMATE_TAI_type$tumour[,2:17]
  
  normal_TAI <- as.numeric(as.character(ESTIMATE_TAI_type$normal[,"TAI"]))
  normal_ESTIMATE <- as.numeric(as.character(ESTIMATE_TAI_type$normal[,"Estimate_score"]))
  normal_phy <- ESTIMATE_TAI_type$normal[,2:17]
  
  all_TAI <- c(normal_TAI, tumour_TAI)
  all_ESTIMATE <- c(normal_ESTIMATE, tumour_ESTIMATE)
  all_phy <- rbind(normal_phy, tumour_phy)
  all_phy <- apply(all_phy, 2, as.numeric, as.character)
  
  loess_model_phy <- rep(NA,16)
  phy_predictions <- matrix(nrow=nrow(all_phy), ncol=16)
  phy_corrected <- matrix(nrow=nrow(all_phy), ncol=16)
  for (i in 1:16)
  {
    loess_model_phy <- loess(all_phy[,i]~all_ESTIMATE)
    phy_predictions[,i] <- predict(loess_model_phy, data.frame(x=all_ESTIMATE))
    phy_corrected[,i] <- all_phy[,i]-phy_predictions[,i]
    phy_corrected[,i] <- phy_corrected[,i] + median(all_phy[,i])
  }
  
  phy_corrected_normal <- phy_corrected[1:nrow(normal_phy),]
  phy_corrected_tumour <- phy_corrected[(nrow(normal_phy)+1):nrow(phy_corrected),]
  
  differences_phy_corrected <- rep(NA,16)
  
  for (i in 1:16)
  {
    differences_phy_corrected[i] <- mean(phy_corrected_tumour[,i])-mean(phy_corrected_normal[,i])
  }
  variance_normal = apply(phy_corrected_normal, 2, var)
  variance_tumour = apply(phy_corrected_tumour, 2, var)
  sum_phy_variances = variance_normal + variance_tumour
  
  df_tumour <- cbind(as.character(ESTIMATE_TAI_type$tumour$Patient), TAI_corrected_tumour, phy_corrected_tumour)
  df_normal <- cbind(as.character(ESTIMATE_TAI_type$normal$Patient), TAI_corrected_normal, phy_corrected_normal)
  colnames(df_tumour) <- c("Patient", "TAI", paste("Phy", 1:16, sep="_"))
  colnames(df_normal) <- c("Patient", "TAI", paste("Phy", 1:16, sep="_"))
  df_tumour <- as.data.frame(df_tumour)
  df_normal <- as.data.frame(df_normal)
  df_tumour[,3:18] <- apply(df_tumour[,3:18], 2, function(x){ as.numeric(as.character(x))})
  df_normal[,3:18] <- apply(df_normal[,3:18], 2, function(x){ as.numeric(as.character(x))})
  return (list(tumour = df_tumour, normal = df_normal))
}

loess_correction_plus_median_new_prostate <- function(ESTIMATE_TAI_all){
  ESTIMATE_TAI_all <- as.data.frame(ESTIMATE_TAI_all)
  ESTIMATE_TAI_all$TAI <- as.numeric(as.character(ESTIMATE_TAI_all$TAI))
  ESTIMATE_TAI_all$Estimate_score <- as.numeric(as.character(ESTIMATE_TAI_all$Estimate_score))
  loess_model <- loess(ESTIMATE_TAI_all$TAI~ESTIMATE_TAI_all$Estimate_score)
  TAI_predictions <- predict(loess_model, data.frame(x=ESTIMATE_TAI_all$Estimate_score))
  TAI_corrected_all <- ESTIMATE_TAI_all$TAI-TAI_predictions
  
  median_all_raw_TAI <- median(ESTIMATE_TAI_all$TAI)
  TAI_corrected_all_median <- TAI_corrected_all+median_all_raw_TAI
  
  df <- cbind(as.character(ESTIMATE_TAI_all$Patient), TAI_corrected_all_median)
  colnames(df) <- c("Patient", "TAI")
  df <- as.data.frame(df)
  return (df)
}

remove_signature_TAI <- function(signature, expression, type){
  print(dim(expression))
  expression_no_signature <- expression[!(expression[,1] %in% signature),]
  print(dim(expression_no_signature))
  TAI_no_signature <- calculate_TAI(expression_no_signature, type)
  
  return(TAI_no_signature)
}

plot_TAI <- function(corrected_TAI_df){
  g_TAI <- ggplot(corrected_TAI_df, aes(x=Tissue_type, y=TAI))+
    theme_bw()+
    geom_boxplot(aes(fill=Tissue_type))+
    #geom_point(position = position_jitter(width = 0.4), size=.2) +
    facet_grid(.~Tumour_type)+
    guides(fill=FALSE)+
    xlab("")+
    ylim(3, 6.4)+
    #ylim(2.5, 5.4)+  #DOMAZET
    theme(strip.background = element_rect(colour="black", fill="white"), 
          strip.text.x = element_text(size=15, face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20),
          axis.text.y= element_text(size=20),
          axis.title.y = element_text(size=20),
          text = element_text(size=18))+
    geom_line(data = significance, aes(x = a, y = b)) +
    annotate("text", x = 1.5, y = 6.25, label = "***", size = 8)
    #annotate("text", x = 1.5, y = 5.25, label = c("***", "***", "***",
    #                                              "", "***", "***", "***"), size = 8) #DOMAZET
    
  return(g_TAI)
}

distribution_TAI_bins <- function(expression_normal, expression_tumour, number_of_bins, type_bins, type){
  #number_of_genes_per_bin <- round(nrow(expression)/number_of_bins)  #100 bins
  TAI_phy_binned <- list()
  corrected_TAI_phy_binned <- list()
  if (type_bins =="specific_ranges"){
    binned_expression_normal <- bin_expression(expression_normal, number_of_bins, type, "specific_ranges")
    binned_expression_tumour <- bin_expression(expression_tumour, number_of_bins, type, "specific_ranges")
  }else {
    binned_expression_normal <- bin_expression(expression_normal, number_of_bins, type, "cummulative")
    binned_expression_tumour <- bin_expression(expression_tumour, number_of_bins, type, "cummulative")
  }
    
  corrected_TAI_phy_binned <- list()
  corrected_TAI_binned_normal <- vector()
  corrected_TAI_binned_tumour <- vector()
  
  for (i in 1:number_of_bins){
    TAI_phy_binned[[as.character(i)]][["normal"]] <- calculate_TAI(binned_expression_normal[[i]], type)
    TAI_phy_binned[[as.character(i)]][["tumour"]] <- calculate_TAI(binned_expression_tumour[[i]], type)      
    corrected_TAI_phy_binned[[as.character(i)]] <- loess_correction_plus_median(TAI_phy_binned[[as.character(i)]])
    corrected_TAI_binned_normal <- rbind(corrected_TAI_binned_normal, as.character(corrected_TAI_phy_binned[[as.character(i)]][["normal"]][,2]))
    corrected_TAI_binned_tumour <- rbind(corrected_TAI_binned_tumour, as.character(corrected_TAI_phy_binned[[as.character(i)]][["tumour"]][,2]))
  }
  corrected_TAI_binned_normal <- as.data.frame(corrected_TAI_binned_normal)
  corrected_TAI_binned_normal <- apply(corrected_TAI_binned_normal, 2, as.numeric)  
  
  corrected_TAI_binned_tumour <- as.data.frame(corrected_TAI_binned_tumour)
  corrected_TAI_binned_tumour <- apply(corrected_TAI_binned_tumour, 2, as.numeric)  
    
    
    
  corrected_TAI_binned_normal_median <- apply(corrected_TAI_binned_normal,1,median)
  corrected_TAI_binned_normal_sd <- apply(corrected_TAI_binned_normal,1,sd)
  
  corrected_TAI_binned_tumour_median <- apply(corrected_TAI_binned_tumour,1,median)
  corrected_TAI_binned_tumour_sd <- apply(corrected_TAI_binned_tumour,1,sd)
  
  return(list(TAI_normal_median=corrected_TAI_binned_normal_median, 
              TAI_normal_sd=corrected_TAI_binned_normal_sd, 
              TAI_tumour_median=corrected_TAI_binned_tumour_median, 
              TAI_tumour_sd=corrected_TAI_binned_tumour_sd))
}


bin_expression <- function(expression, number_of_bins, bin_type, order_genes){
  expression_ordered <- expression[match(order_genes, expression[,1]),]
  number_of_genes_per_bin <- nrow(expression_ordered)/number_of_bins
  bins <- list()
  if(bin_type == "specific_ranges"){
    for (i in 1:number_of_bins){
      bins[[i]] <- expression_ordered[seq((i*number_of_genes_per_bin - number_of_genes_per_bin+1),(i*number_of_genes_per_bin),1),]
    }
  } else if(bin_type == "cummulative"){
    for (i in 1:number_of_bins){
      bins[[i]] <- expression_ordered[seq(1,(i*number_of_genes_per_bin),1),]
    }
  }
  return(bins)
}


calculate_percentage_UC_bins <- function(expression_normal, expression_tumour, number_of_bins, type_bins, order_genes){
  if (type_bins =="specific_ranges"){
    binned_expression_normal <- bin_expression(expression_normal, number_of_bins, "specific_ranges", order_genes)
    binned_expression_tumour <- bin_expression(expression_tumour, number_of_bins, "specific_ranges", order_genes)
  }else {
    binned_expression_normal <- bin_expression(expression_normal, number_of_bins, type, "cummulative")
    binned_expression_tumour <- bin_expression(expression_tumour, number_of_bins, type, "cummulative")
  }
  
  binned_normal_UC <- lapply(binned_expression_normal, calculate_percentage_UC)
  binned_tumour_UC <- lapply(binned_expression_tumour, calculate_percentage_UC)
  median_binned_normal_UC <- unlist(lapply(binned_normal_UC, median))
  median_binned_tumour_UC <- unlist(lapply(binned_tumour_UC, median))
  
  minimum_binned_normal_UC <- unlist(lapply(binned_normal_UC, min))
  minimum_binned_tumour_UC <- unlist(lapply(binned_tumour_UC, min))
  
  maximum_binned_normal_UC <- unlist(lapply(binned_normal_UC, max))
  maximum_binned_tumour_UC <- unlist(lapply(binned_tumour_UC, max))
  
  sd_binned_normal_UC <- unlist(lapply(binned_normal_UC, sd))
  sd_binned_tumour_UC <- unlist(lapply(binned_tumour_UC, sd))
  
  p_values_ts <- vector()
  p_values_greater <- vector()
  p_values_less <- vector()
  for(i in 1:number_of_bins){
    temp_normal <- binned_normal_UC[[i]]
    temp_tumour <- binned_tumour_UC[[i]]
    p_values_ts <- c(p_values_ts, wilcox.test(temp_tumour, temp_normal)$p.value)
    p_values_greater <- c(p_values_greater, wilcox.test(temp_tumour, temp_normal, alternative="greater")$p.value)
    p_values_less <- c(p_values_less, wilcox.test(temp_tumour, temp_normal, alternative="less")$p.value)
    
  }
  
  return(list(normal_median=median_binned_normal_UC, tumour_median=median_binned_tumour_UC,
              normal_max=maximum_binned_normal_UC, tumour_max=maximum_binned_tumour_UC,
              normal_min=minimum_binned_normal_UC, tumour_min=minimum_binned_tumour_UC,
              normal_sd=sd_binned_normal_UC, tumour_sd=sd_binned_tumour_UC,
              p_values_ts = p_values_ts, p_values_greater = p_values_greater,
              p_values_less = p_values_less))
}

calculate_percentage_UC <- function(expression_bin){
  total_expression <- apply(expression_bin[,3:ncol(expression_bin)], 2, sum)
  local_UC_expression <- expression_bin[which(expression_bin$Phylostratum %in% c(1,2,3)),]
  UC_expression <- apply(local_UC_expression[,3:ncol(local_UC_expression)], 2, sum)
  percentage_UC <- UC_expression/total_expression*100
  return(percentage_UC)
}

plot_TAI_distributions <- function(TAI_bins_list, plot_type){
  par(mfrow=c(2,2))
  for(type in tumours){
    print(type)
    median_normal <- TAI_bins_list[[type]]$TAI_normal_median
    median_tumour <- TAI_bins_list[[type]]$TAI_tumour_median
    min_y <- min(c(median_normal, median_tumour), na.rm=TRUE)
    max_y <- max(c(median_normal, median_tumour), na.rm=TRUE)
    
    sd_normal <- TAI_bins_list[[type]]$TAI_normal_sd
    sd_tumour <- TAI_bins_list[[type]]$TAI_tumour_sd
    
    plot(median_normal, type="l", ylim=c(min_y-0.5, max_y+0.5), col="blue", ylab="TAI", xlab="Expression bins", xaxt = "n", main=type)
    lines(median_tumour, ylim=c(min_y, max_y), col="red")
    arrows(1:10, median_normal-(sd_normal^2), 1:10, median_normal+(sd_normal^2), code=3, angle=90, length=0.1, col='blue')
    arrows(1:10, median_tumour-(sd_tumour^2), 1:10, median_tumour+(sd_tumour^2), code=3, angle=90, length=0.1, col='red')
    mtext(text=round(median_normal,2), side=3, at = 1:10, line=0.75, col="blue", cex=0.75)
    mtext(text=round(median_tumour,2), side=3, at = 1:10, line=0, col="red", cex=0.75)
    if (plot_type == "specific_ranges"){
      axis(1, at = 1:10, labels = paste(seq(0,90,10), "-",  seq(10,100,10), "%", sep=""), cex.axis = 0.8)
    } else {
      axis(1, at = 1:10, labels = paste(0, "-",  seq(10,100,10), "%", sep=""), cex.axis = 0.8)
    }
  }
}

