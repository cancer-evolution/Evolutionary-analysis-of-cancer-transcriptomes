#Phy plots all phylostrata

library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape)
library(cowplot)
library(MASS)
library(GSVA)
library(GSEABase)
library(GO.db)
library(gplots)


source("TAI_and_phylostrata_analysis_functions.R")
tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

TAI_phy <- list()
difference_phy_means <- list()
sum_phy_variances <- list()

tumours <- c("LUAD", "LUSC","BRCA", "PRAD", "LIHC", "COAD", "STAD")

normal_expression <- list()
tumour_expression <- list()


for (type in tumours){
  normal_expression[[type]] <- read_expression(type, "normal")
  tumour_expression[[type]] <- read_expression(type, "tumour")
  
  print(paste("Done with", type, "!"))
}

for (type in tumours){
  TAI_phy[[type]]$normal <- as.data.frame(calculate_TAI(normal_expression[[type]], type))
  TAI_phy[[type]]$tumour <- as.data.frame(calculate_TAI(tumour_expression[[type]], type))
  print(type)
  print(nrow(TAI_phy[[type]]$normal))
  print(nrow(TAI_phy[[type]]$tumour))
}


#Calculate TAI using Domazet ages
gene_annotations <- read.csv("Human_ensembl_(GRCh37.p13).txt")
for (type in tumours){
  normal_expression[[type]] <- read_expression_for_Domazet(type, "normal")
  tumour_expression[[type]] <- read_expression_for_Domazet(type, "tumour")
  
  print(paste("Done with", type, "!"))
}

TAI_phy <- list()
for (type in tumours){
  TAI_phy[[type]]$normal <- as.data.frame(calculate_TAI(normal_expression[[type]], type))
  TAI_phy[[type]]$tumour <- as.data.frame(calculate_TAI(tumour_expression[[type]], type))
}


#End of calculate TAI using Domazet ages


corrected_TAI_phy <- list()
difference_phy_means_df <- vector()
sum_phy_variances_df <- vector()
for (tumour_type in tumours){
  corrected_TAI_phy[[tumour_type]] <- loess_correction_plus_median(TAI_phy[[tumour_type]])
  means_normal = apply(corrected_TAI_phy[[tumour_type]]$normal[,3:18], 2, function(x){
    mean(as.numeric(as.character(x)))})
  means_tumour = apply(corrected_TAI_phy[[tumour_type]]$tumour[,3:18], 2,function(x){
    mean(as.numeric(as.character(x)))})
  difference_phy_means_df <- rbind(difference_phy_means_df,
                                   cbind(means_tumour - means_normal))
  var_normal = apply(corrected_TAI_phy[[tumour_type]]$normal[,3:18], 2, function(x){
    var(as.numeric(as.character(x)))})
  var_tumour = apply(corrected_TAI_phy[[tumour_type]]$tumour[,3:18], 2, function(x){
    var(as.numeric(as.character(x)))})
  sum_phy_variances_df <- rbind(sum_phy_variances_df,
                                cbind(means_tumour - means_normal, tumour_type))
  print(tumour_type)
}

##remove cell cycle and recalculate the TAI

#recalculate the TAI

signature_name <- "cell cycle"

load("GOslims_gene_ids.R") 
#list named GOslims with cellular processes as names, and genes as items
signature <- GOslims[[signature_name]]


#cell cycle genes in the TCGA data
sum(normal_expression[[type]][,"GeneID"] %in% signature)

write.table(normal_expression[[type]][,"GeneID"][normal_expression[[type]][,"GeneID"] %in% signature],
            file="Supplementary table 6.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

UC_genes <- normal_expression[[type]][which(normal_expression[[type]]$Phylostratum %in% c(1,2,3)),
                                      "GeneID"]


TAI_phy_ns <- list()
corrected_TAI_phy_ns <- list()

for(type in tumours){
  TAI_phy_ns[[type]][["normal"]] <- remove_signature_TAI(signature, normal_expression[[type]], type)
  TAI_phy_ns[[type]][["tumour"]] <- remove_signature_TAI(signature, tumour_expression[[type]], type)
  corrected_TAI_phy_ns[[type]] <- loess_correction_plus_median(TAI_phy_ns[[type]])
  print(type)
}

corrected_TAI_ns <- vector()
for(tumour_type in tumours){
  print(tumour_type)
  for (tissue in c("normal", "tumour")){
    print(dim(corrected_TAI_phy_ns[[tumour_type]][[tissue]]))
    corrected_TAI_ns <- rbind(corrected_TAI_ns,
                              cbind(as.character(corrected_TAI_phy_ns[[tumour_type]][[tissue]][,2]),
                                    tumour_type, paste(toupper(substr(tissue,1,1)), substr(tissue,2,10), sep="")))    
  }
}

colnames(corrected_TAI_ns) <- c("TAI", "Tumour_type", "Tissue_type")
corrected_TAI_ns <- as.data.frame(corrected_TAI_ns)
corrected_TAI_ns$TAI <- as.numeric(as.character(corrected_TAI_ns$TAI))
corrected_TAI_ns$Tumour_type <- factor(corrected_TAI_ns$Tumour_type, levels=tumours)

p_values <- vector()
for(tumour_type in tumours){
  p_values <- c(p_values, format.pval(wilcox.test(as.numeric(as.character(corrected_TAI_phy_ns[[tumour_type]][["tumour"]][,2])), 
                                                  as.numeric(as.character(corrected_TAI_phy_ns[[tumour_type]][["normal"]][,2])))$p.value, digits=2))
}

ann_text <- data.frame(Tissue_type = rep(c("Normal", "Tumour"),each=7), TAI = rep(c(5,5,5.4,4.5,6.1,4.5,5),each=2),
                       Tissue_type = factor(tumours,levels = tumours))


significance <- data.frame(a = c(1, 1, 2, 2), b = c(6.1, 6.2, 6.2, 6.1))

g_TAI <- plot_TAI(corrected_TAI_ns)

pdf(file="Supp_6.pdf", height=4, width=8)  
g_TAI
dev.off()

#end of cell cycle

phy_diff_gg <- cbind(difference_phy_means_df, sum_phy_variances_df,substr(rownames(difference_phy_means_df), 5,6))
colnames(phy_diff_gg) <- c("Difference", "Variance", "Tumour_type", "Phylostrata")
rownames(phy_diff_gg) <- 1:nrow(phy_diff_gg)
phy_diff_gg <- as.data.frame(phy_diff_gg)
phy_diff_gg$Phylostrata <- factor(phy_diff_gg$Phylostrata, levels=1:16)

phy_diff_gg$Difference <- as.numeric(as.character(phy_diff_gg$Difference))
phy_diff_gg$Variance <- as.numeric(as.character(phy_diff_gg$Variance))
phy_diff_gg$Tumour_type <- factor(phy_diff_gg$Tumour_type, levels=tumours)
phy_diff_gg$U_M <- factor(ifelse(phy_diff_gg$Phylostrata %in% 1:3, "UC", "MC"), levels=c("UC", "MC"))
phy_diff_gg$U <- factor(ifelse(phy_diff_gg$Phylostrata %in% 1:3, "UC", NA), levels=c("UC"))
phy_diff_gg$M <- factor(ifelse(phy_diff_gg$Phylostrata %in% 1:3, NA, "MC"), levels=c("MC"))


#One plot only
gene_age_annotations_UC <- data.frame(a = c(1, 3), b = c(0.055, 0.06, 0.06, 0.055))
gene_age_annotations_MC <- data.frame(a = c(3, 16), b = c(0.055, 0.06, 0.06, 0.055))

g_phy_all <- ggplot(phy_diff_gg, aes(x=Phylostrata, y=Difference))+
  geom_point()+
  #geom_rect(data=NULL,aes(xmin=0.4,xmax=16.6,ymin=-Inf,ymax=0),
  #          fill='lightblue', alpha=0.02)+
  #geom_rect(data=NULL,aes(xmin=0.4,xmax=16.6,ymin=0,ymax=Inf),
  #          fill='lightpink', alpha=0.006)+
  geom_hline(yintercept=0, size=1)+
  annotate("rect", xmin = 0.95, xmax = 3, ymin = 0, ymax = 0.06,
           alpha = .2, fill="pink")+
  annotate("rect", xmin = 3, xmax = 12, ymin = -0.075, ymax = 0,
           alpha = .2, fill="lightblue")+
  annotate("rect", xmin = 12, xmax = 16.05, ymin = -0.075, ymax = 0,
           alpha = .2, fill="yellow")+
  #stat_smooth(aes(group=Tumour_type, color=Tumour_type), # continuous x-axis
  #            se = F, method = "lm", formula = y ~ poly(x, 7)) +
  geom_path(aes(group=Tumour_type, color=Tumour_type))+
  geom_line(data = gene_age_annotations_UC, aes(x = a, y = b))+
  geom_line(data = gene_age_annotations_MC, aes(x = a, y = b))+
  annotate("text", x = 2, y = 0.065, label = "UC genes", size = 7)+
  annotate("text", x = 9, y = 0.065, label = "MC genes", size = 7)+
  scale_color_manual(values=c("red", "blue", "darkgreen", "purple", "deepskyblue", "orange", "black"))+
  annotate("text", x=1:16, y=-0.072, label=c("Cellular organisms", "Eukaryota",
                                             "Opisthokonta", "Metazoa",
                                             "Eumetazoa", "Bilateria",
                                             "Chordata", "Euteleostomi",
                                             "Ammiota", "Mammalia",
                                             "Theria", "Eutheria",
                                             "Euarchontoglires", "Catarrhini",
                                             "Homininae", "Homo sapiens"), size=7,
           angle = 70, hjust = 0)+
  coord_cartesian(ylim = c(-0.068, 0.068), xlim=c(1,16.7))+
  geom_point(size=1) + guides(colour=FALSE)+
  ylab("Difference in\nphylostratum proportion")+
  theme_bw()+
  theme(text = element_text(size=18),axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))

pdf(file="Figure_1c.pdf", height=5.8, width=13)
g_phy_all
dev.off()

#Percentage unicellular and multicellular
p_UC_MC <- vector()
p_UC_MC_sample <- vector()
for (tumour in tumours){
  for (type in c("normal", "tumour")){
    p_UC_MC_sample <- rbind(p_UC_MC_sample,
                            cbind(as.character(TAI_phy[[tumour]][[type]][,1]), 
                                  apply(TAI_phy[[tumour]][[type]][,2:4],1, function(x){
                                    sum(as.numeric(as.character(x)))}), tumour, type))
    t_U <- median(apply(TAI_phy[[tumour]][[type]][,2:4],1, function(x){
      sum(as.numeric(as.character(x)))}))
    t_U_sd <- sd(apply(TAI_phy[[tumour]][[type]][,2:4],1, function(x){
      sum(as.numeric(as.character(x)))}))
  #  t_M <- mean(apply(TAI_phy[[tumour]][[type]][,6:18],1, sum))
    p_UC_MC <- rbind(p_UC_MC, c(t_U, "UC", type, tumour, t_U_sd))
  #  p_UC_MC <- rbind(p_UC_MC, c(t_M, "MC", type, tumour))
  }
} 

p_UC_MC <- as.data.frame(p_UC_MC)
colnames(p_UC_MC) <- c("Percentage", "Age", "Tissue_type", "Tumour_type", "SD")

p_UC_MC_sample <- as.data.frame(p_UC_MC_sample)
colnames(p_UC_MC_sample) <- c("Patient", "Percentage_UC", "Tumour_type", "Tissue_type")

p_UC_MC$Percentage <- as.numeric(as.character(p_UC_MC$Percentage))*100
p_UC_MC$SD <- as.numeric(as.character(p_UC_MC$SD))*100


p_UC_MC_normal_median <- median(subset(p_UC_MC, Tissue_type =="normal")$Percentage)
p_UC_MC_tumour_median <- median(subset(p_UC_MC, Tissue_type =="tumour")$Percentage)

p_UC_MC_median <- data.frame(Median=c(p_UC_MC_normal_median, p_UC_MC_tumour_median), 
                             Tissue_type=c("normal", "tumour"))

ymin <- c(min(subset(p_UC_MC, Tissue_type =="normal")$Percentage),
          min(subset(p_UC_MC, Tissue_type =="tumour")$Percentage))
ymax <-c(max(subset(p_UC_MC, Tissue_type =="normal")$Percentage), 
          max(subset(p_UC_MC, Tissue_type =="tumour")$Percentage))
p_UC_MC_median$ymin <- ymin
p_UC_MC_median$ymax <- ymax

p_UC_MC$Tumour_type <- factor(p_UC_MC$Tumour_type, tumours)

g_per <- ggplot(p_UC_MC, aes(Tissue_type, Percentage))+
  theme_bw()+
  geom_bar(data=p_UC_MC_median, aes(y=Median, x=Tissue_type), stat = "identity", alpha=.2)+
  geom_path(aes(group=Tumour_type, color=Tumour_type), size=.75)+
  scale_color_manual(values=c("red", "blue", "darkgreen", "purple", "deepskyblue", "orange", "black"))+
  geom_errorbar(data=p_UC_MC_median, aes(x=Tissue_type, ymin = ymin, ymax=ymax, y=NULL), width=0, color="grey", size=1)+
  geom_point()+
  coord_cartesian(ylim = c(36, 62))+
  #coord_cartesian(ylim = c(60, 85))+  #DOMAZET
  theme(axis.title.x = element_blank(), legend.title=element_blank(), text = element_text(size=17),
        axis.text.x= element_text(size=17,angle = 45, hjust = 1),
        axis.text.y= element_text(size=17), axis.title.y = element_text(size=17))+
  ylab("Percentage of unicellular\ntranscriptome")+
  scale_x_discrete(labels=c("Normal", "Tumour"))+
  guides(color=guide_legend(title="Tumour types"))

pdf(file="Figure_1b.pdf", height=3.5, width=4)
g_per
dev.off()


#wilcoxon tests
for(type in tumours){
	temp <- p_UC_MC_sample[which(p_UC_MC_sample$Tumour_type == type),]
	temp_n <- as.numeric(as.character(temp[temp$Tissue_type == "normal",2]))
	temp_t <- as.numeric(as.character(temp[temp$Tissue_type == "tumour",2]))
	print(wilcox.test(temp_t, temp_n, alternative="greater")$p.val)
}



#TAI plot
corrected_TAI_df <- vector()

for(tumour_type in tumours){
  for (tissue in c("normal", "tumour")){
    corrected_TAI_df <- rbind(corrected_TAI_df,
                              cbind(as.character(corrected_TAI_phy[[tumour_type]][[tissue]][,2]),
                                    tumour_type, paste(toupper(substr(tissue,1,1)), substr(tissue,2,10), sep="")))    
  }
}
colnames(corrected_TAI_df) <- c("TAI", "Tumour_type", "Tissue_type")
corrected_TAI_df <- as.data.frame(corrected_TAI_df)
corrected_TAI_df$TAI <- as.numeric(as.character(corrected_TAI_df$TAI))
corrected_TAI_df$Tumour_type <- factor(corrected_TAI_df$Tumour_type, levels=tumours)

p_values <- vector()
for(tumour_type in tumours){
  p_values <- c(p_values, format.pval(wilcox.test(as.numeric(as.character(corrected_TAI_phy[[tumour_type]][["tumour"]][,2])), 
              as.numeric(as.character(corrected_TAI_phy[[tumour_type]][["normal"]][,2])))$p.value, digits=2))
}

ann_text <- data.frame(Tissue_type = rep(c("Normal", "Tumour"),each=7), TAI = rep(c(5,5,5.4,4.5,6.1,4.5,5),each=2),
                       Tissue_type = factor(tumours,levels = tumours))


significance <- data.frame(a = c(1, 1, 2, 2), b = c(6.1, 6.2, 6.2, 6.1))

g_TAI <- plot_TAI(corrected_TAI_df)

pdf(file="Figure_1a.pdf", height=4, width=8)  
g_TAI
dev.off()


# #TAI at different ranges

number_of_bins <- 10

order_of_genes <- list()
for (type in tumours){
  total_expression <- cbind(normal_expression[[type]], tumour_expression[[type]][,3:ncol(tumour_expression[[type]])])

  #total_expression <- cbind(tumour_expression[[type]])

  total_expression_corrected_library_size <- apply(total_expression[,3:ncol(total_expression)],2,function(vector){ return(vector/sum(vector))})
  rownames(total_expression_corrected_library_size) <- total_expression[,1]
  
  #Removing genes = 0 in either all normal or tumour samples
  median_normal <- apply(normal_expression[[type]][3:ncol(normal_expression[[type]])],1,median)
  median_tumour <- apply(tumour_expression[[type]][3:ncol(tumour_expression[[type]])],1,median)
  
  remove_tumour <- which(median_normal == 0)
  remove_normal <- which(median_tumour == 0)
  
  remove_total <- sort(unique(c(remove_normal, remove_tumour)))
  
  total_expression_corrected_library_size <- total_expression_corrected_library_size[-remove_total,]
  
  median_expression <- apply(total_expression_corrected_library_size, 1, median)
  names(median_expression) <- rownames(total_expression_corrected_library_size)
  median_expression <- median_expression[which(median_expression != 0)]
  order_of_genes[[type]] <- names(sort(median_expression, decreasing=FALSE))
}

#Distribution by ranges of the percentage of UC
per_UC_binned <- list()
for(type in tumours){
  per_UC_binned[[type]] <- calculate_percentage_UC_bins(normal_expression[[type]], tumour_expression[[type]], number_of_bins, "specific_ranges", order_of_genes[[type]])
  print(type)
}

per_UC_p_val <- vector()
for(type in tumours){
  per_UC_p_val <- rbind(per_UC_p_val, cbind(per_UC_binned[[type]]$p_values_ts,
                                            per_UC_binned[[type]]$p_values_greater,
                                            per_UC_binned[[type]]$p_values_less))
}

colnames(per_UC_p_val) <- c("Two-sided", "Greater", "Less")
rownames(per_UC_p_val) <- rep(tumours, each=10)
per_UC_p_val <- cbind(per_UC_p_val, paste("Bin", 1:10, sep="_"))

bin_vector <- c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")

per_UC_binned_df <- vector()
for(type in tumours){
  t1 <-  cbind(type, "Normal", cbind(bin_vector, per_UC_binned[[type]]$normal_median, per_UC_binned[[type]]$normal_min, per_UC_binned[[type]]$normal_max, per_UC_binned[[type]]$normal_sd))
  t2 <-  cbind(type, "Tumour", cbind(bin_vector, per_UC_binned[[type]]$tumour_median, per_UC_binned[[type]]$tumour_min, per_UC_binned[[type]]$tumour_max, per_UC_binned[[type]]$tumour_sd))
  per_UC_binned_df <- rbind(per_UC_binned_df,t1, t2)
}           
colnames(per_UC_binned_df) <- c("Tumour_type", "Tissue_type", "Bin", "Median_Percentage_UC", "Min_Percentage_UC", "Max_Percentage_UC", "SD_Percentage_UC")

per_UC_binned_df <- as.data.frame(per_UC_binned_df)

per_UC_binned_df$Median_Percentage_UC <- as.numeric(as.character(per_UC_binned_df$Median_Percentage_UC))

per_UC_binned_df$Tumour_type <- factor(per_UC_binned_df$Tumour_type, levels=tumours)

per_UC_binned_df$Bin <- factor(per_UC_binned_df$Bin, levels=bin_vector)
per_UC_binned_df$Tissue_type <- factor(per_UC_binned_df$Tissue_type, levels=c("Tumour", "Normal"))

per_UC_binned_df$Min_Percentage_UC <- as.numeric(as.character(per_UC_binned_df$Min_Percentage_UC))
per_UC_binned_df$Max_Percentage_UC <- as.numeric(as.character(per_UC_binned_df$Max_Percentage_UC))

per_UC_binned_df$SD_Percentage_UC <- as.numeric(as.character(per_UC_binned_df$SD_Percentage_UC))
per_UC_binned_df$SD_Percentage_UC <- as.numeric(as.character(per_UC_binned_df$SD_Percentage_UC))

pdf("Supp_7.pdf", width=10, height=5)
ggplot(per_UC_binned_df, aes(Bin, Median_Percentage_UC))+
  theme_bw()+
  geom_point(aes(colour=Tissue_type))+
  geom_line(aes(group=Tissue_type, colour=Tissue_type)) +
  #geom_path(aes(group=Tumour_type, color=Tumour_type), size=.75)+
  # scale_color_manual(values=c("red", "blue", "darkgreen", "purple", "deepskyblue", "orange", "black"))+
  geom_errorbar(aes(ymin = Median_Percentage_UC-SD_Percentage_UC, ymax=Median_Percentage_UC+SD_Percentage_UC, color=Tissue_type), width=0.25, position = position_dodge(width = 0.3))+
  ylab("Percent unicellular transcripts")+
  xlab("Expression decile")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Tumour_type, nrow=2)
dev.off()


#Density plot
tumours <- c("LUAD", "LUSC","BRCA", "PRAD", "LIHC", "COAD", "STAD")
tT_df <- vector()
gene_phy <- read.csv("geneIDs_entrez_phylostrata.txt")
#File with Gene ID, Entrez and Phylostrata columns

for (type in tumours){
  load(paste(type, "_DEG_tT_all.RData", sep="")) #object resulting from limma analysis, showing the gene IDs/logFC/AveExpr/T/P.Value/adj.P.Val/B for all genes
  t <- expression_tT_raw
  tT_df <- rbind(tT_df,
                 cbind(t[,2], gene_phy[match(rownames(t), gene_phy$Entrez), 3], type))
  
}
colnames(tT_df) <- c("logFC", "Genes", "Tumour_type")
tT_df <- as.data.frame(tT_df)
tT_df <- tT_df[!is.na(tT_df$Genes),]
tT_df$Genes <- ifelse(tT_df$Genes %in% 1:3, "UC", "MC")
tT_df$Genes <- factor(tT_df$Genes, levels=c("UC", "MC"))
tT_df$logFC <- as.numeric(as.character(tT_df$logFC))
tT_df$Tumour_type <- factor(tT_df$Tumour_type, levels=tumours)

g_dens <- ggplot(tT_df, aes(x=logFC, fill=Genes))+
  geom_density(alpha=.5, size=.5)+
  xlim(c(-3,2))+
  facet_grid(.~Tumour_type)+
  theme_bw()

pdf("Supp_8.pdf", width=10, height=5)
g_dens
dev.off()

for(tumour in tumours){
  local <- tT_df[which(tT_df$Tumour_type == tumour),]
  local_UC <- local[which(local$Genes == "UC"),"logFC"]
  local_MC <- local[which(local$Genes == "MC"),"logFC"]
  print(ks.test(local_UC, local_MC)$p.value)
}

#grid.arrange(g_phy, g_dens, g_TAI, ncol=1)


#TAI and differentiation
tumours_diff <- c("PRAD", "LIHC", "STAD")

TAI <- vector()
clinical_data <- data.frame()
for (tumour_type in tumours_diff){
  
  TAI <- rbind(TAI, cbind(read.delim(paste(tumour_type, "_corrected_TAI_phy_normal.txt", sep="")), #File with patient ID and TAI score (normal)
                      Tissue_type="Normal", Tumour_type=tumour_type),
               cbind(read.delim(paste(tumour_type, "_corrected_TAI_phy_tumour.txt", sep="")),  #File with patient ID and TAI score (normal)
                     Tissue_type="Tumour", Tumour_type=tumour_type))
  
  clinical_data_TCGA <- read.delim(paste("nationwidechildrens.org_clinical_patient_", tumour_type, ".txt", sep=""))  #Metadata files from TCGA
  
  colnames(clinical_data_TCGA) <- unlist(clinical_data_TCGA[1,])
  clinical_data_TCGA <- clinical_data_TCGA[-c(1,2),]
  patients_clinical <- data.frame(do.call('rbind', strsplit(as.character(clinical_data_TCGA$bcr_patient_barcode),'-',fixed=TRUE)))$X3
  
  if ("neoplasm_histologic_grade" %in% colnames(clinical_data_TCGA))
  {
    clinical_data <- rbind(clinical_data, data.frame(Patient=patients_clinical,
                                                     Gleason=rep(NA, length(patients_clinical)),
                                                     Grade=paste("Grade", substr(clinical_data_TCGA$neoplasm_histologic_grade,2,2)), Tumour_type = tumour_type))
  } else if ("gleason_score" %in% colnames(clinical_data_TCGA)){
    clinical_data <- rbind(clinical_data, data.frame(Patient=patients_clinical,
                                                     Gleason=paste("Gleason", clinical_data_TCGA$gleason_score),
                                                     Grade=rep(NA, length(patients_clinical)), Tumour_type = tumour_type))
  }
} 
colnames(TAI)[2] <-"TAI"
TAI_clinical <- cbind(TAI, clinical_data[match(interaction(TAI$Patient, TAI$Tumour_type),
                                               interaction(clinical_data$Patient, clinical_data$Tumour_type)), 2:3])

TAI_clinical <- TAI_clinical[!(TAI_clinical$Grade %in% c("[Not Available]", "GX")),]

TAI_clinical_glea <- TAI_clinical[which(TAI_clinical$Tumour_type == "PRAD"),]
TAI_clinical_glea$Gleason <- factor(TAI_clinical_glea$Gleason, 
                                    levels=c("Normal", "Gleason 6", "Gleason 7", "Gleason 8", "Gleason 9", "Gleason 10"))
TAI_clinical_glea$Gleason[TAI_clinical_glea$Tissue_type == "Normal"] <- "Normal"
TAI_clinical_glea <- TAI_clinical_glea[!is.na(TAI_clinical_glea$Gleason),]

g_glea <- ggplot(TAI_clinical_glea, aes(x=Gleason, y=TAI, fill=Gleason))+
          geom_boxplot()+
          geom_point(position = position_jitter(width = 0.2)) +
          guides(fill=FALSE)+
          #facet_grid(.~Tumour_type)+
          theme_bw()+
          theme(axis.title.x = element_blank(),
          strip.background = element_rect(colour="black", fill="white"), 
          strip.text.x = element_text(size=15, face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=19),
          axis.text.y= element_text(size=17),
          text = element_text(size=18),
          axis.title.y = element_text(size=17))

pdf("Figure_1d.pdf", height=3)
g_glea
dev.off()

#Jonckheere-Terpstra test to check decreasing trend
library(clinfun)

TAI_clinical_glea$Gleason <- factor(TAI_clinical_glea$Gleason, levels=c("Normal", "Gleason 6", "Gleason 7", "Gleason 8", "Gleason 9", "Gleason 10"), ordered=TRUE)
jonckheere.test(x=TAI_clinical_glea$TAI, g=TAI_clinical_glea$Gleason, alternative="decreasing")

Glea_6 <- TAI_clinical_glea[TAI_clinical_glea$Gleason == "Gleason 6",]
Glea_8_9_10 <- TAI_clinical_glea[TAI_clinical_glea$Gleason %in% c("Gleason 8", "Gleason 9", "Gleason 10"),]

wilcox.test(Glea_8_9_10$TAI,Glea_6$TAI,alternative="less")

boxplot(Glea_6$TAI, Glea_8_9_10$TAI)

p <- polr(Gleason ~ TAI, data=TAI_clinical_glea, Hess=TRUE)
ctable <- coef(summary(p))
p_value <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2



TAI_clinical_grade <- TAI_clinical[which(TAI_clinical$Tumour_type != "PRAD"),]

TAI_clinical_grade$Grade <- factor(TAI_clinical_grade$Grade, 
                                    levels=c("Normal", "Grade 1", "Grade 2", "Grade 3", "Grade 4"))
TAI_clinical_grade$Grade[TAI_clinical_grade$Tissue_type == "Normal"] <- "Normal"
TAI_clinical_grade <- TAI_clinical_grade[!is.na(TAI_clinical_grade$Grade),]


g_grade_l <- ggplot(subset(TAI_clinical_grade, Tumour_type=="LIHC"), aes(x=Grade, y=TAI, fill=Grade))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2)) +
  facet_grid(.~Tumour_type, scales="free")+
  guides(fill=FALSE)+
  scale_y_continuous("TAI")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=10, face="bold"))




g_grade_s <- ggplot(subset(TAI_clinical_grade, Tumour_type=="STAD"), aes(x=Grade, y=TAI, fill=Grade))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2)) +
  facet_grid(.~Tumour_type, scales="free")+
  guides(fill=FALSE)+
  scale_y_continuous("TAI")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=10, face="bold"))

pdf("Supp_9.pdf", height=5, width=10)
grid.arrange(g_grade_l, g_grade_s, ncol=2)
dev.off()

TAI_clinical_grade$Grade <- factor(TAI_clinical_grade$Grade, levels=c("Normal", "Grade 1", "Grade 2", "Grade 3", "Grade 4"), ordered=TRUE)

LIHC_grade <- subset(TAI_clinical_grade, Tumour_type=="LIHC")
jonckheere.test(x=LIHC_grade$TAI, g=LIHC_grade$Grade, alternative="decreasing")

STAD_grade <- subset(TAI_clinical_grade, Tumour_type=="STAD")
jonckheere.test(x=STAD_grade$TAI, g=STAD_grade$Grade, alternative="decreasing")


#Clinical data and percentage UC
p_UC_clinical <- data.frame(p_UC_MC_sample[match(clinical_data$Patient,
                                                 p_UC_MC_sample[,1]),],
                  clinical_data[,c("Gleason", "Grade")])          

p_UC_clinical <- p_UC_clinical[!(p_UC_clinical$Grade %in% c("Grade X", "Grade N")),] 

p_UC_clinical$Gleason <- factor(p_UC_clinical$Gleason, levels=c("Normal", "Gleason 6", "Gleason 7", "Gleason 8",
                                                                "Gleason 9", "Gleason 10"))
p_UC_clinical$Grade <- factor(p_UC_clinical$Grade, levels=c("Normal", "Grade 1", "Grade 2", "Grade 3",
                                                                "Grade 4"))

p_UC_clinical[which(p_UC_clinical$Tissue_type == "normal"),c("Gleason", "Grade")] <- "Normal" 

p_UC_clinical <- p_UC_clinical[!is.na(p_UC_clinical$Patient),]

p_UC_glea <- p_UC_clinical[p_UC_clinical$Tumour_type == "PRAD",c(1:5)]
p_UC_glea$Percentage_UC <- as.numeric(as.character(p_UC_glea$Percentage_UC))

ggplot(p_UC_glea, aes(x=Gleason, y=Percentage_UC, fill=Gleason))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2)) +
  guides(fill=FALSE)+
  facet_grid(.~Tumour_type)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=10, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_UC_grade <- p_UC_clinical[p_UC_clinical$Tumour_type != "PRAD",c(1:4, 6)]
p_UC_grade$Percentage_UC <- as.numeric(as.character(p_UC_grade$Percentage_UC))

ggplot(p_UC_grade, aes(x=Grade, y=Percentage_UC, fill=Grade))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2)) +
  guides(fill=FALSE)+
  facet_grid(.~Tumour_type)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=10, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))


#Clone number and proliferation marker
CNV_clonality_table <- read.delim("CNV_clonality_table.txt")  #TCGA metadata
tumours <- c("LUAD", "LUSC", "BRCA","PRAD","LIHC", "COAD", "STAD") 

CNV_clonality_table[,1] <- substr(CNV_clonality_table[,1],9,12)

TAI_t <- vector()

for (tumour_type in tumours){
  
  TAI_t <- rbind(TAI_t, cbind(read.delim(paste(tumour_type, "_corrected_TAI_phy_tumour.txt", sep="")),  #table of patient ID, TAI
                     Tissue_type="Tumour", Tumour_type=tumour_type))
}
colnames(TAI_t)[2] <- "TAI"

TAI_pro <- data.frame(TAI_t, CNV_clonality_table[match(TAI_t$Patient, CNV_clonality_table$TCGA.ID), "MKI67_mRNAexpression"])
colnames(TAI_pro)[5] <- "MKI67_mRNAexpression"
TAI_pro <- TAI_pro[!is.na(TAI_pro$MKI67_mRNAexpression),]

cor_pro <- vector()
cor_test <- vector()
for (tumour_type in c("LUAD", "LUSC", "PRAD")){
  cor_pro <- c(cor_pro, round(cor(subset(TAI_pro,  Tumour_type==tumour_type)$TAI, subset(TAI_pro,  Tumour_type==tumour_type)$MKI67_mRNAexpression, method="sp"),3))
  cor_test <- c(cor_test, cor.test(subset(TAI_pro,  Tumour_type==tumour_type)$TAI, subset(TAI_pro,  Tumour_type==tumour_type)$MKI67_mRNAexpression, method="sp")$p.value)
  
}
cor_pro <- paste("Cor:", cor_pro)
g_pro <- ggplot(TAI_pro, aes(x=TAI, y=MKI67_mRNAexpression))+
  geom_point()+
  facet_grid(.~Tumour_type)+
  annotate("text", x=4.2, y=9, label=cor_pro)+
  theme_bw()

g_pro_p <- ggplot(subset(TAI_pro, Tumour_type=="PRAD"), aes(y=TAI, x=MKI67_mRNAexpression))+
  theme_bw()+
  geom_point()+
  xlab("MKI67 mRNA expression")+
  #facet_grid(.~Tumour_type)+
  annotate("text", y=4, x=7.3, label=cor_pro[3], size=7)+
  theme(strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=15, face="bold"),
        text = element_text(size=18),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20))
  
pdf("Figure_1e.pdf", height=3, width=5)
g_pro_p
dev.off()

g_pro_other <- ggplot(subset(TAI_pro, Tumour_type!="PRAD"), aes(y=TAI, x=MKI67_mRNAexpression))+
  theme_bw()+
  geom_point()+
  xlab("MKI67 mRNA expression")+
  facet_grid(.~Tumour_type)+
  annotate("text", x=8.5, y=4.3, label=cor_pro[1:2], size=4)+
  theme(strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size=10, face="bold"))

pdf("Supp_10.pdf", height=5, width=10)
g_pro_other
dev.off()


TAI_clone <- data.frame(TAI_t, CNV_clonality_table[match(TAI_t$Patient, CNV_clonality_table$TCGA.ID), "CloneNumber.PurityNormalized."])
colnames(TAI_clone)[5] <- "Clone_number"
TAI_clone <- TAI_clone[!is.na(TAI_clone$Clone_number),]
cor_clone <- vector()
for (tumour_type in c("LUAD", "LUSC", "PRAD", "STAD")){
  cor_clone <- c(cor_clone, round(cor(subset(TAI_clone,  Tumour_type==tumour_type)$TAI, subset(TAI_clone,  Tumour_type==tumour_type)$Clone_number, method="sp"),3))
}
cor_clone <- paste("Cor:", cor_clone)
g_clone <- ggplot(TAI_clone, aes(x=TAI, y=Clone_number))+
  geom_point()+
  facet_grid(.~Tumour_type)+
  annotate("text", x=4.2, y=17, label=cor_clone)+
  theme_bw()

TAI_CNV <- data.frame(TAI_t, CNV_clonality_table[match(TAI_t$Patient, CNV_clonality_table$TCGA.ID), "lowORhigh_cnvAbundance"])
colnames(TAI_CNV)[5] <- "CNV_abundance"
TAI_CNV <- TAI_CNV[!is.na(TAI_CNV$CNV_abundance),]
cor_CNV <- vector()
for (tumour_type in c("LUAD", "LUSC", "PRAD", "STAD")){
  cor_CNV <- c(cor_CNV, round(cor(subset(TAI_CNV,  Tumour_type==tumour_type)$TAI, subset(TAI_CNV,  Tumour_type==tumour_type)$CNV_abundance, method="sp"),3))
}
cor_CNV <- paste("Cor:", cor_CNV)
g_CNV <- ggplot(TAI_CNV, aes(x=TAI, y=CNV_abundance))+
  geom_point()+
  facet_grid(.~Tumour_type)+
  annotate("text", x=4.2, y=-6.5, label=cor_clone)+
  theme_bw()

grid.arrange(g_pro, g_clone, g_CNV, ncol=1)

