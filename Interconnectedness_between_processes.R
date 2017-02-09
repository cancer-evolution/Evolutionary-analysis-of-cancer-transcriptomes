#Overlap and connections between pathways
library(GSEABase)
library(GO.db)
library(mGSZ)
library(igraph)
library(ggplot2)
def.par <- par(no.readonly = TRUE)
library("gridExtra")
library(reshape2)
library(gplots)
library(scales)
library(ggrepel)
library(clinfun)

source("Interconnectedness_between_processes_functions.R")
pathway_database <- "GOslims"
network_database <- "PathwayCommons"
#network_database <- "HPRD"


pathway_ages <- #Data frame of GOslim and GOslim ages


###Excluding the interactions found within a pathway from the measure of the Degree_of_interactions
if(network_database == "FI"){
  network <- read.delim("FIsInGene_121514_with_annotations.txt")  #Downloaded from Reactome's website
  network <- unique(network[,c(1,2)])
  gene_annotations <- #Data frame with columns: Ensembl.Gene.ID,Ensembl.Protein.ID,Associated.Gene.Name,Associated.Gene.DB,EntrezGene.ID
  network[,1] <- gene_annotations[match(network[,1], gene_annotations$Associated.Gene.Name), "EntrezGene.ID"]  
  network[,2] <- gene_annotations[match(network[,2], gene_annotations$Associated.Gene.Name), "EntrezGene.ID"]  
  network <- network[!is.na(network[,1]),]
  network <- network[!is.na(network[,2]),]
} else if (network_database == "biogrid"){
  network <- read.delim("BIOGRID-ORGANISM-Homo_sapiens-3.4.128.tab2.txt")  #Downloaded from BioGrid's website
  network <- network[,c(2,3)]
  colnames(network) <- c("Gene1", "Gene2")  
} else if (network_database == "HPRD"){
  network <- read.delim("BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", header=FALSE) #Downloaded from the HPRD's website
  network <- network[,c(2,5)]
  colnames(network) <- c("Gene1", "Gene2")     
} else if (network_database == "PathwayCommons"){
  network <- read.delim("Pathway Commons.7.All.BINARY_SIF.hgnc.sif", header=FALSE, dec=",")  #Downloaded from Pathway Commons' website
  gene_annotations <-  #Data frame with columns: Ensembl.Gene.ID,Ensembl.Protein.ID,Associated.Gene.Name,Associated.Gene.DB,EntrezGene.ID
  network[,1] <- gene_annotations[match(network[,1], gene_annotations[,3]),5]
  network[,3] <- gene_annotations[match(network[,3], gene_annotations[,3]),5]
  network <- network[!is.na(network[,1]),]
  network <- network[!is.na(network[,3]),]
  network <- network[,c(1,3)]  
  colnames(network) <- c("Gene1", "Gene2")
  rownames(network) <- NULL
  network <- unique(network)
}

connection_distance_1_between_pathways2 <- calculate_degree_interaction(network)

connection_distance_1_between_pathways2 <- add_ages(connection_distance_1_between_pathways2)

#Leave only one-sided interactions (A-B only, not A-B and B-A)
t <- as.data.frame(t(apply(connection_distance_1_between_pathways2[,c(1,2)], 1, function(x){
                                                               as.character(sort(x))
                                                               })))
t <- unique(t)
connection_distance_1_between_pathways2 <- connection_distance_1_between_pathways2[match(interaction(t), interaction(connection_distance_1_between_pathways2[,c(1,2)])),]

connection_distance_1_between_pathways2_for_print <- connection_distance_1_between_pathways2[,c("Pathway1", "Pathway2", "Degree_of_interaction")]
colnames(connection_distance_1_between_pathways2_for_print) <- c("Pathway1", "Pathway2", "Interconnectedness")
write.table(connection_distance_1_between_pathways2_for_print, file="Supplementary table 5.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

#Heatmap of interconnectedness
pdf("Supp_15.pdf", width=12, height=9.5)
plot_heatmap_of_interconnectedness(connection_distance_1_between_pathways2)
dev.off()


#Distribution of genes across GOslims

GOslims <- #List, with names GOslims and items gene entrez

pdf("Supp_11.pdf")
barplot(sort(table(unname(unlist(GOslims)))), xlab="Genes", ylab="Number of GOslims in which a gene is found", ylim=c(0,25), xaxt='n')
dev.off()

#Choose the top 10% of slims with the strongest interactions
#Find cutoffs for each:
interactions <- list()
cutoff <- list()
for(type in c("UC_UC", "UC_MC", "MC_MC")){
  interactions[[type]] <- connection_distance_1_between_pathways2[which(connection_distance_1_between_pathways2$Interaction_type == type),]
  interactions[[type]]$"Degree_of_interaction" <- as.numeric(as.character(interactions[[type]]$"Degree_of_interaction"))
  interactions[[type]] <- interactions[[type]][order(interactions[[type]]$"Degree_of_interaction", decreasing=TRUE),]
  interactions[[type]] <- interactions[[type]][!is.nan(interactions[[type]]$Degree_of_interaction),]
  print(type)
  cutoff[[type]] <- interactions[[type]]$"Degree_of_interaction"[floor(nrow(interactions[[type]])*0.1)]
  print(cutoff[[type]])
}


t <- connection_distance_1_between_pathways2
t$Degree_of_interaction <- as.numeric(as.character(t$Degree_of_interaction))
t <- t[!is.na(t$Interaction_type),]

pdf("Figure_3b.pdf", height=5)
ggplot(t, aes(Degree_of_interaction, fill = Interaction_type))+
  geom_density(alpha = 0.5)+
  scale_x_continuous(limits = c(0, 0.25))+
  #scale_x_continuous(limits = c(0, 0.01))+
  xlab("Interconnectedness")+
  theme_bw()
dev.off()


#Compare median interaction of interactions

median(t[which(t$Interaction_type == "UC_UC"),"Degree_of_interaction"], na.rm=TRUE)
median(t[which(t$Interaction_type == "UC_MC"),"Degree_of_interaction"], na.rm=TRUE)
median(t[which(t$Interaction_type == "MC_MC"),"Degree_of_interaction"], na.rm=TRUE)


tumours <- c("LUAD", "LUSC" ,"BRCA", "PRAD", "LIHC", "COAD", "STAD")

ssGSEA_results <- load_ssGSEA_object(pathway_database, paired=FALSE)


cor_ssGSEA <- calculate_correlation_of_ssGSEA_scores(ssGSEA_results, tumours)

cor_number_edges <- list()

for(tumour in tumours){
  cor_number_edges[[tumour]] <- data.frame(connection_distance_1_between_pathways2[,c(1:3,6:8)], 
                                           cor_ssGSEA[[tumour]][["normal"]][match(interaction(connection_distance_1_between_pathways2[,c(1,2)]), interaction(cor_ssGSEA[[tumour]][["normal"]][,c(1,2)])),3],
                                           cor_ssGSEA[[tumour]][["tumour"]][match(interaction(connection_distance_1_between_pathways2[,c(1,2)]), interaction(cor_ssGSEA[[tumour]][["tumour"]][,c(1,2)])),3])
  cor_number_edges[[tumour]] <- cbind(cor_number_edges[[tumour]], as.character(pathway_ages[match(cor_number_edges[[tumour]][,1], pathway_ages[,1]),2]))
  cor_number_edges[[tumour]] <- cbind(cor_number_edges[[tumour]], as.character(pathway_ages[match(cor_number_edges[[tumour]][,2], pathway_ages[,1]),2]))
  cor_number_edges[[tumour]] <- cbind(cor_number_edges[[tumour]], paste(cor_number_edges[[tumour]][,9], cor_number_edges[[tumour]][,10], sep="_"))
  cor_number_edges[[tumour]] <- as.data.frame(cor_number_edges[[tumour]])
  cor_number_edges[[tumour]][,6] <- as.numeric(as.character(cor_number_edges[[tumour]][,6]))
  cor_number_edges[[tumour]][,7] <- as.numeric(as.character(cor_number_edges[[tumour]][,7]))
  cor_number_edges[[tumour]][,8] <- as.numeric(as.character(cor_number_edges[[tumour]][,8]))
  colnames(cor_number_edges[[tumour]]) <- c("Pathway1", "Pathway2", "Number_edges", "Number_genes_pathway_1", "Number_genes_pathway_2", "Degree_of_interaction", "Correlation_normal", "Correlation_tumour", "Age_Pathway1", "Age_Pathway2", "Interaction_type")
  cor_number_edges[[tumour]]$"Interaction_type" <- replace(cor_number_edges[[tumour]]$"Interaction_type", which(cor_number_edges[[tumour]]$"Interaction_type" == "MC_UC"), "UC_MC")
  cor_number_edges[[tumour]] <- cor_number_edges[[tumour]][!is.na(cor_number_edges[[tumour]][,9]),]
  cor_number_edges[[tumour]] <- cor_number_edges[[tumour]][!is.na(cor_number_edges[[tumour]][,10]),]
  cor_number_edges[[tumour]] <- cor_number_edges[[tumour]][intersect_edge_type_with_cutoff(cor_number_edges[[tumour]], "all"),]
}

for(tumour in tumours){
  print(cor_number_edges[[tumour]][intersect(which(cor_number_edges[[tumour]]$Pathway1 == "cell junction organization"),
                                       which(cor_number_edges[[tumour]]$Pathway2 == "chromosome organization")), c("Correlation_normal", "Correlation_tumour")])
}

if(network_database == "PathwayCommons"){
  #Barplot of correlation between processes
  bar_df <- vector()
  for(tumour in tumours){
    t <- cor_number_edges[[tumour]][intersect(which(cor_number_edges[[tumour]]$Pathway1 == "cell death"),
                                              which(cor_number_edges[[tumour]]$Pathway2 == "cellular component assembly")),]
    t1 <- c(t$Correlation_normal, "Normal", tumour)
    t2 <- c(t$Correlation_tumour, "Tumour", tumour)
    
    bar_df <- rbind(bar_df, t1, t2)
    
  }
  
  bar_df <- as.data.frame(bar_df)
  colnames(bar_df) <- c("Correlation", "Type", "Tumour")
  rownames(bar_df) <- 1:nrow(bar_df)
  
  bar_df$Correlation <- as.numeric(as.character(bar_df$Correlation))
  bar_df$Tumour <- factor(bar_df$Tumour, levels=c("LUSC", "STAD", "BRCA", "LUAD", "COAD", "PRAD", "LIHC"))
  
  pdf("Supp_19.pdf", height=5)
  ggplot(bar_df, aes(x=Tumour, y=Correlation, fill=Type))+
    geom_bar(stat="identity",position="dodge")+
    xlab("")+
    #ggtitle("Correlation of expression between cell death\n and cellular component assembly")+
    theme_bw()+
    theme(legend.title=element_blank())+
    geom_hline(aes(yintercept=0))
  dev.off()
}

#Median plot
median_cor_number_edges <- vector()
for(tumour in tumours){
  median_cor_number_edges <- rbind(median_cor_number_edges, cor_number_edges[[tumour]])
}
median_cor_number_edges_n <- aggregate(Correlation_normal ~ Pathway1+Pathway2+Degree_of_interaction+Interaction_type, median_cor_number_edges, median)
median_cor_number_edges_t <- aggregate(Correlation_tumour ~ Pathway1+Pathway2+Degree_of_interaction+Interaction_type, median_cor_number_edges, median)
median_cor_number_edges <- data.frame(median_cor_number_edges_n, Correlation_tumour = median_cor_number_edges_t$Correlation_tumour)

median_cor_number_edges[intersect(which(median_cor_number_edges$Pathway1 == "cell junction organization"),
                                  which(median_cor_number_edges$Pathway2 == "chromosome organization")),]

pdf("Median_interactions_between_pathways_PC.pdf", height=7, width=12)
  all_nodes <- unique(as.vector(c(as.character(cor_number_edges[[tumour]][,1]), as.character(cor_number_edges[[tumour]][,2]))))
  for (type in c("normal", "tumour")){
    if(type == "normal"){
      edge_table <- median_cor_number_edges_n[,c(1,2,3,5,4)]
    }else{
      edge_table <- median_cor_number_edges_t[,c(1,2,3,5,4)]
    }
    plot_circle_of_correlations_PNAS(edge_table, "tumour-median", type, all_nodes, 0.001, "qusage")
    legend(x=-1.3,y=-1.35, c("UC-UC interaction","UC-MC interaction", "MC-MC interaction"),
           lty=1, col=c("red", "green", "blue"),  bg="black", text.col="white", cex = 0.7)
    legend(x=-.6,y=-1.35, c("Positive correlation","Negative correlation"),
           lty=1:2, col="white", bg="black", text.col="white", cex = 0.7)
    legend(x=.2,y=-1.35, c("Upregulated","Downregulated"),
           fill=c("red", "blue"), border="white", bty="o", bg="black", text.col="white", cex = 0.7)
    legend(x=.7,y=-1.35, c("UC process","MC process"),pt.cex=1.7,
           col="white", pch=21:22, border="white", bg="black", text.col="white", cex = 0.7)
  }
dev.off()



#Median
wilcox.test(median_cor_number_edges_t$Correlation_tumour[which(median_cor_number_edges_t$Interaction_type == "UC_MC")],
            median_cor_number_edges_n$Correlation_normal[which(median_cor_number_edges_n$Interaction_type == "UC_MC")], alternative="less")

#Fisher test to see whether positive or negative edges are enriched in tumours in each edge type
median_cor_number_edges <- data.frame(median_cor_number_edges_n, Correlation_tumour = median_cor_number_edges_t$Correlation_tumour)
median_cor_number_edges$Direction_correlation_normal <- ifelse(median_cor_number_edges$Correlation_normal>0, "Positive", "Negative")
median_cor_number_edges$Direction_correlation_tumour <- ifelse(median_cor_number_edges$Correlation_tumour>0, "Positive", "Negative")

for(edge_type in c("UC_UC", "UC_MC", "MC_MC")){
  print(edge_type)
  temp <- median_cor_number_edges[which(median_cor_number_edges$Interaction_type == edge_type),]
  normal_positive <- sum(temp$Direction_correlation_normal == "Positive")
  normal_negative <- sum(temp$Direction_correlation_normal == "Negative")
  tumour_positive <- sum(temp$Direction_correlation_tumour == "Positive")
  tumour_negative <- sum(temp$Direction_correlation_tumour == "Negative")    
    
  a <- fisher.test(cbind(c(tumour_positive, normal_positive), c(tumour_negative, normal_negative)), alternative="greater")$p.val
  b <- fisher.test(cbind(c(tumour_negative, normal_negative),c(tumour_positive, normal_positive)), alternative="greater")$p.val
    
  print(paste("Greater positive in tumour", a))
  print(paste("Greater negative in tumour", b))
}


#UC-MC edges that flip

ks.test(median_cor_number_edges[median_cor_number_edges$Interaction_type == "UC_UC","Correlation_normal"],
        median_cor_number_edges[median_cor_number_edges$Interaction_type == "UC_UC","Correlation_tumour"], alternative="greater")

ks.test(median_cor_number_edges[median_cor_number_edges$Interaction_type == "MC_MC","Correlation_normal"],
        median_cor_number_edges[median_cor_number_edges$Interaction_type == "MC_MC","Correlation_tumour"], alternative="greater")

ks.test(median_cor_number_edges[median_cor_number_edges$Interaction_type == "UC_MC","Correlation_normal"],
        median_cor_number_edges[median_cor_number_edges$Interaction_type == "UC_MC","Correlation_tumour"], alternative="less")

median_cor_number_edges_melt <- melt(median_cor_number_edges[,c("Pathway1", "Pathway2", "Correlation_normal", "Correlation_tumour", "Interaction_type")])

colnames(median_cor_number_edges_melt)[4:5] <- c("Tissue","Correlation")
median_cor_number_edges_melt$Tissue <- as.character(median_cor_number_edges_melt$Tissue)
median_cor_number_edges_melt$Tissue <- replace(median_cor_number_edges_melt$Tissue,
                                               median_cor_number_edges_melt$Tissue == "Correlation_normal",
                                               "Normal")
median_cor_number_edges_melt$Tissue <- replace(median_cor_number_edges_melt$Tissue,
                                               median_cor_number_edges_melt$Tissue == "Correlation_tumour",
                                               "Tumour")

median_cor_number_edges_melt$Interaction_type <- factor(median_cor_number_edges_melt$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))

pdf("Supp_16.pdf", height=5)
ggplot(median_cor_number_edges_melt, aes(x=Correlation))+
  geom_density(aes(fill=Tissue), alpha=.5)+
  facet_grid(~Interaction_type)+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme_bw()
dev.off()


median_cor_number_edges$Interaction_type <- factor(median_cor_number_edges$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))
coeff <- coef(lm(Correlation_tumour ~ Correlation_normal, data = median_cor_number_edges))

pdf("Figure_3e.pdf", heigh=6)
ggplot(median_cor_number_edges, aes(x=Correlation_normal, y=Correlation_tumour))+
  xlim(-1,1) + ylim(-1,1)+
  geom_point(aes(color=Interaction_type))+
  ###stat_density2d(aes(color=Interaction_type, alpha= ..level..), geom="density2d",
  ###               size=2, contour=TRUE)+
  geom_abline(intercept=0, slope=1, col="black", lty="dashed", size=.25)+
  ###geom_abline(intercept=coeff[1], slope=coeff[2], col="red", lty="dashed", size=.25)+
  stat_smooth(data=median_cor_number_edges, method="lm", se=FALSE, colour="blue", size=.5)+
  #stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "UC_MC"), method="lm", se=FALSE, colour="darkgreen")+
  #stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "UC_UC"), method="lm", se=FALSE, colour="red")+
  #stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "MC_MC"), method="lm", se=FALSE, colour="blue")+
  geom_vline(xintercept = 0, col="black", size=.25)+
  geom_hline(yintercept=0, col="black", size=.25)+
  #geom_text_repel(aes(label = paste(Pathway1, Pathway2, sep="\n")), size=3) +
  theme_bw()+
  xlab("Correlation in normals")+
  ylab("Correlation in tumours")+
  theme(text = element_text(size=18),axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.title=element_blank(),
        legend.position=c(.85,.15))
dev.off()


median_cor_number_edges$Interaction_type <- factor(median_cor_number_edges$Interaction_type,
                                                   levels=c("UC_UC", "UC_MC", "MC_MC"))
pdf("Supp_17a.pdf", heigh=5)
ggplot(median_cor_number_edges, aes(x=Correlation_normal, y=Correlation_tumour))+
  xlim(-1,1) + ylim(-1,1)+
  geom_point(aes(color=Interaction_type), size=1)+
  ###stat_density2d(aes(color=Interaction_type, alpha= ..level..), geom="density2d",
  ###               size=2, contour=TRUE)+
  geom_abline(intercept=0, slope=1, col="black", lty="dashed", size=.25)+
  ###geom_abline(intercept=coeff[1], slope=coeff[2], col="red", lty="dashed", size=.25)+
  #stat_smooth(data=median_cor_number_edges, method="lm", se=FALSE, colour="blue", size=.5)+
  stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "UC_MC"), method="lm", se=FALSE, colour="darkgreen")+
  stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "UC_UC"), method="lm", se=FALSE, colour="red")+
  stat_smooth(data=subset(median_cor_number_edges, Interaction_type == "MC_MC"), method="lm", se=FALSE, colour="blue")+
  geom_vline(xintercept = 0, col="black", size=.25)+
  geom_hline(yintercept=0, col="black", size=.25)+
  #geom_text_repel(aes(label = paste(Pathway1, Pathway2, sep="\n")), size=3) +
  theme_bw()
dev.off()

pos_to_neg <- median_cor_number_edges[intersect(which(median_cor_number_edges$Correlation_normal>0),
                                  which(median_cor_number_edges$Correlation_tumour<0)),]
neg_to_pos <- median_cor_number_edges[intersect(which(median_cor_number_edges$Correlation_normal<0),
                                                which(median_cor_number_edges$Correlation_tumour>0)),]

write.table(median_cor_number_edges, file="Table_of_correlation.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)


edges_flip <- rbind(pos_to_neg, neg_to_pos)

write.table(edges_flip, file="Table_of_correlation_flip.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

#Percentage of edge types with positive and negative correlations
cor_number_edges_df <- vector()
for(tumour in tumours){
  cor_number_edges_df <- rbind(cor_number_edges_df,
                               data.frame(cor_number_edges[[tumour]], tumour_type=tumour))
}

cor_number_edges_df <- cor_number_edges_df[c(intersect_edge_type_with_cutoff(cor_number_edges_df, "UC_UC"),
                                             intersect_edge_type_with_cutoff(cor_number_edges_df, "UC_MC"),
                                             intersect_edge_type_with_cutoff(cor_number_edges_df, "MC_MC")),]
cor_number_edges_df$Interaction_type <- factor(cor_number_edges_df$Interaction_type, levels=c("UC_UC",
                                                                                              "UC_MC",
                                                                                              "MC_MC"))
cor_number_edges_df$Difference_correlation <- cor_number_edges_df$Correlation_tumour-cor_number_edges_df$Correlation_normal

edges_positive_negative <- cor_number_edges_df
edges_positive_negative <- cbind(cor_number_edges_df, 
                                           ifelse(edges_positive_negative$Correlation_normal>0, "Positive", "Negative"),
                                           ifelse(edges_positive_negative$Correlation_tumour>0, "Positive", "Negative"))
colnames(edges_positive_negative)[14:15] <- c("Direction_correlation_normal", "Direction_correlation_tumour")

edges_positive_negative$Interaction_type <- factor(edges_positive_negative$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))

#Barplots for all tumours
ag_cor_tumour <- aggregate(Direction_correlation_tumour ~ Interaction_type+tumour_type, edges_positive_negative, table)
ag_cor_normal <- aggregate(Direction_correlation_normal ~ Interaction_type+tumour_type, edges_positive_negative, table)

range_cor_tumour <- aggregate(Direction_correlation_tumour ~ Interaction_type, ag_cor_tumour, range)
range_cor_normal <- aggregate(Direction_correlation_normal ~ Interaction_type, ag_cor_normal, range)

median_cor_tumour <- aggregate(Direction_correlation_tumour ~ Interaction_type, ag_cor_tumour, median)
median_cor_normal <- aggregate(Direction_correlation_normal ~ Interaction_type, ag_cor_normal, median)

min_of_negative_counts_t <- range_cor_tumour[,2][,1]+median_cor_tumour[,3]
min_of_negative_counts_n <- range_cor_normal[,2][,1]+median_cor_normal[,3]
max_of_negative_counts_t <- range_cor_tumour[,2][,2]+median_cor_tumour[,3]
max_of_negative_counts_n <- range_cor_normal[,2][,2]+median_cor_normal[,3]

edges_positive_negative_all <- data.frame(Interaction_type=rep(c("UC_UC", "UC_MC", "MC_MC"),4),
                                          Median_count=c(median_cor_tumour[,3], median_cor_tumour[,2],
                                                         median_cor_normal[,3], median_cor_normal[,2]),
                                          Min_count = c(range_cor_tumour[,3][,1],min_of_negative_counts_t,
                                                        range_cor_normal[,3][,1],min_of_negative_counts_n),
                                          Max_count = c(range_cor_tumour[,3][,2],max_of_negative_counts_t,
                                                         range_cor_normal[,3][,2],max_of_negative_counts_n),
                                          Direction = rep(c(rep("Positive",3), rep("Negative",3)),2),
                                          Tissue_type = c(rep("Tumour",6), rep("Normal",6)))


edges_positive_negative_all$Interaction_type <- factor(edges_positive_negative_all$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))

t_n <- melt(cbind(ag_cor_normal[,1:2],unlist(ag_cor_normal[,3])))
t_t <- melt(cbind(ag_cor_tumour[,1:2],unlist(ag_cor_tumour[,3])))
ag_cor <- rbind(cbind(t_n, name="Normal"), cbind(t_t, name="Tumour"))
colnames(ag_cor)[3:5] <- c("Direction", "Number", "Tissue_type") 

ag_cor$Median <- edges_positive_negative_all$Median_count[match(interaction(ag_cor$Interaction_type, ag_cor$Direction, ag_cor$Tissue_type), 
                                                                interaction(edges_positive_negative_all$Interaction_type,
                                                                            edges_positive_negative_all$Direction, edges_positive_negative_all$Tissue_type))]
ag_cor$Number[which(ag_cor$Direction=="Negative")] <- ag_cor$Number[which(ag_cor$Direction=="Negative")]+
                                                      ag_cor$Median[which(ag_cor$Direction=="Positive")]

pdf("Figure_3d.pdf", width=5)
ggplot() +
  geom_bar(data=edges_positive_negative_all, aes(x=Interaction_type, y=Median_count, fill=Direction), stat="identity")+
  scale_fill_manual(values = c("orange","blue")) +
  geom_point(data=ag_cor, aes(x=Interaction_type, y=Number, fill=Direction), colour="black", shape=21, position=position_dodge(.2))+
  #geom_errorbar(aes(ymin=Min_count, ymax=Max_count), width=.1, position=position_dodge(.1))+
  facet_grid(.~Tissue_type)+
  ylab("Median count")+
  xlab("")+
  theme_bw()+
  theme(strip.background = element_rect(colour="black", fill="white"),
        text = element_text(size=18),axis.text.x= element_text(size=15, angle=45, hjust=1),
        axis.text.y= element_text(size=20),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        strip.text.x = element_text(size=15, face="bold"))
dev.off()




#Enrichment in positive and negative edges
for(tumour in tumours){
  temp <- edges_positive_negative[which(edges_positive_negative$tumour_type == tumour),]
  normal_positive_global <- sum(temp$Direction_correlation_normal == "Positive")
  normal_negative_global <- sum(temp$Direction_correlation_normal == "Negative")
  tumour_positive_global <- sum(temp$Direction_correlation_tumour == "Positive")
  tumour_negative_global <- sum(temp$Direction_correlation_tumour == "Negative")
  
  print(tumour)
  for(type in c("UC_UC", "UC_MC", "MC_MC")){
    print(type)
    temp_local <- temp[which(temp$Interaction_type == type),]
    normal_positive_local <- sum(temp_local$Direction_correlation_normal == "Positive")
    normal_negative_local <- sum(temp_local$Direction_correlation_normal == "Negative")
    
    tumour_positive_local <- sum(temp_local$Direction_correlation_tumour == "Positive")
    tumour_negative_local <- sum(temp_local$Direction_correlation_tumour == "Negative")
    
    print(paste("Enrichment normal positive:", round(phyper(q=normal_positive_local, m=normal_positive_global, n=normal_negative_global, k=normal_positive_local+normal_negative_local, lower.tail=FALSE),4)))
    print(paste("Enrichment normal negative:", round(phyper(q=normal_negative_local, m=normal_negative_global, n=normal_positive_global, k=normal_positive_local+normal_negative_local, lower.tail=FALSE),4)))
    
    print(paste("Enrichment tumour positive:", round(phyper(q=tumour_positive_local, m=tumour_positive_global, n=tumour_negative_global, k=tumour_positive_local+tumour_negative_local, lower.tail=FALSE),4)))
    print(paste("Enrichment tumour negative:", round(phyper(q=tumour_negative_local, m=tumour_negative_global, n=tumour_positive_global, k=tumour_positive_local+tumour_negative_local, lower.tail=FALSE),4)))
    
  }
}

#Fisher test to see whether positive or negative edges are enriched in tumours in each edge type
for(tumour in tumours){
  local <- edges_positive_negative[which(edges_positive_negative$tumour_type == tumour),]
  print(tumour)
  for(edge_type in c("UC_UC", "UC_MC", "MC_MC")){
    print(edge_type)
    temp <- local[which(local$Interaction_type == edge_type),]
    normal_positive <- sum(temp$Direction_correlation_normal == "Positive")
    normal_negative <- sum(temp$Direction_correlation_normal == "Negative")
    tumour_positive <- sum(temp$Direction_correlation_tumour == "Positive")
    tumour_negative <- sum(temp$Direction_correlation_tumour == "Negative")
  
   
    a <- fisher.test(cbind(c(tumour_positive, normal_positive), c(tumour_negative, normal_negative)), alternative="greater")$p.val
    b <- fisher.test(cbind(c(tumour_negative, normal_negative),c(tumour_positive, normal_positive)), alternative="greater")$p.val
    
    print(paste("Greater positive in tumour", a))
    print(paste("Greater negative in tumour", b))
  }
}
  
####Variance analysis

genes_slim <- #List with entrez IDs as names, GOslims as items

#Differences in the variances of the edges across normals and tumours
var_cor_edges_n <- aggregate(Correlation_normal ~ Pathway1+Pathway2+Degree_of_interaction+Interaction_type, cor_number_edges_df, var)
var_cor_edges_t <- aggregate(Correlation_tumour ~ Pathway1+Pathway2+Degree_of_interaction+Interaction_type, cor_number_edges_df, var)

var_cor_edges <- data.frame(cbind(var_cor_edges_n, var_cor_edges_t[,5]))
colnames(var_cor_edges)[5:6] <- c("Var_normal", "Var_tumour")

var.test(var_cor_edges$Var_tumour, var_cor_edges$Var_normal, alternative="less") #Variance of the correlation of the edges in tumours is less than in normals
wilcox.test(var_cor_edges$Var_tumour, var_cor_edges$Var_normal, alternative="less")
var_cor_edges$Difference_var <- var_cor_edges$Var_tumour - var_cor_edges$Var_normal

var_cor_edges_m <- as.matrix(var_cor_edges)
var_cor_edges_gg <- data.frame(rbind(cbind(var_cor_edges_m[,1:5], "Normal"),
                                     cbind(var_cor_edges_m[,c(1:4,6)], "Tumour")))
colnames(var_cor_edges_gg) <- c(colnames(var_cor_edges)[1:4], "Variance", "Tissue_type")                              

var_cor_edges_gg$Variance <- as.numeric(as.character(var_cor_edges_gg$Variance))
var_cor_edges_gg$Interaction_type <- factor(var_cor_edges_gg$Interaction_type, levels=c("UC_UC", "UC_MC", "MC_MC"))

pdf("Figure_3f.pdf", height=5)
ggplot(var_cor_edges_gg, aes(x=Tissue_type, y=Variance, fill=Tissue_type)) + 
  geom_boxplot()+
  theme_bw()+
  geom_point(position=position_jitter(.2))+
  ylab("Variance of edges across tissues")+
  facet_grid(.~Interaction_type)+
  guides(fill=FALSE)+
  theme(axis.title.x=element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        text = element_text(size=18),axis.text.x= element_text(size=20, angle=45, hjust=1),
        axis.text.y= element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size=15, face="bold"))
dev.off()  

cor_number_edges_df_cor <- rbind(data.frame(cor_number_edges_df[c(1:6)], Correlation=cor_number_edges_df[,7], Tissue_type="Normal", cor_number_edges_df[c(9:13)]),
                                 data.frame(cor_number_edges_df[c(1:6)], Correlation=cor_number_edges_df[,8], Tissue_type="Tumour", cor_number_edges_df[c(9:13)]))
maximum <- aggregate(Correlation~Pathway1+Pathway2, cor_number_edges_df_cor, max)
maximum <- maximum[order(maximum$Correlation, decreasing=TRUE),]

cor_number_edges_df_cor$Pathway_interaction <- factor(interaction(cor_number_edges_df_cor$Pathway1, cor_number_edges_df_cor$Pathway2),
                                                                                          levels=interaction(maximum[,1], maximum[,2]))

pdf("Supp_18.pdf", height=5, width=10)
ggplot(cor_number_edges_df_cor, aes(x=Pathway_interaction, y=Correlation)) + 
  theme_bw()+
  geom_point(y =1.01, aes(colour=Interaction_type))+
  geom_path(data=cor_number_edges_df_cor[cor_number_edges_df_cor$Tissue_type=="Normal",], aes(group=interaction(interaction(Pathway1, Pathway2), Tissue_type)),
                colour="lightblue")+
  geom_point(data=cor_number_edges_df_cor[cor_number_edges_df_cor$Tissue_type=="Tumour",], color="black", size=.8)+
  geom_hline(yintercept=0, size=.1)+
  xlab("Cellular processes")+
  #facet_grid(.~Interaction_type)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()



#Searching for another pair of ME processes (UC-MC)

edges_flip_UC_MC <- edges_flip[edges_flip$Interaction_type == "UC_MC",]
edges_flip_UC_MC <- edges_flip_UC_MC[edges_flip_UC_MC$Correlation_normal > 0,]

edges_flip_UC_MC

local_pos_to_neg <- list()
for(i in 1:nrow(edges_flip_UC_MC)){
  pair <- edges_flip_UC_MC[i, 1:2]
  pair_name <- paste(as.character(unlist(pair[1,])), collapse="-")
  local_pos_to_neg[[pair_name]] <- vector()
  for(tumour in tumours){
    cor_normal <- cor_number_edges[[tumour]][intersect(which(cor_number_edges[[tumour]]$Pathway1 == as.character(pair[1,1])),
                                                       which(cor_number_edges[[tumour]]$Pathway2 == as.character(pair[1,2]))), "Correlation_normal"]
    
    cor_tumour <- cor_number_edges[[tumour]][intersect(which(cor_number_edges[[tumour]]$Pathway1 == as.character(pair[1,1])),
                                                       which(cor_number_edges[[tumour]]$Pathway2 == as.character(pair[1,2]))), "Correlation_tumour"]
    local_pos_to_neg[[pair_name]] <- c(local_pos_to_neg[[pair_name]], cor_normal > 0 && cor_tumour < 0)
  }
  local_pos_to_neg[[pair_name]] <- sum(local_pos_to_neg[[pair_name]])
}
