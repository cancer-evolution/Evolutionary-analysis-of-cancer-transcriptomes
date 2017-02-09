library(qusage)
library(mGSZ)
library(limma)
library(edgeR)
library(GSEABase)
library(GO.db)
library(reshape2)
library(ggplot2)
library(scales)
library(GSVA)
library(grid)
library(gplots)
library(ggdendro)
library(gridExtra)
library(ggrepel)
library(igraph)

source("Activity_levels_processes_and_age_functions.R")

tumours <- c("LUAD", "LUSC", "BRCA","PRAD","LIHC", "COAD","STAD")

for (i in tumours)
{
  qusage_script(i, "goslims")
}

qusage_results <- list()
qusage_results_matrix <- vector()
for (type in tumours){
  qusage_results[[type]] <- read.delim(paste(type,"_goslims_selfcontained_qusage_results.txt", sep=""))  #Results from qusage. Columns: pathway.name, log.fold.change, p.Value, FDR
  qusage_results_names <- qusage_results[[type]][,1]
  qusage_results[[type]] <- t(apply(qusage_results[[type]],1,function(local){
    if (is.na(local[4])){
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

GO_ages <- ##Data frame with GOterm (ID), GOname, and Age
colnames(GO_ages) <- c("GOterm", "GOname", "Age")


#Continuous tile plot
qusage_results_melt <- melt(qusage_results_matrix)
colnames(qusage_results_melt) <- c("Tumour", "GOSlim", "FC")
qusage_results_melt <- data.frame(qusage_results_melt, Age=GO_ages[match(qusage_results_melt$GOSlim, GO_ages[,2]), 3])

#ordered ggplot
t_N <- c("plasma membrane organization","cytoskeleton organization","aging","embryo development",
"vesicle-mediated transport","secondary metabolic process","transposition","growth",
"homeostatic process","pigmentation","cell division","protein maturation")

t_M <- c("reproduction","cell morphogenesis","immune system process","circulatory system process",
         "cell adhesion","signal transduction","cell-cell signaling","cell death","cell proliferation",
         "developmental maturation","cell differentiation","extracellular matrix organization",
         "cell junction organization","locomotion","anatomical structure formation involved in morphogenesis",
         "anatomical structure development","cell motility","neurological system process",
         "response to stress - multicellular")

t_U <- c("transmembrane transport","sulfur compound metabolic process","lipid metabolic process",
         "transport","protein complex assembly","vacuolar transport","cellular component assembly",
         "small molecule metabolic process","carbohydrate metabolic process","cofactor metabolic process",
         "autophagy","catabolic process",
         "generation of precursor metabolites and energy",
         "membrane organization",
         "cellular protein modification process",
         "macromolecular complex assembly",
         "cellular amino acid metabolic process",
         "cytoskeleton-dependent intracellular transport","biosynthetic process","response to stress - general",
         "cellular nitrogen compound metabolic process","nitrogen cycle metabolic process","chromosome segregation",
         "translation","protein targeting","nucleobase-containing compound catabolic process","mitotic nuclear division",
         "cell cycle","chromosome organization","mitochondrion organization","protein folding",
         "nucleocytoplasmic transport","symbiosis, encompassing mutualism through parasitism","DNA metabolic process",
         "mRNA processing","ribonucleoprotein complex assembly","tRNA metabolic process","ribosome biogenesis")

qusage_results_melt_ordered <- qusage_results_melt
qusage_results_melt_ordered$GOSlim <- factor(qusage_results_melt_ordered$GOSlim, levels=c(t_N, t_M, t_U))
qusage_results_melt_ordered$Age <- factor(qusage_results_melt_ordered$Age, levels=c("UC", "MC"))

qusage_results_melt_ordered <- qusage_results_melt_ordered[!is.na(qusage_results_melt_ordered$Age),]
g_heat <- ggplot(qusage_results_melt_ordered, aes(x=Tumour, y =GOSlim))+
  geom_tile(aes(fill=FC)) +
  scale_fill_gradientn(colours=c("darkblue", "blue", "royalblue1",
                                 "white",
                                 "red", "red", "darkred"),
                       values=rescale(c(-1.5, -.7, -0.3,
                                        0,
                                        0.3, 0.5)),
                       guide="colorbar", name="logFC")+
  scale_y_discrete(expand = c(0, 0.25))+
  #geom_point(data=qusage_results_melt_ordered, aes(x=7.55, y=1:10, color=Age), size=3)+
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.text=element_text(size=17), legend.title=element_text(size=17),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_text(colour="black", size=20),
        #axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=16, colour="black"),
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.y = element_text(size=18, face="bold"))+
  facet_grid(Age~., shrink = FALSE, scales='free', space='free')

pdf("Figure_2a1.pdf", width=9.5, height=15)
g_heat
dev.off()


#Difference between UC and MC components
gene_ages <- #Data framew with GeneID, Entrez, and Phylostrata
slim_terms <- getOBOCollection("goslim_generic.obo") #downloaded from GO website
go_bp_children <- GOBPOFFSPRING$"GO:0008150"

gene_ages[,3] <- as.numeric(gene_ages[,3])

slim_terms <- ids(slim_terms)[which(ids(slim_terms) %in% go_bp_children)]

gonames <- Term(GOTERM)

slim_genes <- #list of GOIDs as names, entrez as items

slim_genes_ages <- lapply(names(slim_genes), function(slim){
  local_genes <- slim_genes[[slim]]
  local_genes_ages <- cbind(local_genes, gene_ages[match(local_genes, gene_ages[,2]),3])
  UC_genes <- local_genes_ages[which(local_genes_ages[,2] %in% 1:3),1]
  MC_genes <- local_genes_ages[which(local_genes_ages[,2] %in% 4:16),1]
  return(list(UC=UC_genes, MC=MC_genes))
})
names(slim_genes_ages) <- names(slim_genes)

qusage_slim_UC_MC_genes <- list()
qusage_slim_UC_MC_genes_object <- list()
tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")
for (type in tumours){
  load(paste(type, "_voom_expression.RData", sep="")) #objects created after voom normalization in limma pipeline. See limma guide
  load(paste(type, "_design_matrix.RData", sep="")) #design matrix after voom normalization in limma pipeline. See limma guide

  qusage_slim_UC_MC_genes[[type]] <- sapply(names(slim_genes_ages), function(GOslim){
    print(GOslim)
    return(list(qusage_script_explicit(slim_genes_ages[[GOslim]], as.matrix(expression_v$E), expression_design)))
  })
  print(paste("Done with", type))
}
save(qusage_slim_UC_MC_genes, file="qusage_slim_UC_MC_genes.Rdata")
save(qusage_slim_UC_MC_genes_object, file="qusage_slim_UC_MC_genes_object.Rdata")

load("qusage_slim_UC_MC_genes_stress_separated.Rdata")

qusage_slim_UC_MC_genes_df <- vector()
for(type in tumours){
  for(slim in names(qusage_slim_UC_MC_genes[[type]])){
    local_results <- qusage_slim_UC_MC_genes[[type]][[slim]]
    if(!is.null(local_results)){
      for(i in 1:nrow(local_results)){
        if(local_results[i,4]<0.05){
          qusage_slim_UC_MC_genes_df <- rbind(qusage_slim_UC_MC_genes_df,
                                              c(slim,paste(as.character(local_results[i,1]), "genes"),
                                                type,local_results[i,2], as.character(local_results[i,1])))      
        }
      }
    }
  }
}

colnames(qusage_slim_UC_MC_genes_df) <- c("GOslims", "GOslim_children", "Tumour", "FC", "Age")

qusage_slim_UC_MC_genes_df <- as.data.frame(qusage_slim_UC_MC_genes_df) 
qusage_slim_UC_MC_genes_df[,4] <- as.numeric(as.character(qusage_slim_UC_MC_genes_df[,4]))

qusage_slim_UC_MC_genes_df$Tumour <- factor(qusage_slim_UC_MC_genes_df$Tumour, levels=tumours)

median_logFC <- aggregate(FC~GOslims+Age,qusage_slim_UC_MC_genes_df,median)
min_log_FC <- aggregate(FC~GOslims+Age,qusage_slim_UC_MC_genes_df,min)
max_log_FC <- aggregate(FC~GOslims+Age,qusage_slim_UC_MC_genes_df,max)


#Difference in the logFC of UC and MC components
all_slims <- as.character(unique(qusage_slim_UC_MC_genes_df[,1]))
all_slims_p_val <- vector()
for(slim in all_slims){
  temp <- qusage_slim_UC_MC_genes_df[which(qusage_slim_UC_MC_genes_df[,1] == slim),]
  temp_UC <- temp[temp[,2]=="UC genes",]
  temp_MC <- temp[temp[,2]=="MC genes",]
  if(length(temp_UC[,4]) != 0 && length(temp_MC[,4]) != 0){
    p <- wilcox.test(temp_UC[,4], temp_MC[,4])$p.value
    p_up <- wilcox.test(temp_UC[,4], temp_MC[,4], alternative="greater")$p.value
    p_down <- wilcox.test(temp_UC[,4], temp_MC[,4], alternative="less")$p.value
    all_slims_p_val <- rbind(all_slims_p_val,
                             c(slim, p, p_up, p_down))
  } 
}
all_slims_p_val <- as.data.frame(all_slims_p_val)
colnames(all_slims_p_val) <- c("slim", "p.ts", "p.UCup", "p.UCdown")
all_slims_p_val$p.ts <- as.numeric(as.character(all_slims_p_val$p.ts))
all_slims_p_val$p.UCup <- as.numeric(as.character(all_slims_p_val$p.UCup))
all_slims_p_val$p.UCdown <- as.numeric(as.character(all_slims_p_val$p.UCdown))

all_slims_p_val$p.ts.adj <- p.adjust(all_slims_p_val[,"p.ts"], method="BH")
all_slims_p_val$p.up.adj <- p.adjust(all_slims_p_val[,"p.UCup"], method="BH")
all_slims_p_val$p.down.adj <- p.adjust(all_slims_p_val[,"p.UCdown"], method="BH")


median_logFC <- data.frame(median_logFC, Min=min_log_FC$FC, Max=max_log_FC$FC)


median_logFC$SlimAges <- GO_ages[match(median_logFC$GOslims, GO_ages[,1]),"Age"]


#Order by heatmap
median_logFC$GOslims <- factor(median_logFC$GOslims, levels=c(t_U, t_M, t_N))

median_logFC$SlimAges <- factor(median_logFC$SlimAges, levels=c("UC", "MC"))
median_logFC$Age <- factor(median_logFC$Age, levels=c("UC", "MC"))

median_logFC$p.adj <- all_slims_p_val[match(median_logFC$GOslims, all_slims_p_val$slim), "p.ts.adj"]
median_logFC$Significance <- ifelse(median_logFC$p.adj<0.05, 4, NA)

median_logFC$p.adj.up <- all_slims_p_val[match(median_logFC$GOslims, all_slims_p_val$slim), "p.up.adj"]
median_logFC$Significance_up <- ifelse(median_logFC$p.adj.up<0.05, 100, NA)

median_logFC$p.adj.down <- all_slims_p_val[match(median_logFC$GOslims, all_slims_p_val$slim), "p.down.adj"]
median_logFC$Significance_down <- ifelse(median_logFC$p.adj.down<0.05, 4, NA)


median_logFC$Significance <- as.numeric(as.character(median_logFC$Significance))

median_logFC <- median_logFC[!is.na(median_logFC$SlimAges),]

#Calculate percentage of processes where UC is up
temp_median_logFC <- median_logFC[median_logFC$Age == "UC",]  #select those with a UC component
temp_median_logFC <- temp_median_logFC[!is.na(temp_median_logFC$p.adj),] #exclude those with a p value = NA. These do not have a MC component
sum(!is.na(temp_median_logFC$Significance_up)) 
length(temp_median_logFC$Significance_up)

g_diff_B <- ggplot(median_logFC, aes(y=GOslims, x=FC))+
  geom_point(aes(colour=Age), size=3)+
  geom_hline(yintercept=0, size=.2)+
  geom_errorbarh(aes(xmax = Max, xmin = Min, colour=Age), height=.05)+
  scale_colour_manual(values=c("red", "blue"))+
  theme_bw()+
  theme(legend.position="bottom")+
  geom_point(aes(size=Significance_up), colour="black", x=-1.9, shape=4)+
  #geom_point(aes(size=Significance), x=-2.2)+
  scale_size(range=c(2,0), guide=FALSE)+
  geom_vline(xintercept=0, col="black", size=.25)+
  xlab("logFC")+
  theme(axis.text.y = element_text(size=20),
        text = element_text(size=18),axis.text.x= element_text(size=18),
        axis.title.x = element_text(size=20), axis.title.y = element_blank())+
  facet_grid(SlimAges~., shrink = FALSE, scales='free', space='free')


pdf("Figure_2a2.pdf", width=10, height=16)
g_diff_B
dev.off()

#Difference in the GOslim components
qusage_slim_UC_MC_genes_df2 <- data.frame(qusage_slim_UC_MC_genes_df, apply(qusage_slim_UC_MC_genes_df,1,function(row){
  slim <- as.character(row[1])
  age <- as.character(row[5])
  return(c(length(slim_genes_ages[[slim]][[age]])))
}))
colnames(qusage_slim_UC_MC_genes_df2)[6] <- "Number_of_genes"

qusage_slim_UC_MC_genes_diff_df <- calculate_difference_between_UC_MC_components(qusage_slim_UC_MC_genes_df2, TRUE)
error_bars_genes <- create_df_of_error_bars(qusage_slim_UC_MC_genes_diff_df)

error_bars_genes$Direction <- ifelse(error_bars_genes$Median_FC>0, "Upregulated", "Downregulated")
error_bars_genes$Direction <- factor(error_bars_genes$Direction, levels=c("Upregulated", "Downregulated"))
error_bars_genes$Driven <- ifelse(error_bars_genes$Median_logFC_difference>0, "Driven by UC", "Driven by MC")
error_bars_genes$Driven <- factor(error_bars_genes$Driven, levels=c("Driven by UC", "Driven by MC"))


g_abs_diff <- plot_error_FC_diff_error_bars(error_bars_genes)

pdf("Figure_2b.pdf", width=8, height=5.3)
g_abs_diff
dev.off()
