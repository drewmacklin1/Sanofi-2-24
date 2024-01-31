setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
library(DESeq2)
library(tidyverse)

studylist<-c("AD_EC_SFG_PMID33432193", "ALS_Spinalcord_Internal_Zelic", "ALS_FTLD_PrimMotorCortex_GSE174332","MS_Cortex_Internal_Proto", "MS_CorticalSubcorticalLesions_PRJNA544731", "MS_Lesions_PMID34497421", "PD_Midbrain_PMID34919646", "PD_PreFCortex_PPR453641", "FTD_FrontTempCortex_PMID35879464", "PD_Putamen_Internal_Unknown", "PD_SubNigra_PMID35513515", "AD_OC_OTC_PMID33609158")
pseudo_locations <- NULL
for(i in 1:length(studylist)){
  files <- list.files(path=paste0('datasets/', studylist[i], '/pseudobulk_counts_matrices'))
  for(j in 1:length(files)){
    if(grepl("icrogli", files[j])){
      pseudo_locations<-c(pseudo_locations, paste0('datasets/', studylist[i], '/pseudobulk_counts_matrices/', files[j]))
    }
  }
}
pseudo_locations<-unique(pseudo_locations)

odds <- function(x) x[ x %% 2 == 1 ]
modifier <- odds(1:length(pseudo_locations))
for (j in modifier){
  
  #read in files, and convert to dataframes
  data_counts <- as.data.frame(read.csv(pseudo_locations[j]))
  rownames(data_counts)<-data_counts$X
  data_counts<-data_counts[-1]
  data_counts <- round(data_counts)
  metadata <- as.data.frame(read.csv(pseudo_locations[j+1]))
  rownames(metadata)<-metadata$X
  metadata<-metadata[-1]
  rownames(metadata) <- colnames(data_counts)
  metadata$sample <- rownames(metadata)
  metadata$condition <- as.factor(metadata$condition)
  
  #run DESeq2
  if(length(unique(metadata$gender))>1){
    dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition + gender)
  }else{
    dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
  }
  dds <- estimateSizeFactors(dds)
  #keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
  #dds <- dds[keep,]
  dds <- DESeq(dds)
  
  celltype<- gsub("_pseudobulk_data.csv", "", gsub("datasets/AD_EC_SFG_PMID33432193/pseudobulk_counts_matrices/", "", pseudo_locations[j]))
  
  counts <- as.data.frame(counts(dds, normalized = TRUE))
  control<-rownames(metadata[metadata$condition=='control',])
  disease<-rownames(metadata[metadata$condition!='control',])
  micro_expression<- data.frame(matrix(ncol = 2, nrow = length(rownames(counts))))
  colnames(micro_expression) <- c("disease_exp", "control_exp")
  rownames(micro_expression) <- rownames(counts)
  
  for(a in 1:length(rownames(counts))){
    micro_expression[rownames(counts)[a],"control_exp"] <- sum(counts[a,control])/length(control)
    micro_expression[rownames(counts)[a],"disease_exp"] <- sum(counts[a,disease])/length(disease)
  }
  micro_expression$dis_over_ctrl <- micro_expression$disease_exp/micro_expression$control_exp
  write.csv(micro_expression, paste0("microglia_markers/", gsub("datasets/", "", gsub("/pseudobulk_counts_matrices/" , "", gsub("pseudobulk_data", "", pseudo_locations[j])))))
}
  
files<-list.files("microglia_markers")
a <- read.csv(paste0("microglia_markers/", files[1]))
b <- read.csv(paste0("microglia_markers/", files[2]))
c <- read.csv(paste0("microglia_markers/", files[3]))
d <- read.csv(paste0("microglia_markers/", files[4]))
e <- read.csv(paste0("microglia_markers/", files[6]))
f <- read.csv(paste0("microglia_markers/", files[7]))
g <- read.csv(paste0("microglia_markers/", files[8]))
h <- read.csv(paste0("microglia_markers/", files[9]))
i <- read.csv(paste0("microglia_markers/", files[10]))
j <- read.csv(paste0("microglia_markers/", files[11]))
k <- read.csv(paste0("microglia_markers/", files[12]))

merged<-left_join(a, b, by = 'X')
merged<-left_join(merged, c, by = 'X')
merged<-left_join(merged, d, by = 'X')
merged<-left_join(merged, e, by = 'X')
merged<-left_join(merged, f, by = 'X')
merged<-left_join(merged, g, by = 'X')
merged<-left_join(merged, h, by = 'X')
merged<-left_join(merged, i, by = 'X')
merged<-left_join(merged, j, by = 'X')
merged<-left_join(merged, k, by = 'X')

dis_over_ctrl_list<-colnames(merged)[grepl("dis_over_ctrl", colnames(merged))]
ctrl_list<-colnames(merged)[grepl("control", colnames(merged))]
dis_list<-colnames(merged)[grepl("disease", colnames(merged))]
merged$mean_control<-rowMeans(merged[,ctrl_list], na.rm = TRUE)
merged$mean_disease<-rowMeans(merged[,dis_list], na.rm = TRUE)
merged$mean_ratio<-merged$mean_disease/merged$mean_control
merged<-merged[,c("X", "mean_ratio", "mean_disease", "mean_control", dis_list, ctrl_list)]

write.csv(merged, 'microglia_markers/results.csv')

View(merged[merged$mean_disease>40 & merged$mean_control<11,]) #PNLDC1, CRYM-AS1
View(merged[merged$mean_disease>140 & merged$mean_control<35,]) #HSPA4L, HSPE1-MOB4
View(merged[merged$mean_disease>300 & merged$mean_control<150,]) #NUPR1, FAM46A, HSPA6, MB21D1