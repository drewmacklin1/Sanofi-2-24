#compares the similarity of all differential expression files in the database and creates a heatmap
library(tidyverse)
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")
studylist<-c("AD_EC_SFG_PMID33432193", "ALS_Spinalcord_Internal_Zelic", "ALS_FTLD_PrimMotorCortex_GSE174332", "FTD_FrontTempCortex_PMID35879464", "MS_Cortex_Internal_Proto", "MS_CorticalSubcorticalLesions_PRJNA544731", "MS_Lesions_PMID34497421", "PD_Midbrain_PMID34919646", "PD_PreFCortex_PPR453641", "PD_Putamen_Internal_Unknown", "PD_SubNigra_PMID35513515", "AD_OC_OTC_PMID33609158")
all_genes <- read.csv("results_target_credentialing/all_genes.csv", row.names = 1)
all_genes<-all_genes$x

#find all diffex paths
all_dataset_locations <-NULL
for(a in 1:length(studylist)){
  files<-list.files(paste0("datasets/", studylist[a], "/", "differential_expression_results"))
  for(b in 1:length(files)){
      all_dataset_locations<-c(all_dataset_locations, paste0("datasets/", studylist[a], "/", "differential_expression_results/", files[b]))
  }
}

all_cells_heatmap<- data.frame(matrix(ncol = length(all_dataset_locations), nrow = length(all_genes)))
colnames(all_cells_heatmap) <- gsub("vs.Control.csv", "", gsub("__vs.__Control.csv", "", gsub("/differential_expression_results/diffex_results_", "_", gsub("datasets/", "", all_dataset_locations))))
rownames(all_cells_heatmap) <- all_genes

#create matrix with all the pvalues of all the genes from all the diffex files
for(a in 1:length(all_dataset_locations)){
  data<- read.csv(all_dataset_locations[a])
  for(b in 1:length(rownames(data))){
    all_cells_heatmap[data[b, "X"],colnames(all_cells_heatmap)[a]] <- data[b,"pvalue"]
  }
}

# keep<-rowSums(is.na(all_cells_heatmap))<2
# all_cells_heatmap<-all_cells_heatmap[keep,]
# keep<-rowSums(abs(all_cells_heatmap))>6
# all_cells_heatmap<-all_cells_heatmap[keep,]
# all_cells_heatmap<-as.matrix(all_cells_heatmap)
# all_cells_heatmap[is.na(all_cells_heatmap)]<-0
# plot<-heatmap(t(all_cells_heatmap))

#find the euclidean distance between each diffex file's pvalues in 63,000 dimensional space and create a heatmap
library(DESeq2)
library(pheatmap)
library("RColorBrewer")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDists <- dist(t(all_cells_heatmap))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(all_cells_heatmap)
colnames(sampleDistMatrix) <- NULL
ph<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, 
         col=colors,
         legend_breaks = c(0, round(max(sampleDistMatrix)/4), round(2*max(sampleDistMatrix)/4), round(3*max(sampleDistMatrix)/4), round(max(sampleDistMatrix)))-1, 
         main = "Raw P All Cells No Filter", 
         legend_labels = c(0, as.character(round(max(sampleDistMatrix)/4)), as.character(round(2*max(sampleDistMatrix)/4)), as.character(round(3*max(sampleDistMatrix)/4)), "N-Dimensional Distance"),
         legend = TRUE, 
         fontsize_row = 5, 
         cellheight = 5)
ggsave(filename = "RawP_All_Cells_No_Filter.png",
       plot = ph,
       device = 'png',
       path = "datasets",
       width = 20,
       height = 20,
       units = "in")