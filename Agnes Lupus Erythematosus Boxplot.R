setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
library(tidyverse)
library(Seurat)

pckd <- readRDS("Agnes/lupus.rds")
pckd_AU <- read.table(file = 'Agnes/lupus_431genes.tsv', sep = '\t', header = TRUE)

metadata <- as.data.frame(pckd@meta.data)
metadata$Barcodes<-rownames(metadata)
rownames<-rownames(metadata)
metadata<-left_join(metadata, pckd_AU, by = "Barcodes")
rownames(metadata)<-rownames

metadata$cell_type <- metadata$Author.s.cell.type..Level.1.
metadata$sample_id <- metadata$Sample.ID
metadata$condition <- as.character(metadata$Condition)
metadata$gender <- metadata$Gender
metadata$enrichment<-metadata$AUCell.score...AUCell.score..1.
metadata$condition[metadata$condition=='normal'] <- 'ctrl'
metadata$condition[metadata$condition=='systemic lupus erythematosus'] <- 'SLE'
metadata$cell_type[metadata$cell_type=='CD4-positive, alpha-beta T cell'] <- 'CD4T'
metadata$cell_type[metadata$cell_type=='classical monocyte'] <- 'Mono'
metadata$cell_type[metadata$cell_type=='B cell'] <- 'B'
metadata$cell_type[metadata$cell_type=='CD8-positive, alpha-beta T cell'] <- 'CD8T'
metadata$cell_type[metadata$cell_type=='natural killer cell'] <- 'NK'
metadata$cell_type[metadata$cell_type=='non-classical monocyte'] <- 'NC.Mono'
metadata$cell_type[metadata$cell_type=='prolif'] <- 'Prolif'
metadata$cell_type[metadata$cell_type=='conventional dendritic cell'] <- 'CDC'
metadata$cell_type[metadata$cell_type=='plasmacytoid dendritic cell'] <- 'PDC'
metadata$cell_type[metadata$cell_type=='plasmablast'] <- 'PB'
metadata$cell_type[metadata$cell_type=='progenitor cell'] <- 'Progen'
metadata$cell_type_condition<-paste(metadata$cell_type,metadata$condition)

print(unique(metadata$cell_type))
print(unique(metadata$condition))
print(unique(metadata$gender))
print(unique(metadata$sample_id))
print(unique(metadata$cell_type_condition))


#cell_condition_list<-sort(unique(metadata$cell_type_condition), decreasing = FALSE)
cell_condition_list<-c("NK ctrl", "NK SLE","CD4T ctrl", "CD4T SLE", "Mono ctrl", "Mono SLE", "NC.Mono ctrl", "NC.Mono SLE",  "PB ctrl", "PB SLE", "Progen ctrl", "Progen SLE","PDC ctrl", "PDC SLE", "CDC ctrl", "CDC SLE", "Prolif ctrl", "Prolif SLE", "CD8T ctrl","CD8T SLE","B ctrl", "B SLE")

data<- data.frame(matrix(ncol = length(cell_condition_list), nrow = 20000))
colnames(data)<-cell_condition_list

for(a in 1:length(cell_condition_list)){
  rownames<-rownames(subset(metadata, metadata$cell_type_condition==cell_condition_list[a]))
  cell_sample_enrichment <- metadata[rownames(metadata) %in% rownames, "enrichment"]
  data[1:length(cell_sample_enrichment),cell_condition_list[a]]<-cell_sample_enrichment
}

colors = rep(c("blue", "red"),14)
par(cex.axis=.75)
boxplot(data, las=2, boxfill = colors ,outline=FALSE, at = c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17, 19,20, 22,23, 25,26, 28,29, 31,32))
title(ylab="Enrichment",
      # xlab="Cell type", 
      main=paste("Systemic Lupus Erythematosus"))

ranking<- data.frame(matrix(ncol = length(cell_condition_list), nrow = 1))
colnames(ranking)<-cell_condition_list
means<-colMeans(data, na.rm = TRUE)
#means<-colSums(data, na.rm = TRUE)/(20000-(colSums(is.na(data))))
ranking[1,1]<-means[2]/means[1]
ranking[1,3]<-means[4]/means[3]
ranking[1,5]<-means[6]/means[5]
ranking[1,7]<-means[8]/means[7]
ranking[1,9]<-means[10]/means[9]
ranking[1,11]<-means[12]/means[11]
ranking[1,13]<-means[14]/means[13]
ranking[1,15]<-means[16]/means[15]
ranking[1,17]<-means[18]/means[17]
ranking[1,19]<-means[20]/means[19]
ranking[1,21]<-means[22]/means[21]
ranking<-as.data.frame(t(ranking))
View(ranking)
write.csv(ranking, "Agnes/Lupus_ratios.csv")










cell_condition_list<-sort(unique(metadata$cell_type_condition), decreasing = FALSE)
data<- data.frame(matrix(ncol = length(cell_condition_list), nrow = 20000))
colnames(data)<-cell_condition_list

for(a in 1:length(cell_condition_list)){
  rownames<-rownames(subset(metadata, metadata$cell_type_condition==cell_condition_list[a]))
  cell_sample_enrichment <- metadata[rownames(metadata) %in% rownames, "enrichment"]
  data[1:length(cell_sample_enrichment),cell_condition_list[a]]<-cell_sample_enrichment
}

par(cex.axis=.75)
boxplot(data, las=2,at = c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17, 19,20, 22,23, 25,26, 28,29, 31,32))
title(ylab="Enrichment",
      # xlab="Cell type", 
      main=paste("Systemic Lupus Erythematosus"))

significance<-data.frame(matrix(ncol = 22, nrow = 20))
odds <- function(x) x[ x %% 2 == 1 ]
modifier <- odds(1:length(colnames(data)))
for (j in modifier){
  q<-shapiro.test(na.omit(sample(data[,j], 5000, replace = FALSE, prob = NULL)))
  significance[2,j]<-paste(q[1], "W score")
  significance[3,j]<-paste(q[2], "P value")
  significance[1,j]<-paste(colnames(data)[j],q[3])

  z<-shapiro.test(na.omit(sample(data[,j+1], 5000, replace = FALSE, prob = NULL)))
  significance[6,j]<-paste(z[1], "W score")
  significance[7,j]<-paste(z[2], "P value")
  significance[5,j]<-paste(colnames(data)[j+1],z[3])

  #check variance of control vs disease var()
  significance[9,j]<-paste(var(data[,j], na.rm = TRUE)[1],colnames(data)[j], "Variance") 
  significance[10,j]<-paste(var(data[,j+1], na.rm = TRUE)[1], colnames(data)[j+1], "Variance")
  
  #two sample t-test t.test()
  t<-t.test(data[,j], data[,j+1], alternative = "two.sided", var.equal = FALSE)
  significance[12,j]<-paste(t[1], "statistic")
  significance[13,j]<-paste(t[2], "parameter")
  significance[14,j]<-paste(t[3], "T-test Pvalue")
  significance[15,j]<-paste(t[4], "Confidence interval")
  significance[16,j]<-paste(t[5], "Estimate")
  significance[17,j]<-paste(t[6], "Null value")
  significance[18,j]<-paste(t[7], "stderr")
  significance[19,j]<-paste(t[8], 'alternative')
  significance[20,j]<-paste(t[9], "method")
}
significance<-significance[,c(1,3,5,7,9,11,13,15,17,19,21)]


# condition_list<-c("ctrl","SLE")
# data2<-data.frame(matrix(ncol=2, nrow=5000))
# colnames(data2)<-condition_list
# for(a in 1:length(condition_list)){
#   rownames<-rownames(subset(metadata, metadata$condition==condition_list[a]))
#   cell_sample_enrichment <- metadata[rownames(metadata) %in% rownames, "enrichment"]
#   data2[1:length(cell_sample_enrichment),condition_list[a]]<-cell_sample_enrichment
# }
# #boxplot(data2, las=2)
# 
# #check normality of data shapiro.test()
# shapiro.test(na.omit(sample(data2$ctrl, 5000, replace = FALSE, prob = NULL)))
# shapiro.test(na.omit(sample(data2$SLE, 5000, replace = FALSE, prob = NULL)))
# 
# #check variance of control vs disease var()
# var(data2$ctrl, na.rm = TRUE)
# var(data2$SLE, na.rm = TRUE)
# 
# #two sample t-test t.test()
# t.test(data2$ctrl, data2$SLE, alternative = "two.sided", var.equal = FALSE)