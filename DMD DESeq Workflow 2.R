setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(DESeq2)
library(tidyverse)

studylist<-list.files("DMD_datasets")
odds <- function(x) x[ x %% 2 == 1 ]

for(s in 1:length(studylist)){
  study <- studylist[s]
  files <- list.files(path= paste0('DMD_datasets/', study, '/pseudobulk_counts_matrices'))
  
  modifier <- odds(1:length(files))
  for (j in modifier){
    data_counts<-as.data.frame(read.csv(paste0('DMD_datasets/', study, '/pseudobulk_counts_matrices/', files[j])))
    rownames(data_counts)<-data_counts$X
    data_counts<-data_counts[-1]
    data_counts <- round(data_counts)
    metadata<-as.data.frame(read.csv(paste0('DMD_datasets/', study, '/pseudobulk_counts_matrices/', files[j+1])))
    rownames(metadata) <- colnames(data_counts)
    metadata$sample <- rownames(metadata)
    
    if(length(metadata$sample)>3){
      #run DESeq2
      if(length(unique(metadata$gender))>1){
        dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition + gender)
      }else{
        dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
      }
      dds <- estimateSizeFactors(dds)
      keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
      dds <- dds[keep,]
      dds <- DESeq(dds)
      
      #define disease/treatment states and compare each one to control
      for (i in 1:length(unique(metadata$condition))){
        disease <- as.character(unique(metadata$condition)[i])
        if(disease!='control'){
          result <- results(dds, contrast = c('condition', disease, 'control')) 
          
          #Convert Results to dataframe and add celltype identifier and comparison identifier
          final_results <- data.frame(result)
          final_results$celltype <- gsub(".csv", "", gsub("_pseudobulk_data", "", files[j]))
          #final_results$celltype <- "stalk_like"
          final_results$gene_name <- rownames(final_results)
          final_results$comparison <- gsub(" ", "_", paste(disease, "vs.",'Control'))
          final_results$dataset.origin <- study
          
          #Export results as csv to folder
          write.csv(final_results, paste0("DMD_datasets/", study, "/differential_expression_results/diffex_results_", final_results$celltype[1], final_results$comparison[1], ".csv"))
        }
      }
      counts <- as.data.frame(counts(dds, normalized = TRUE))
      counts$gene_name <- rownames(counts)
      write.csv(counts, paste0("DMD_datasets/", study, "/normalized_pseudobulk_counts_matrices/", final_results$celltype[1], ".csv"))
    }
  }
}

#All Cells CSV
a<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
b<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
c<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
merge<-rbind(a,b)
merge<-rbind(merge,c)
merge<-data.frame(merge)
significant <- subset(filter(merge, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
significant <- subset(filter(significant, baseMean>150))
significant$abs<-abs(significant$log2FoldChange)
significant<-significant[order(-significant$abs), ]
write.csv(significant, "DMD_datasets/All_cells_results.csv")

#All cells Volcano Plots
library(ggplot2)
library(ggrepel)
a<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
b<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
c<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")

a$delabel <- ifelse(a$gene_name %in% head(a[order(a$padj), "gene_name"], 50), a$gene_name, NA)
a$diffex <- ifelse(abs(a$log2FoldChange) > .585 & a$padj<.05, "diffex", NA)

b$delabel <- ifelse(b$gene_name %in% head(b[order(b$padj), "gene_name"], 50), b$gene_name, NA)
b$diffex <- ifelse(abs(b$log2FoldChange) > .585 & b$padj<.05, "diffex", NA)

c$delabel <- ifelse(c$gene_name %in% head(c[order(c$padj), "gene_name"], 50), c$gene_name, NA)
c$diffex <- ifelse(abs(c$log2FoldChange) > .585 & c$padj<.05, "diffex", NA)

myvolcanoplota <- ggplot(data = a, aes(x = log2FoldChange, y = -log10(padj),label = delabel, col = diffex)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(a$dataset.origin[1])

myvolcanoplotb <- ggplot(data = b, aes(x = log2FoldChange, y = -log10(padj),label = delabel, col = diffex)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(b$dataset.origin[1])

myvolcanoplotc <- ggplot(data = c, aes(x = log2FoldChange, y = -log10(padj),label = delabel, col = diffex)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(c$dataset.origin[1])

myvolcanoplota
myvolcanoplotb
myvolcanoplotc

tripleoverlap<-c("LGALS3" , "ACKR3",   "APOE"  ,  "PDGFD" ,  "JAM2"   , "CDK6"  ,  "TIMP3"  , "RUNX1"  , "PTPRC"  , "DOCK8"  , "EBF2" ,  "VEGFA" ,  "CALCRL" , "SLC7A2",  "CD74"  ,  "COL4A5" , "CD44" ,   "PIP4K2A", "RETREG1", "COL1A1" , "SPTAN1",  "ATP1A1") 

a$GOI <- ifelse(a$gene_name %in% tripleoverlap, a$gene_name, NA)
b$GOI <- ifelse(b$gene_name %in% tripleoverlap, b$gene_name, NA)
c$GOI <- ifelse(c$gene_name %in% tripleoverlap, c$gene_name, NA)

myvolcanoplota2 <- ggplot(data = a, aes(x = log2FoldChange, y = -log10(padj),label = GOI, col = GOI)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(a$dataset.origin[1])

myvolcanoplotb2 <- ggplot(data = b, aes(x = log2FoldChange, y = -log10(padj),label = GOI, col = GOI)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(b$dataset.origin[1])

myvolcanoplotc2 <- ggplot(data = c, aes(x = log2FoldChange, y = -log10(padj),label = GOI, col = GOI)) +
  geom_point()+
  geom_text_repel(max.overlaps = Inf)+
  ggtitle(c$dataset.origin[1])

myvolcanoplota2
myvolcanoplotb2
myvolcanoplotc2

#All cells Venn Diagram results
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")

a<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
b<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")
c<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_vs._Control.csv")

a<-subset(filter(a, abs(log2FoldChange)>.585))
a<-subset(filter(a, padj<.05))
a<-subset(filter(a, baseMean>150))
a$abs<-abs(a$log2FoldChange)
a<-a[order(-a$abs), ]
alist<-as.character(a$gene_name)

apos<-subset(filter(a, log2FoldChange>.585))
apos<-apos[order(-apos$log2FoldChange), ]
apos<-as.character(apos$gene_name)
aneg<-subset(filter(a, log2FoldChange< -.585))
aneg<-aneg[order(aneg$log2FoldChange), ]
aneg<-as.character(aneg$gene_name)

b<-subset(filter(b, abs(log2FoldChange)>.585))
b<-subset(filter(b, padj<.05))
b<-subset(filter(b, baseMean>150))
b$abs<-abs(b$log2FoldChange)
b<-b[order(-b$abs), ]
blist<-as.character(b$gene_name)

bpos<-subset(filter(b, log2FoldChange>.585))
bpos<-bpos[order(-bpos$log2FoldChange), ]
bpos<-as.character(bpos$gene_name)
bneg<-subset(filter(b, log2FoldChange< -.585))
bneg<-bneg[order(bneg$log2FoldChange), ]
bneg<-as.character(bneg$gene_name)

c<-subset(filter(c, abs(log2FoldChange)>.585))
c<-subset(filter(c, padj<.05))
c<-subset(filter(c, baseMean>150))
c$abs<-abs(c$log2FoldChange)
c<-c[order(-c$abs), ]
clist<-as.character(c$gene_name)

cpos<-subset(filter(c, log2FoldChange>.585))
cpos<-cpos[order(-cpos$log2FoldChange), ]
cpos<-as.character(cpos$gene_name)
cneg<-subset(filter(c, log2FoldChange< -.585))
cneg<-cneg[order(cneg$log2FoldChange), ]
cneg<-as.character(cneg$gene_name)

x <- list(
  GSE218201 = alist, 
  PRJNA771932 = blist, 
  PRJNA772047 = clist)

xpos <- list(
  GSE218201 = apos, 
  PRJNA771932 = bpos, 
  PRJNA772047 = cpos)

xneg <- list(
  GSE218201 = aneg, 
  PRJNA771932 = bneg, 
  PRJNA772047 = cneg)

venn.plot <- ggvenn(x)
venn.pos <- ggvenn(xpos)
venn.neg <- ggvenn(xneg)

library(gplots)
v.table <- venn(x)
v.neg<-venn(xneg)
v.pos<-venn(xpos)
print(v.table)
print(v.pos)
print(v.neg)