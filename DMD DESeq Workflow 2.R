setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(DESeq2)
library(tidyverse)

studylist<-list.files("DMD_datasets")[1:4]
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

#All Cells (Severe and Moderate CSV)
a<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
b<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
c<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
merge<-rbind(a,b)
merge<-rbind(merge,c)
merge<-data.frame(merge)
significant <- subset(filter(merge, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
significant <- subset(filter(significant, baseMean>100))
significant$abs<-abs(significant$log2FoldChange)
significant<-significant[order(-significant$abs), ]
write.csv(significant, "DMD_datasets/X temp/All_cells_results (Severe).csv")

a<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
b<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
c<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
merge<-rbind(a,b)
merge<-rbind(merge,c)
merge<-data.frame(merge)
significant <- subset(filter(merge, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
significant <- subset(filter(significant, baseMean>100))
significant$abs<-abs(significant$log2FoldChange)
significant<-significant[order(-significant$abs), ]
write.csv(significant, "DMD_datasets/X temp/All_cells_results (Moderate).csv")

#All cells Volcano Plots
library(ggplot2)
library(ggrepel)
severe213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe218<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
moderate213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate771<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/differential_expression_results/PRJNA658118_results.csv")
moderate604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/differential_expression_results/PRJNA604175_results.csv")
human_moderate<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
human_severe<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
  
volcanoplot<-function(a){
  a$delabel <- ifelse(a$gene_name %in% head(a[order(a$padj), "gene_name"], 20), a$gene_name, NA)
  a$diffex <- ifelse(abs(a$log2FoldChange) > .585 & a$padj<.05, "diffex", NA)
  return(ggplot(data = a, aes(x = log2FoldChange, y = -log10(padj),label = delabel, col = diffex)) +
    geom_point()+
    geom_text_repel(max.overlaps = Inf)+
    ggtitle(a$dataset.origin[1]))
}
volcanoplot(severe213)
volcanoplot(severe218)
volcanoplot(moderate213)
volcanoplot(moderate771)
volcanoplot(moderate604)
volcanoplot(moderate658)
volcanoplot(human_moderate)
volcanoplot(human_severe)

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

##########################################################################################################################################################
#UpSetR plots
library(UpSetR)

severe213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe218<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
moderate213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate771<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/differential_expression_results/PRJNA658118_results.csv")
moderate604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/differential_expression_results/PRJNA604175_results.csv")

significant_fun_up <- function (x){
  significant <- subset(filter(x, log2FoldChange>.585))
  significant <- subset(filter(significant, padj<.05))
  return(significant)
}

severe213<-significant_fun_up(severe213)
severe218<-significant_fun_up(severe218)
moderate213<-significant_fun_up(moderate213)
moderate771<-significant_fun_up(moderate771)
moderate604<-significant_fun_up(moderate604)
moderate658<-significant_fun_up(moderate658)

diffexfiles<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658")
allDEGs<-unique(c(severe213$X, severe218$X, moderate213$X, moderate771$X, moderate658$X, moderate604$X))
overlap <- data.frame(matrix(data = 0, ncol = length(diffexfiles), nrow = length(allDEGs)))
rownames(overlap)<-allDEGs
colnames(overlap)<-diffexfiles
for (a in 1:length(allDEGs)){
  if(allDEGs[a] %in% severe213$X){
    overlap$severe213[a]<-1
  }
  if(allDEGs[a] %in% severe218$X){
    overlap$severe218[a]<-1
  }
  if(allDEGs[a] %in% moderate213$X){
    overlap$moderate213[a]<-1
  }
  if(allDEGs[a] %in% moderate771$X){
    overlap$moderate771[a]<-1
  }
  if(allDEGs[a] %in% moderate604$X){
    overlap$moderate604[a]<-1
  }
  if(allDEGs[a] %in% moderate658$X){
    overlap$moderate658[a]<-1
  }
}

plot_up <- upset(overlap, sets = colnames(overlap), order.by = "freq")
plot_up

combos<-overlap%>%group_by_all%>%count
combos<-subset(combos, rowSums(combos[1:6])>1)
combonames<-paste0(combos$severe213,combos$severe218, combos$moderate213, combos$moderate771, combos$moderate604, combos$moderate658)
for(x in 1:length(combonames)){
  combonamex<-as.data.frame(t(str_split(combonames[x], "")[[1]]))
  colnames(combonamex)<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658")
  genelist<-NULL
  for(y in 1:dim(overlap)[1]){
    if(all(combonamex == overlap[y,])){
      genelist<-c(genelist, rownames(overlap[y,]))
    }
  }
  filename<-colnames(combonamex)[combonamex == 1]
  if(length(filename) == 2){
    filename<-paste(filename[1],filename[2])
  }
  if(length(filename) == 3){
    filename<-paste(filename[1],filename[2], filename[3])
  }
  if(length(filename) == 4){
    filename<-paste(filename[1],filename[2], filename[3], filename[4])
  }
  if(length(filename) == 5){
    filename<-paste(filename[1],filename[2], filename[3], filename[4], filename[5])
  }
  cat(genelist, file = paste0("DMD_datasets/X temp/UP", filename, ".txt"))
  df <- data.frame(matrix(data = 0, ncol = 2*sum(as.numeric(combonamex))+1, nrow = length(genelist)))
  df_colnames<-"Gene"
  for(i in 1:length(colnames(combonamex)[combonamex == 1])){
    df_colnames<-c(df_colnames, paste0(colnames(combonamex)[combonamex == 1][i], "_L2FC"), paste0(colnames(combonamex)[combonamex == 1][i], "_padj"))
  }
  colnames(df)<-df_colnames
  df$Gene[1:length(genelist)]<-genelist
  for(j in 1:length(genelist)){
    if("severe213_padj" %in% colnames(df)){
      df[j,"severe213_padj"]<-severe213[severe213$X==genelist[j],"padj"]
      df[j,"severe213_L2FC"]<-severe213[severe213$X==genelist[j],"log2FoldChange"]
    }
    if("moderate213_padj" %in% colnames(df)){
      df[j,"moderate213_padj"]<-moderate213[moderate213$X==genelist[j],"padj"]
      df[j,"moderate213_L2FC"]<-moderate213[moderate213$X==genelist[j],"log2FoldChange"]
    }
    if("severe218_padj" %in% colnames(df)){
      df[j,"severe218_padj"]<-severe218[severe218$X==genelist[j],"padj"]
      df[j,"severe218_L2FC"]<-severe218[severe218$X==genelist[j],"log2FoldChange"]
    }
    if("moderate771_padj" %in% colnames(df)){
      df[j,"moderate771_padj"]<-moderate771[moderate771$X==genelist[j],"padj"]
      df[j,"moderate771_L2FC"]<-moderate771[moderate771$X==genelist[j],"log2FoldChange"]
    }
    if("moderate604_padj" %in% colnames(df)){
      df[j,"moderate604_padj"]<-moderate604[moderate604$X==genelist[j],"padj"]
      df[j,"moderate604_L2FC"]<-moderate604[moderate604$X==genelist[j],"log2FoldChange"]
    }
    if("moderate658_padj" %in% colnames(df)){
      df[j,"moderate658_padj"]<-moderate658[moderate658$X==genelist[j],"padj"]
      df[j,"moderate658_L2FC"]<-moderate658[moderate658$X==genelist[j],"log2FoldChange"]
    }
  }
  write.csv(df, paste0("DMD_datasets/X temp/genelists/UP", filename, "(w-L2FC&padj).csv"))
}

overlap$total<-rowSums(overlap[,1:5])
up_list<-rownames(overlap)[overlap$total>2]
cat(up_list, file = paste0("DMD_datasets/X temp/Up_regulated_genelist.txt"))

allDMD_up<-rownames(overlap[rowSums(overlap)==4,])




moderate213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/differential_expression_results/PRJNA604175_All_results.csv")
moderate658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/differential_expression_results/PRJNA658118_All_results.csv")
moderate771<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate772<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
severe213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe218<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe772<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
moderate213$X<-toupper(moderate213$X)
moderate658$X<-toupper(moderate658$X)
moderate771$X<-toupper(moderate771$X)
moderate772$X<-toupper(moderate772$X)
moderate604$X<-toupper(moderate604$X)
severe213$X<-toupper(severe213$X)
severe218$X<-toupper(severe218$X)
severe772$X<-toupper(severe772$X)

significant_fun_down <- function (x){
  significant <- subset(filter(x, log2FoldChange>.585))
  significant <- subset(filter(significant, padj<.05))
  return(significant)
}

severe213<-significant_fun_down(severe213)
severe218<-significant_fun_down(severe218)
moderate213<-significant_fun_down(moderate213)
moderate771<-significant_fun_down(moderate771)
moderate604<-significant_fun_down(moderate604)
moderate658<-significant_fun_down(moderate658)
moderate772<-significant_fun_down(moderate772)
severe772<-significant_fun_down(severe772)

diffexfiles<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658", "moderate772", "severe772")
allDEGs<-unique(c(severe213$X, severe218$X, moderate213$X, moderate771$X, moderate658$X, moderate604$X, moderate772$X, severe772$X))
overlap <- data.frame(matrix(data = 0, ncol = length(diffexfiles), nrow = length(allDEGs)))
rownames(overlap)<-allDEGs
colnames(overlap)<-diffexfiles
for (a in 1:length(allDEGs)){
  if(allDEGs[a] %in% severe213$X){
    overlap$severe213[a]<-1
  }
  if(allDEGs[a] %in% severe218$X){
    overlap$severe218[a]<-1
  }
  if(allDEGs[a] %in% moderate213$X){
    overlap$moderate213[a]<-1
  }
  if(allDEGs[a] %in% moderate771$X){
    overlap$moderate771[a]<-1
  }
  if(allDEGs[a] %in% moderate604$X){
    overlap$moderate604[a]<-1
  }
  if(allDEGs[a] %in% moderate658$X){
    overlap$moderate658[a]<-1
  }
  if(allDEGs[a] %in% moderate772$X){
    overlap$moderate772[a]<-1
  }
  if(allDEGs[a] %in% severe772$X){
    overlap$severe772[a]<-1
  }
}

#library(grid)
plot_down <- upset(overlap, sets = colnames(overlap), order.by = "freq")#+grid.text("Down Regulated genes",x = 0.65, y=0.95, gp=gpar(fontsize=20))
plot_down






combos<-overlap%>%group_by_all%>%count
combos<-subset(combos, rowSums(combos[1:6])>1)
combonames<-paste0(combos$severe213,combos$severe218, combos$moderate213, combos$moderate771, combos$moderate604, combos$moderate658)
for(x in 1:length(combonames)){
  combonamex<-as.data.frame(t(str_split(combonames[x], "")[[1]]))
  colnames(combonamex)<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658")
  genelist<-NULL
 
  for(y in 1:dim(overlap)[1]){
    if(all(combonamex == overlap[y,])){
      genelist<-c(genelist, rownames(overlap[y,]))
    }
  }
  filename<-colnames(combonamex)[combonamex == 1]
  if(length(filename) == 2){
    filename<-paste(filename[1],filename[2])
  }
  if(length(filename) == 3){
    filename<-paste(filename[1],filename[2], filename[3])
  }
  if(length(filename) == 4){
    filename<-paste(filename[1],filename[2], filename[3], filename[4])
  }
  if(length(filename) == 5){
    filename<-paste(filename[1],filename[2], filename[3], filename[4], filename[5])
  }
  cat(genelist, file = paste0("DMD_datasets/X temp/DOWN", filename, ".txt"))
 
  df <- data.frame(matrix(data = 0, ncol = 2*sum(as.numeric(combonamex))+1, nrow = length(genelist)))
  df_colnames<-"Gene"
  for(i in 1:length(colnames(combonamex)[combonamex == 1])){
    df_colnames<-c(df_colnames, paste0(colnames(combonamex)[combonamex == 1][i], "_L2FC"), paste0(colnames(combonamex)[combonamex == 1][i], "_padj"))
  }
  colnames(df)<-df_colnames
  df$Gene[1:length(genelist)]<-genelist
  for(j in 1:length(genelist)){
    if("severe213_padj" %in% colnames(df)){
      df[j,"severe213_padj"]<-severe213[severe213$X==genelist[j],"padj"]
      df[j,"severe213_L2FC"]<-severe213[severe213$X==genelist[j],"log2FoldChange"]
    }
    if("moderate213_padj" %in% colnames(df)){
      df[j,"moderate213_padj"]<-moderate213[moderate213$X==genelist[j],"padj"]
      df[j,"moderate213_L2FC"]<-moderate213[moderate213$X==genelist[j],"log2FoldChange"]
    }
    if("severe218_padj" %in% colnames(df)){
      df[j,"severe218_padj"]<-severe218[severe218$X==genelist[j],"padj"]
      df[j,"severe218_L2FC"]<-severe218[severe218$X==genelist[j],"log2FoldChange"]
    }
    if("moderate771_padj" %in% colnames(df)){
      df[j,"moderate771_padj"]<-moderate771[moderate771$X==genelist[j],"padj"]
      df[j,"moderate771_L2FC"]<-moderate771[moderate771$X==genelist[j],"log2FoldChange"]
    }
    if("moderate604_padj" %in% colnames(df)){
      df[j,"moderate604_padj"]<-moderate604[moderate604$X==genelist[j],"padj"]
      df[j,"moderate604_L2FC"]<-moderate604[moderate604$X==genelist[j],"log2FoldChange"]
    }
    if("moderate658_padj" %in% colnames(df)){
      df[j,"moderate658_padj"]<-moderate658[moderate658$X==genelist[j],"padj"]
      df[j,"moderate658_L2FC"]<-moderate658[moderate658$X==genelist[j],"log2FoldChange"]
    }
  }
  write.csv(df, paste0("DMD_datasets/X temp/genelists/DOWN", filename, "(w-L2FC&padj).csv"))
}

overlap$total<-rowSums(overlap[,1:5])
up_list<-rownames(overlap)[overlap$total>2]
cat(up_list, file = paste0("DMD_datasets/X temp/Down_regulated_genelist.txt"))

allDMD_down<-rownames(overlap[rowSums(overlap)==4,])


#############################################################################   BULK ANALYSIS   ######################### BULK ANALYSIS #############################
#preparing RareMuscle_DMD_External_PRJNA604175_Bulk_20240110
combine_duplicates<-function(df, duplist){
  new_df<-NULL
  for(a in 1:length(duplist)){
    combo<-c(colSums(df[df$GeneName == duplist[a], c(1,2,3,4)]), duplist[a])#add counts for duplicated genes
    new_df<-rbind(new_df, combo)
  }
  colnames(new_df)[5]<-"GeneName"
  return(new_df[!duplicated(new_df),])
}
data604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/counts_metadata/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110_counts.csv")
meta604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/counts_metadata/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110_metadata.csv")
data604<-data604[,c(1,2,7,8,9)] #remove treatment groups and samples with poor quality
colnames(data604)[c(2,3,4,5)]<-c("DMD1", "control2", "control1", "DMD2") #rename columns
rownames(data604)<-data604$GeneID
data604<-data604[,-1]
data604[,c(1,2,3,4)]<-round(data604[,c(1,2,3,4)]) #round data for DESeq2
data604$GeneName<-meta604$GeneName #will convert ensembl to geneID later
data604$source<-meta604$Source #use to filter for protein coding genes only?
data604<-data604[data604$source == "protein_coding",]
data604<-data604[rowSums(data604[,c(1,2,3,4)])>0,] #filter out genes with no counts
data604<-data604[,-6]
duplicates<-data604$GeneName[duplicated(data604$GeneName)]
dup_df<-data604[data604$GeneName %in% duplicates,]
dup_df<-combine_duplicates(dup_df, duplicates)
data604<-data604[!(data604$GeneName %in% duplicates),]
data604<-rbind(dup_df, data604)
rownames(data604)<-data604$GeneName
data604<-data604[,-5]
i <- c(1, 2, 3, 4)
data604[ , i] <- apply(data604[ , i], 2,function(x) as.numeric(as.character(x)))
write.csv(data604, "DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/counts_metadata/cleaned_counts.csv")
metadata604 <- data.frame(matrix(data = 0, ncol = 3, nrow = 4))
rownames(metadata604)<-colnames(data604)
colnames(metadata604)<-c("condition", "Age (months)","tissue")
metadata604["DMD1",]<-c("DMD (moderate)", 6, "Left Ventricle")
metadata604["DMD2",]<-c("DMD (moderate)", 6, "Left Ventricle")
metadata604["control1",]<-c("control", 6, "Left Ventricle")
metadata604["control2",]<-c("control", 6, "Left Ventricle")
write.csv(metadata604, "DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/counts_metadata/cleaned_metadata.csv")
data_counts<-data604
metadata<-metadata604
dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
dds <- dds[keep,]
dds <- DESeq(dds)
result <- results(dds, contrast = c('condition', 'DMD (moderate)', 'control'))
final_results <- data.frame(result)
final_results$gene_name <- rownames(final_results)
final_results$comparison <- "DMD (Moderate) vs Control"
write.csv(final_results, paste0("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/differential_expression_results/PRJNA604175_results.csv"))



#preparing RareMuscle_DMD_External_PRJNA658118_Bulk_20240122
combine_duplicates2<-function(df, duplist){
  new_df<-NULL
  for(a in 1:length(duplist)){
    combo<-c(colSums(df[df$GeneName == duplist[a], c(1,2,3,4,5,6)]), duplist[a]) #make this colMeans
    new_df<-rbind(new_df, combo)
  }
  colnames(new_df)[7]<-"GeneName"
  return(new_df[!duplicated(new_df),])
}
data658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/counts_metadata/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122_counts.csv")
meta658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/counts_metadata/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122_metadata.csv")
colnames(data658)[c(2,3,4,5,6,7)]<-c("TA_WT_2 RNA-seq", "TA_WT_1 RNA-seq", "TA_D51_4 RNA-seq", "TA_D51_5 RNA-seq", "TA_D51_6 RNA-seq", "TA_WT_3 RNA-seq") #rename columns
rownames(data658)<-data658$GeneID
data658<-data658[,-1]
data658[,c(1,2,3,4,5,6)]<-round(data658[,c(1,2,3,4,5,6)]) #round data for DESeq2
data658$GeneName<-meta658$GeneName #will convert ensembl to geneID later
data658$source<-meta658$Source #use to filter for protein coding genes only?
data658<-data658[data658$source == "protein_coding",]
data658<-data658[rowSums(data658[,c(1,2,3,4,5,6)])>0,] #filter out genes with no counts
data658<-data658[,-8]
duplicates<-data658$GeneName[duplicated(data658$GeneName)]
dup_df<-data658[data658$GeneName %in% duplicates,]
dup_df<-combine_duplicates2(dup_df, duplicates)
data658<-data658[!(data658$GeneName %in% duplicates),]
data658<-rbind(dup_df, data658)
rownames(data658)<-data658$GeneName
data658<-data658[,-7]
i <- c(1, 2, 3, 4,5,6)
data658[ , i] <- apply(data658[ , i], 2,function(x) as.numeric(as.character(x)))
write.csv(data658, "DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/counts_metadata/cleaned_counts.csv")
metadata658 <- data.frame(matrix(data = 0, ncol = 3, nrow = 6))
rownames(metadata658)<-colnames(data658)
colnames(metadata658)<-c("condition", "Age (months)","tissue")
metadata658["TA_WT_2 RNA-seq",]<-c("control", 1, "Tibialis")
metadata658["TA_WT_1 RNA-seq",]<-c("control", 1, "Tibialis")
metadata658["TA_WT_3 RNA-seq",]<-c("control", 1, "Tibialis")
metadata658["TA_D51_4 RNA-seq",]<-c("DMD (moderate)", 1, "Tibialis")
metadata658["TA_D51_5 RNA-seq",]<-c("DMD (moderate)", 1, "Tibialis")
metadata658["TA_D51_6 RNA-seq",]<-c("DMD (moderate)", 1, "Tibialis")
write.csv(metadata658, "DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/counts_metadata/cleaned_metadata.csv")
data_counts<-data658
metadata<-metadata658
dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
dds <- dds[keep,]
dds <- DESeq(dds)
result <- results(dds, contrast = c('condition', 'DMD (moderate)', 'control'))
final_results <- data.frame(result)
final_results$gene_name <- rownames(final_results)
final_results$comparison <- "DMD (Moderate) vs Control"
write.csv(final_results, paste0("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/differential_expression_results/PRJNA658118_results.csv"))



############################### MICROARRAY DATA #################################################
files<-list.files("DMD_datasets/X temp/MicroArray")
moderate109<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[1]))
severe109<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[2]))
moderate136<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[3]))
severe136<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[4]))
moderate384<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[5]))
severe384<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[6]))
moderate601<-read_tsv(paste0("DMD_datasets/X temp/MicroArray/", files[7]))

clean<-function(x){
  x <- subset(x, !(is.na(x$Gene.symbol)))
  x <- subset(filter(x, abs(logFC)>.585))
  x <- subset(filter(x, adj.P.Val<.05))
  return(x[,c("adj.P.Val", "logFC", "Gene.symbol")])
}

moderate109<-clean(moderate109)
severe109<-clean(severe109)
moderate136<-clean(moderate136)
severe136<-clean(severe136)
moderate384<-clean(moderate384)
severe384<-clean(severe384)
moderate601<-clean(moderate601)

remove_duplicates<-function(x){
  x$Gene.symbol<-toupper(x$Gene.symbol)
  if(TRUE %in% (duplicated(x$Gene.symbol))){
    duplicates<-unique(x$Gene.symbol[duplicated(x$Gene.symbol)])
    dup_df<-x[x$Gene.symbol %in% duplicates,]
    new_df<-NULL
    for(a in 1:length(duplicates)){
      combo<-c(colMeans(dup_df[dup_df$Gene.symbol == duplicates[a], c(1,2)]), duplicates[a])
      new_df<-rbind(new_df, combo)
    }
    colnames(new_df)[3]<-"Gene.symbol"
    x<-x[!(x$Gene.symbol %in% duplicates),]
    x<-rbind(new_df, x)
    rownames(x)<-x$Gene.symbol
    i <- c(1,2)
    x[ , i] <- apply(x[ , i], 2,function(x) (as.numeric(as.character(x)))) 
    x[,2]<- -x[,2]#flip l2fc so it is DMD vs control
    return(x)
  }
  return(x)
}

moderate109<-remove_duplicates(moderate109)
severe109<-remove_duplicates(severe109)
moderate136<-remove_duplicates(moderate136)
severe136<-remove_duplicates(severe136)
moderate384<-remove_duplicates(moderate384)
severe384<-remove_duplicates(severe384)
moderate601<-remove_duplicates(moderate601)

write.csv(moderate109, "DMD_datasets/X temp/MicroArray/moderate109_clean_diffex.csv")
write.csv(moderate136, "DMD_datasets/X temp/MicroArray/moderate136_clean_diffex.csv")
write.csv(moderate384, "DMD_datasets/X temp/MicroArray/moderate384_clean_diffex.csv")
write.csv(moderate601, "DMD_datasets/X temp/MicroArray/moderate601_clean_diffex.csv")
write.csv(severe109, "DMD_datasets/X temp/MicroArray/severe109_clean_diffex.csv")
write.csv(severe136, "DMD_datasets/X temp/MicroArray/severe136_clean_diffex.csv")
write.csv(severe384, "DMD_datasets/X temp/MicroArray/severe384_clean_diffex.csv")


##### microArray heatmap
library(pheatmap)
allDEGs<-unique(c(rownames(moderate109),rownames(moderate136), rownames(moderate384), rownames(moderate601), rownames(severe109), rownames(severe136), rownames(severe384)))
diffexfiles<-c("moderate109", "moderate136", "moderate384", "moderate601", "severe109", "severe136", "severe384")
overlap <- data.frame(matrix(data = 0, ncol = length(diffexfiles), nrow = length(allDEGs)))
rownames(overlap)<-allDEGs
colnames(overlap)<-diffexfiles
for (a in 1:length(allDEGs)){
  if(allDEGs[a] %in% moderate109$Gene.symbol){
    overlap$moderate109[a]<-moderate109$logFC[moderate109$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate136$Gene.symbol){
    overlap$moderate136[a]<-moderate136$logFC[moderate136$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate384$Gene.symbol){
    overlap$moderate384[a]<-moderate384$logFC[moderate384$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate601$Gene.symbol){
    overlap$moderate601[a]<-moderate601$logFC[moderate601$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe109$Gene.symbol){
    overlap$severe109[a]<-severe109$logFC[severe109$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe136$Gene.symbol){
    overlap$severe136[a]<-severe136$logFC[severe136$Gene.symbol == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe384$Gene.symbol){
    overlap$severe384[a]<-severe384$logFC[severe384$Gene.symbol == allDEGs[a]]
  }
}
overlap<-overlap[rowSums(overlap==0)<=3,]
overlap$keep<-apply(overlap[,1:7], 1, median, na.rm=T)
overlap<-overlap[abs(overlap$keep)>2.5,]
overlap<-overlap[,1:7]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(overlap), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(overlap)/paletteLength, max(overlap), length.out=floor(paletteLength/2)))
if(TRUE %in% duplicated(myBreaks)){ # If the lowest L2FC is positive or the highest is negative
  duplicate_length<-sum(duplicated(myBreaks))+1
  if(myBreaks[length(myBreaks)] == 0){ #if the highest is negative
    myBreaks[duplicated(myBreaks)]<-seq(from=.1/sum(duplicated(myBreaks)), to=.1, length.out=sum(duplicated(myBreaks)))
  }
  if(myBreaks[1] == 0){ #if the lowest is positive
    myBreaks[duplicated(myBreaks)]<-seq(from=-.1, to=-.1/sum(duplicated(myBreaks)), length.out=sum(duplicated(myBreaks)))
    hold<-duplicate_length+1
    myBreaks<-c(myBreaks[2:duplicate_length],0,myBreaks[hold:length(myBreaks)])
  }
}
ph <- pheatmap(overlap,
               color=myColor,
               breaks=myBreaks,
               main = paste('Log2FoldChange'))
ggsave(filename = paste0("MicroArray_heatmap.png"),
       plot = ph,
       device = 'png',
       path = paste0("DMD_datasets/X temp/Heatmaps"),
       width = 8,
       height = 15,
       units = "in",
       limitsize = FALSE)