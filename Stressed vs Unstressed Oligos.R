setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(Seurat)
library(tidyverse)
library(DESeq2)

seurat.obj<-readRDS("OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj@meta.data$keep<-FALSE
seurat.obj@meta.data$keep[seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte. %in% c("stressed oligodendrocyte", "normal oligodendrocyte")]<-TRUE
Idents(seurat.obj)<-seurat.obj@meta.data$keep
seurat.obj<-subset(seurat.obj, idents = TRUE)
seurat.obj@meta.data$patient_info<-paste(seurat.obj@meta.data$Age..years., seurat.obj@meta.data$gender..standardized., seurat.obj@meta.data$Clinical.notes)
seurat.obj@meta.data$sample_id<-paste(seurat.obj@meta.data$Sample.ID, seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.)
seurat.obj@meta.data$gender<-seurat.obj@meta.data$Gender
seurat.obj@meta.data$condition<-seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.
seurat.obj@meta.data$cell_type<-seurat.obj@meta.data$Condition #Tissue but is put as cell type to be compatible with the pseudobulk pipeline

print(unique(seurat.obj@meta.data$sample_id))
print(unique(seurat.obj@meta.data$patient_info))
print(unique(seurat.obj@meta.data$gender))
print(unique(seurat.obj@meta.data$condition))
print(unique(seurat.obj@meta.data$cell_type))

metadata2<-seurat.obj@meta.data
celltypes_list<-unique(metadata2$cell_type)

#split by cell type
for(x in 1:length(celltypes_list)){
  celltype_rownames <- rownames(subset(metadata2, metadata2$cell_type==celltypes_list[x]))
  if(length(celltype_rownames)>200){
    celltype_counts <- seurat.obj@assays$RNA@counts[,colnames(seurat.obj@assays$RNA@counts) %in% celltype_rownames]
    celltype_metadata <- metadata2[rownames(metadata2) %in% celltype_rownames,]
    sample_ids <- unique(celltype_metadata$sample_id)
    
    #split by sample,
    output_counts <- data.frame()
    output_metadata <- data.frame()
    for(y in 1:length(sample_ids)){
      sample_ids_rownames <- rownames(subset(celltype_metadata, celltype_metadata$sample_id==sample_ids[y]))
      grouped_celltype_counts <- as.data.frame(celltype_counts[,colnames(celltype_counts) %in% sample_ids_rownames])
      grouped_celltype_metadata <- celltype_metadata[rownames(celltype_metadata) %in% sample_ids_rownames,]
      
      #combine within sample
      if(sum(grouped_celltype_counts[1])>1){
        grouped_celltype_counts <- as.data.frame(rowSums(grouped_celltype_counts))
        colnames(grouped_celltype_counts) <-sample_ids[y]
        grouped_celltype_metadata <- grouped_celltype_metadata[1,]
        rownames(grouped_celltype_metadata) <- grouped_celltype_metadata$sample_id
        if(length(output_counts)==0){
          output_counts <- grouped_celltype_counts
        }else{
          output_counts <- cbind(output_counts,grouped_celltype_counts)
        }
        if(length(output_metadata)==0){
          output_metadata <- grouped_celltype_metadata
        }else{
          output_metadata <- rbind(output_metadata,grouped_celltype_metadata)
        }
      }
    }
    data_counts<-output_counts
    data_counts <- round(data_counts)
    metadata<-output_metadata
    metadata$sample <- rownames(metadata)
    dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition + patient_info)
    dds <- estimateSizeFactors(dds)
    keep <- rowSums(dds@colData$condition == 'normal oligodendrocyte' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'normal oligodendrocyte') & rowSums(dds@colData$condition != 'normal oligodendrocyte' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'normal oligodendrocyte')
    dds <- dds[keep,]
    dds <- DESeq(dds)
    result <- results(dds, contrast = c('condition', 'stressed oligodendrocyte', 'normal oligodendrocyte'))
    final_results <- data.frame(result)
    final_results$gene_name <- rownames(final_results)
    final_results$comparison <- "Stressed vs normal"
    final_results$tissue<-celltypes_list[x]
    write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_(stressed.vs.normal oligos gse180759 - patient controlled " , final_results$tissue[1], ").csv"))
  }
}


####################################################ALL CELLS#############################################################################
sample_ids <- unique(metadata2$sample_id)

#split by sample
output_counts <- data.frame()
output_metadata <- data.frame()
for(y in 1:length(sample_ids)){
  sample_ids_rownames <- rownames(subset(metadata2, metadata2$sample_id==sample_ids[y]))
  grouped_sample_counts <- as.matrix(seurat.obj$RNA@counts[,colnames(seurat.obj$RNA@counts) %in% sample_ids_rownames])
  grouped_sample_metadata <- metadata2[rownames(metadata2) %in% sample_ids_rownames,]
  
  #combine within sample
  grouped_sample_counts <- as.data.frame(rowSums(grouped_sample_counts))
  colnames(grouped_sample_counts) <-sample_ids[y]
  grouped_sample_metadata <- grouped_sample_metadata[1,]
  rownames(grouped_sample_metadata) <- grouped_sample_metadata$sample_id
  if(length(output_counts)==0){
    output_counts <- grouped_sample_counts
  }else{
    output_counts <- cbind(output_counts,grouped_sample_counts)
  }
  if(length(output_metadata)==0){
    output_metadata <- grouped_sample_metadata
  }else{
    output_metadata <- rbind(output_metadata,grouped_sample_metadata)
  }
}
data_counts<-output_counts
data_counts <- round(data_counts)
metadata<-output_metadata
metadata$sample <- rownames(metadata)
dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition + patient_info)
dds <- estimateSizeFactors(dds)
keep <- rowSums(dds@colData$condition == 'normal oligodendrocyte' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'normal oligodendrocyte') & rowSums(dds@colData$condition != 'normal oligodendrocyte' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'normal oligodendrocyte')
dds <- dds[keep,]
dds <- DESeq(dds)
result <- results(dds, contrast = c('condition', 'stressed oligodendrocyte', 'normal oligodendrocyte'))
final_results <- data.frame(result)
final_results$gene_name <- rownames(final_results)
final_results$comparison <- "Stressed vs normal"
write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_(stressed.vs.normal oligos gse180759 - patient controlled All tissues).csv"))

#stressed vs normal quantification
healthy_normal<-sum(metadata2$Clinical.notes == "Unassigned" & metadata2$condition == "normal oligodendrocyte")
MS_normal<-sum(metadata2$Clinical.notes != "Unassigned" & metadata2$condition == "normal oligodendrocyte")
healthy_stressed<-sum(metadata2$Clinical.notes == "Unassigned" & metadata2$condition != "normal oligodendrocyte")
MS_stressed<-sum(metadata2$Clinical.notes != "Unassigned" & metadata2$condition != "normal oligodendrocyte")
healthy_normal
healthy_stressed
MS_normal
MS_stressed

#MS normal vs healthy normal or stressed vs stressed
seurat.obj<-readRDS("OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj@meta.data$keep<-FALSE
seurat.obj@meta.data$keep[seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.== "stressed oligodendrocyte"]<-TRUE  #change this for stressed vs stressed
Idents(seurat.obj)<-seurat.obj@meta.data$keep
seurat.obj<-subset(seurat.obj, idents = TRUE)
seurat.obj@meta.data$patient_info<-paste(seurat.obj@meta.data$Age..years., seurat.obj@meta.data$gender..standardized., seurat.obj@meta.data$Clinical.notes)
seurat.obj@meta.data$sample_id<-paste(seurat.obj@meta.data$Sample.ID, seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.)
seurat.obj@meta.data$gender<-seurat.obj@meta.data$Gender
seurat.obj@meta.data$condition[seurat.obj@meta.data$Clinical.notes == "Unassigned"]<-"control"
seurat.obj@meta.data$condition[seurat.obj@meta.data$Clinical.notes != "Unassigned"]<-"MS"
seurat.obj@meta.data$cell_type<-seurat.obj@meta.data$Condition #Tissue but is put as cell type to be compatible with the pseudobulk pipeline

print(unique(seurat.obj@meta.data$sample_id))
print(unique(seurat.obj@meta.data$patient_info))
print(unique(seurat.obj@meta.data$gender))
print(unique(seurat.obj@meta.data$condition))
print(unique(seurat.obj@meta.data$cell_type))

metadata2<-seurat.obj@meta.data

sample_ids <- unique(metadata2$sample_id)

#split by sample
output_counts <- data.frame()
output_metadata <- data.frame()
for(y in 1:length(sample_ids)){
  sample_ids_rownames <- rownames(subset(metadata2, metadata2$sample_id==sample_ids[y]))
  grouped_sample_counts <- as.matrix(seurat.obj$RNA@counts[,colnames(seurat.obj$RNA@counts) %in% sample_ids_rownames])
  grouped_sample_metadata <- metadata2[rownames(metadata2) %in% sample_ids_rownames,]
  
  #combine within sample
  grouped_sample_counts <- as.data.frame(rowSums(grouped_sample_counts))
  colnames(grouped_sample_counts) <-sample_ids[y]
  grouped_sample_metadata <- grouped_sample_metadata[1,]
  rownames(grouped_sample_metadata) <- grouped_sample_metadata$sample_id
  if(length(output_counts)==0){
    output_counts <- grouped_sample_counts
  }else{
    output_counts <- cbind(output_counts,grouped_sample_counts)
  }
  if(length(output_metadata)==0){
    output_metadata <- grouped_sample_metadata
  }else{
    output_metadata <- rbind(output_metadata,grouped_sample_metadata)
  }
}
data_counts<-output_counts
data_counts <- round(data_counts)
metadata<-output_metadata
metadata$sample <- rownames(metadata)
dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition) #+ patient_info
dds <- estimateSizeFactors(dds)
keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
dds <- dds[keep,]
dds <- DESeq(dds)
result <- results(dds, contrast = c('condition', 'MS', 'control'))
final_results <- data.frame(result)
final_results$gene_name <- rownames(final_results)
final_results$comparison <- "stressed vs stressed oligos in MS vs control"
write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_(stressed vs stressed in MS vs control oligos gse180759 All tissues).csv"))




#Volcano plots for slides
volcanoplot<-function(a){
  a$delabel <- ifelse(a$gene_name %in% head(a[order(a$padj), "gene_name"], 20), a$gene_name, NA)
  a$diffex <- ifelse(abs(a$log2FoldChange) > .585 & a$padj<.05, "diffex", NA)
  return(ggplot(data = a, aes(x = log2FoldChange, y = -log10(padj),label = delabel, col = diffex)) +
           geom_point()+
           geom_text_repel(max.overlaps = Inf)+
           ggtitle(a$dataset.origin[1]))
}

stressedvsnormal<-list.files("OligoAnalysis/differential_expression_results/")[4:10]
for(x in 1:length(stressedvsnormal)){
  a<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[x]))
  volcanoplot(a)
  print(stressedvsnormal[x])
}

#pairwise group info
library(Seurat)
seurat.obj<-readRDS("OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj@meta.data$keep<-FALSE
seurat.obj@meta.data$keep[seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte %in% c("normal oligodendrocyte", "stressed oligodendrocyte")]<-TRUE  #change this for stressed vs stressed
Idents(seurat.obj)<-seurat.obj@meta.data$keep
seurat.obj<-subset(seurat.obj, idents = TRUE)
seurat.obj@meta.data$stressed<-seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.
seurat.obj@meta.data$condition[seurat.obj@meta.data$Clinical.notes == "Unassigned"]<-"control"
seurat.obj@meta.data$condition[seurat.obj@meta.data$Clinical.notes != "Unassigned"]<-"MS"
seurat.obj@meta.data$tissue<-seurat.obj@meta.data$Condition

print(unique(seurat.obj@meta.data$stressed))
print(unique(seurat.obj@meta.data$condition))
print(unique(seurat.obj@meta.data$tissue))

print(unique(seurat.obj@meta.data$stressed[seurat.obj@meta.data$condition == "control"]))
print(unique(seurat.obj@meta.data$condition[seurat.obj@meta.data$condition == "control"]))
print(unique(seurat.obj@meta.data$tissue[seurat.obj@meta.data$condition == "control"]))

stressed_MS_ChronicActive <-sum(seurat.obj@meta.data$stressed == "stressed oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Chronic active MS lesion edge")
stressed_MS_ChronicInactive <-sum(seurat.obj@meta.data$stressed == "stressed oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Chronic inactive MS lesion edge")
stressed_MS_MSPeriplaqueWM <-sum(seurat.obj@meta.data$stressed == "stressed oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "MS periplaque white matter")
stressed_MS_NormalWM <-sum(seurat.obj@meta.data$stressed == "stressed oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Normal control white matter")
stressed_Ctrl_NormalWM <-sum(seurat.obj@meta.data$stressed == "stressed oligodendrocyte" & seurat.obj@meta.data$condition == "control" & seurat.obj@meta.data$tissue == "Normal control white matter")
Normal_MS_ChronicActive <-sum(seurat.obj@meta.data$stressed == "normal oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Chronic active MS lesion edge")
Normal_MS_ChronicInactive <-sum(seurat.obj@meta.data$stressed == "normal oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Chronic inactive MS lesion edge")
Normal_MS_MSPeriplaqueWM <-sum(seurat.obj@meta.data$stressed == "normal oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "MS periplaque white matter")
Normal_MS_NormalWM <-sum(seurat.obj@meta.data$stressed == "normal oligodendrocyte" & seurat.obj@meta.data$condition == "MS" & seurat.obj@meta.data$tissue == "Normal control white matter")
Normal_Ctrl_NormalWM <-sum(seurat.obj@meta.data$stressed == "normal oligodendrocyte" & seurat.obj@meta.data$condition == "control" & seurat.obj@meta.data$tissue == "Normal control white matter")

#upsetr plots
library(UpSetR)

significant_fun_up <- function (x){
  significant <- subset(filter(x, log2FoldChange>.585)) #< -.585
  significant <- subset(filter(significant, padj<.05))
  return(significant)
}

stressedvsnormal<-list.files("OligoAnalysis/differential_expression_results/")[17:21]
a<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[1]))
a<-significant_fun_up(a)
b<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[2]))
b<-significant_fun_up(b)
c<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[3]))
c<-significant_fun_up(c)
d<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[4]))
d<-significant_fun_up(d)
e<-read.csv(paste0("OligoAnalysis/differential_expression_results/", stressedvsnormal[5]))
e<-significant_fun_up(e)

diffexfiles<-stressedvsnormal

allDEGs<-unique(c(a$X, b$X, c$X, d$X, e$X))
overlap <- data.frame(matrix(data = 0, ncol = length(diffexfiles), nrow = length(allDEGs)))
rownames(overlap)<-allDEGs
for (x in 1:length(allDEGs)){
  if(allDEGs[x] %in% a$X){
    overlap$X1[x]<-1
  }
  if(allDEGs[x] %in% b$X){
    overlap$X2[x]<-1
  }
  if(allDEGs[x] %in% c$X){
    overlap$X3[x]<-1
  }
  if(allDEGs[x] %in% d$X){
    overlap$X4[x]<-1
  }
  if(allDEGs[x] %in% e$X){
    overlap$X5[x]<-1
  }
}
colnames(overlap)<-diffexfiles

plot_up <- upset(overlap, sets = colnames(overlap), order.by = "freq")
plot_up

overlap$total<-rowSums(overlap[,1:5])
up_list<-rownames(overlap)[overlap$total==5]
cat(up_list, file = paste0("OligoAnalysis/Up_regulated_genelist.txt"))


Idents(seurat.obj)<-seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.
seurat.obj <- FindVariableFeatures(object = seurat.obj)
seurat.obj <- ScaleData(object = seurat.obj)
seurat.obj <- RunPCA(object = seurat.obj)
seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:15)
DimPlot(seurat.obj, reduction = 'umap', group.by = "Author.s.cell.subtype..Oligodendrocyte.", raster = FALSE)

#looking for stressed profile in other datasets
seurat.obj2<-readRDS("OligoAnalysis/seurat_object/GSE118257.rds")
obj2meta<-as.data.frame(read_tsv("OligoAnalysis/seurat_object/GSE118257.metadata (2).tsv"))
seurat.obj2@meta.data$Up_regulated_Enrichment<-obj2meta$`Up_Regulated_Stressed_Oligos - AUCell score`
seurat.obj2@meta.data$Down_regulated_Enrichment<-obj2meta$`Down_Regulated_Stressed_Oligos - AUCell score`
Idents(seurat.obj2)<-seurat.obj2@meta.data$cell.type...subgroup..standardized.
seurat.obj2<-seurat.obj2<-subset(seurat.obj2, ident = "oligodendrocyte")
mean(seurat.obj2@meta.data$Up_regulated_Enrichment[seurat.obj2@meta.data$Condition == "healthy control"])
mean(seurat.obj2@meta.data$Up_regulated_Enrichment[seurat.obj2@meta.data$Condition == "Multiple sclerosis (MS)"])
mean(seurat.obj2@meta.data$Down_regulated_Enrichment[seurat.obj2@meta.data$Condition == "healthy control"])
mean(seurat.obj2@meta.data$Down_regulated_Enrichment[seurat.obj2@meta.data$Condition == "Multiple sclerosis (MS)"])
seurat.obj2 <- FindVariableFeatures(object = seurat.obj2)
seurat.obj2 <- ScaleData(object = seurat.obj2)
seurat.obj2 <- RunPCA(object = seurat.obj2)
seurat.obj2 <- RunUMAP(object = seurat.obj2, dims = 1:15)
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj2@meta.data$Up_regulated_Enrichment)))
FeaturePlot(seurat.obj2, cols= myColor, reduction = "umap", features = "Up_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj2@meta.data$Down_regulated_Enrichment)))
FeaturePlot(seurat.obj2, cols= myColor, reduction = "umap", features = "Down_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj2@meta.data$Condition)))
DimPlot(seurat.obj2, cols = myColor, reduction = 'umap', group.by = "Condition", raster = FALSE)

seurat.obj3<-readRDS("OligoAnalysis/seurat_object/PMID31316211.rds")
obj3meta<-as.data.frame(read_tsv("OligoAnalysis/seurat_object/PMID31316211.metadata (2).tsv"))
seurat.obj3@meta.data$Up_regulated_Enrichment<-obj3meta$`Up_Regulated_Stressed_Oligos - AUCell score`
seurat.obj3@meta.data$Down_regulated_Enrichment<-obj3meta$`Down_Regulated_Stressed_Oligos - AUCell score`
Idents(seurat.obj3)<-seurat.obj3@meta.data$cell.type...subgroup..standardized.
seurat.obj3<-seurat.obj3<-subset(seurat.obj3, ident = "oligodendrocyte")
mean(seurat.obj3@meta.data$Up_regulated_Enrichment[seurat.obj3@meta.data$Condition == "Normal"])
mean(seurat.obj3@meta.data$Up_regulated_Enrichment[seurat.obj3@meta.data$Condition == "Multiple sclerosis (MS)"])
mean(seurat.obj3@meta.data$Down_regulated_Enrichment[seurat.obj3@meta.data$Condition == "Normal"])
mean(seurat.obj3@meta.data$Down_regulated_Enrichment[seurat.obj3@meta.data$Condition == "Multiple sclerosis (MS)"])
seurat.obj3 <- FindVariableFeatures(object = seurat.obj3)
seurat.obj3 <- ScaleData(object = seurat.obj3)
seurat.obj3 <- RunPCA(object = seurat.obj3)
seurat.obj3 <- RunUMAP(object = seurat.obj3, dims = 1:15)
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj3@meta.data$Up_regulated_Enrichment)))
FeaturePlot(seurat.obj3, cols= myColor, reduction = "umap", features = "Up_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj3@meta.data$Down_regulated_Enrichment)))
FeaturePlot(seurat.obj3, cols= myColor, reduction = "umap", features = "Down_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj3@meta.data$Condition)))
DimPlot(seurat.obj3, cols = myColor, reduction = 'umap', group.by = "Condition", raster = FALSE)

seurat.obj4<-readRDS("mouse_datasets/EAE_Internal_Hammond/seurat_object/CSF1r inhibitor mouse EAE.rds")
obj4meta<-as.data.frame(read_tsv("OligoAnalysis/seurat_object/EAE metadata (AUCellEnrichment).tsv"))
seurat.obj4@meta.data$Up_regulated_Enrichment<-obj4meta$`AUCell scores - Up_Regulated_Stressed_Oligos`
seurat.obj4@meta.data$Down_regulated_Enrichment<-obj4meta$`AUCell scores - Down_Regulated_Stressed_Oligos`
Idents(seurat.obj4)<-seurat.obj4@meta.data$cell_type
seurat.obj4<-subset(seurat.obj4, ident = "Oligodendrocyte")
seurat.obj4@meta.data$condition[seurat.obj4@meta.data$condition == "Naive"]<-"control"
seurat.obj4@meta.data$condition[seurat.obj4@meta.data$condition == "EAE+vehicle"]<-"MS"
seurat.obj4@meta.data$condition[seurat.obj4@meta.data$condition == "EAE+CSF1Ri"]<-"treated"
seurat.obj4@meta.data$keep<-"keep"
seurat.obj4@meta.data$keep[seurat.obj4@meta.data$condition == "treated"] <-"toss"
Idents(seurat.obj4)<-seurat.obj4@meta.data$keep
seurat.obj4<-subset(seurat.obj4, ident = "keep")
mean(seurat.obj4@meta.data$Up_regulated_Enrichment[seurat.obj4@meta.data$condition == "control"])
mean(seurat.obj4@meta.data$Up_regulated_Enrichment[seurat.obj4@meta.data$condition == "MS"])
mean(seurat.obj4@meta.data$Down_regulated_Enrichment[seurat.obj4@meta.data$condition == "control"])
mean(seurat.obj4@meta.data$Down_regulated_Enrichment[seurat.obj4@meta.data$condition == "MS"])
seurat.obj4 <- FindVariableFeatures(object = seurat.obj4)
seurat.obj4 <- ScaleData(object = seurat.obj4)
seurat.obj4 <- RunPCA(object = seurat.obj4)
seurat.obj4 <- RunUMAP(object = seurat.obj4, dims = 1:15)
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj4@meta.data$Up_regulated_Enrichment)))
FeaturePlot(seurat.obj4, cols= myColor, reduction = "umap", features = "Up_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj4@meta.data$Down_regulated_Enrichment)))
FeaturePlot(seurat.obj4, cols= myColor, reduction = "umap", features = "Down_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj4@meta.data$condition)))
DimPlot(seurat.obj4, cols = myColor, reduction = 'umap', group.by = "condition", raster = FALSE)

seurat.obj5<-readRDS("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/PMCB/SCB/hOligo_MS_data/EBPi_snRNAseq_full_seurat.RDS")
obj5meta<-as.data.frame(read_tsv("OligoAnalysis/seurat_object/Cuprizone metadata (AUCellEnrichment).tsv"))
seurat.obj5@meta.data$Up_regulated_Enrichment<-obj5meta$`AUCell scores - Up_Regulated_Stressed_Oligos`
seurat.obj5@meta.data$Down_regulated_Enrichment<-obj5meta$`AUCell scores - Down_Regulated_Stressed_Oligos`
Idents(seurat.obj5)<-seurat.obj5@meta.data$sargent_cellstates
seurat.obj5<-subset(seurat.obj5, ident = "OLIGODENDROCYTES")
seurat.obj5@meta.data$condition<-"control"
seurat.obj5@meta.data$condition[seurat.obj5@meta.data$cuprizone == "Cuprizone"]<-"MS"
seurat.obj5@meta.data$condition[seurat.obj5@meta.data$treatment == "RA742"]<-"treated"
seurat.obj5@meta.data$keep<-"keep"
seurat.obj5@meta.data$keep[seurat.obj5@meta.data$condition == "treated"] <-"toss"
Idents(seurat.obj5)<-seurat.obj5@meta.data$keep
seurat.obj5<-subset(seurat.obj5, ident = "keep")
mean(seurat.obj5@meta.data$Up_regulated_Enrichment[seurat.obj5@meta.data$condition == "control"])
mean(seurat.obj5@meta.data$Up_regulated_Enrichment[seurat.obj5@meta.data$condition == "MS"])
mean(seurat.obj5@meta.data$Down_regulated_Enrichment[seurat.obj5@meta.data$condition == "control"])
mean(seurat.obj5@meta.data$Down_regulated_Enrichment[seurat.obj5@meta.data$condition == "MS"])
seurat.obj5 <- FindVariableFeatures(object = seurat.obj5)
seurat.obj5 <- ScaleData(object = seurat.obj5)
seurat.obj5 <- RunPCA(object = seurat.obj5)
seurat.obj5 <- RunUMAP(object = seurat.obj5, dims = 1:15)
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj5@meta.data$Up_regulated_Enrichment)))
FeaturePlot(seurat.obj5, cols= myColor, reduction = "umap", features = "Up_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj5@meta.data$Down_regulated_Enrichment)))
FeaturePlot(seurat.obj5, cols= myColor, reduction = "umap", features = "Down_regulated_Enrichment")
myColor <- colorRampPalette(c("blue", "white", "red"))(length(unique(seurat.obj5@meta.data$condition)))
DimPlot(seurat.obj5, cols = myColor, reduction = 'umap', group.by = "condition", raster = FALSE)