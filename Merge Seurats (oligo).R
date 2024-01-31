setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
library(Seurat)

files<-list.files("OligoAnalysis/seurat_object")

a<-readRDS(paste0("OligoAnalysis/seurat_object/", files[1]))
unique(a@meta.data$cell.type...subgroup..standardized.)
unique(a@meta.data$cell.type..standardized.)
unique(a@meta.data$Author.s.cell.type)
unique(a@meta.data$detailed.cell.type..standardized.)
unique(a@meta.data$detailed.cell.type...subgroup..standardized.)
Idents(a)<-a@meta.data$cell.type...subgroup..standardized.
a<-subset(a, idents = "oligodendrocyte")
a@meta.data$cell_type<-a@meta.data$cell.type...subgroup..standardized.
a@meta.data$gender<-a@meta.data$gender..standardized.
a@meta.data$condition<-a@meta.data$Condition
a@meta.data$sample_id<-paste(gsub(".rds", "", files[1]), a@meta.data$Sample.ID)
print(unique(a@meta.data$cell_type))
print(unique(a@meta.data$gender))
print(unique(a@meta.data$condition))
print(unique(a@meta.data$sample_id))
saveRDS(a, file = paste0("OligoAnalysis/seurat_object/Oligo", files[1]))


b<-readRDS(paste0("OligoAnalysis/seurat_object/", files[2]))
unique(b@meta.data$cell.type...subgroup..standardized.)
unique(b@meta.data$cell.type..standardized.)
unique(b@meta.data$Author.s.cell.type)
unique(b@meta.data$detailed.cell.type..standardized.)
unique(b@meta.data$detailed.cell.type...subgroup..standardized.)
Idents(b)<-b@meta.data$cell.type...subgroup..standardized.
b<-subset(b, idents = "oligodendrocyte")
b@meta.data$cell_type<-b@meta.data$cell.type...subgroup..standardized.
b@meta.data$gender<-b@meta.data$gender..standardized.
b@meta.data$condition<-b@meta.data$Condition
b@meta.data$sample_id<-paste(gsub(".rds", "", files[2]), b@meta.data$Sample.ID)
print(unique(b@meta.data$cell_type))
print(unique(b@meta.data$gender))
print(unique(b@meta.data$condition))
print(unique(b@meta.data$sample_id))
saveRDS(b, file = paste0("OligoAnalysis/seurat_object/Oligo", files[2]))

c<-readRDS(paste0("OligoAnalysis/seurat_object/", files[9]))
unique(c@meta.data$cell.type...subgroup..standardized.)
unique(c@meta.data$cell.type..standardized.)
unique(c@meta.data$Author.s.cell.type)
unique(c@meta.data$detailed.cell.type..standardized.)
unique(c@meta.data$detailed.cell.type...subgroup..standardized.)
Idents(c)<-c@meta.data$cell.type...subgroup..standardized.
c<-subset(c, idents = "oligodendrocyte")
c@meta.data$cell_type<-c@meta.data$cell.type...subgroup..standardized.
c@meta.data$gender<-c@meta.data$gender..standardized.
c@meta.data$condition<-c@meta.data$Condition
c@meta.data$sample_id<-paste(gsub(".rds", "", files[3]), c@meta.data$Sample.ID)
print(unique(c@meta.data$cell_type))
print(unique(c@meta.data$gender))
print(unique(c@meta.data$condition))
print(unique(c@meta.data$sample_id))
saveRDS(c, file = paste0("OligoAnalysis/seurat_object/Oligo", files[3]))

library(ggplot2)
library(tidyverse)
library(gridExtra)

merged<-merge(a, c(b,c))

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged_seurat_filtered <- subset(merged, subset = nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 5)
rm(list=setdiff(ls(), "merged_seurat_filtered"))

merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:15)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.1)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:15)

p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = "sample_id", raster = FALSE)
p1

saveRDS(merged_seurat_filtered, file = "OligoAnalysis/seurat_object/Oligo_merged_seurat_filtered.rds")
merged_seurat_filtered <- readRDS(file = "OligoAnalysis/seurat_object/Oligo_merged_seurat_filtered.rds")

#do the following if you need to do batch correction
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'sample_id')
obj.list<-obj.list[c(1:12,16:22,24:30, 32:39, 41, 43:46, 48:59)]#objects 13, 14, 15, 23, 31, 40, 42, 47 contain fewer than 30 cells
################     EDIT DIDNT SAVE. RE-ADD HERE    ########### #objects _________________ contain fewer than 100 cells
for(i in 1:length(obj.list)){ 
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])f
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, k.weight = 100) #many samples contain fewer than 100 anchors "Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length" (fewest anchors is 37)
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15)
saveRDS(seurat.integrated, file = "OligoAnalysis/seurat_object/Oligo_integrated(1).rds")
seurat.integrated<-readRDS("OligoAnalysis/seurat_object/Oligo/integrated(1).rds")

p2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'sample_id', raster = FALSE)
p3 <- DimPlot(seurat.integrated, reduction = 'umap', raster = FALSE)
grid.arrange(p2, p3, ncol = 2)

#add dataset column and quantify
seurat.obj<-readRDS("OligoAnalysis/seurat_object/Oligo_integrated(1).rds")
unique(seurat.obj@meta.data$cell.type...subgroup..standardized.)

seurat.obj@meta.data$dataset<-seurat.obj@meta.data$sample_id
x<-NULL
y<-NULL
z<-NULL
for(w in 1:dim(seurat.obj@meta.data)[1]){
  x<-c(x,grepl("GSE118257", seurat.obj@meta.data$dataset[w]))
  y<-c(y,grepl("GSE180759", seurat.obj@meta.data$dataset[w]))
  z<-c(z,grepl("PMID31316211", seurat.obj@meta.data$dataset[w]))
}
seurat.obj@meta.data$dataset[x]<-"GSE118257"
seurat.obj@meta.data$dataset[y]<-"GSE180759"
seurat.obj@meta.data$dataset[z]<-"PMID31316211"
seurat.obj@meta.data$condition[seurat.obj@meta.data$condition=="Normal"]<-"control"
seurat.obj@meta.data$condition[seurat.obj@meta.data$condition=="Normal control white matter"]<-"control"
seurat.obj@meta.data$condition[seurat.obj@meta.data$condition=="healthy control"]<-"control"
seurat.obj@meta.data$condition[seurat.obj@meta.data$condition!="control"]<-"MS"

datasets<-data.frame(matrix(0, nrow =5, ncol =1))
colnames(datasets)<-"cells"
rownames(datasets)<-c("GSE118257", 'GSE180759', 'PMID31316211', "healthy", "MS")
datasets[1,1] <-sum(x==TRUE)
datasets[2,1] <-sum(y==TRUE)
datasets[3,1] <-sum(z==TRUE)
datasets[4,1] <-sum(seurat.obj@meta.data$condition == "control")
datasets[5,1] <-sum(seurat.obj@meta.data$condition == "MS")
View(datasets)
print(paste("GSE118257", sum(x==TRUE)))
print(paste("GSE180759", sum(y==TRUE)))
print(paste("PMID31316211", sum(z==TRUE)))

DimPlot(seurat.obj, reduction = 'umap', group.by = "dataset", raster = FALSE)
DimPlot(seurat.obj, reduction = 'umap', group.by = "condition", raster = FALSE)

saveRDS(seurat.obj, file = "OligoAnalysis/seurat_object/Oligo_integrated(1).rds")
