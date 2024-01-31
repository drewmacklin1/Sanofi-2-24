setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
library(Seurat)

studylist<-c('MS_Cortex_Internal_Proto', 'PD_Midbrain_PMID34919646', 'AD_EC_SFG_PMID33432193', 'PD_PreFCortex_PPR453641', 'FTD_FrontTempCortex_PMID35879464', 'AD_OC_OTC_PMID33609158', 'ALS_FTLD_PrimMotorCortex_GSE174332', 'PD_SubNigra_PMID35513515', 'ALS_Spinalcord_Internal_Zelic') 
study <- studylist[6]
seurat.obj_name <- list.files(paste0("datasets/", study, "/seurat_object"))
seurat.obj <- readRDS(paste0('datasets/', study, "/seurat_object/", seurat.obj_name))

metadata <- as.data.frame(seurat.obj@meta.data)
View(metadata)

seurat.obj@meta.data$cell_type <- seurat.obj@meta.data$Sub.cell.type..version.3.
seurat.obj@meta.data$sample_id <- seurat.obj@meta.data$orig.ident
seurat.obj@meta.data$condition <- NA
seurat.obj@meta.data$gender <- NA
seurat.obj@meta.data$condition[seurat.obj@meta.data$sample_id %in% c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p' )] <- 'control'
seurat.obj@meta.data$condition[seurat.obj@meta.data$sample_id %in% c('q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'aa', 'ab', 'ac', 'ad', 'ae', 'af','ag', 'ah','ai', 'aj')] <- 'AD'
seurat.obj@meta.data$gender[seurat.obj@meta.data$sample_id %in% c('a', 'b', 'c', 'd', 'e', 'f', 'k', 'l', 'o', 'p', 'u', 'v', 'ai', 'aj')] <- 'male'
seurat.obj@meta.data$gender[seurat.obj@meta.data$sample_id %in% c('g', 'h', 'i', 'j', 'm', 'n', 'q', 'r', 's', 't', 'w', 'x', 'y', 'z', 'aa', 'ab','ac', 'ad','ae', 'af', 'ag', 'ah')] <- 'female'
seurat.obj@meta.data$sample_id <- paste(seurat.obj@meta.data$sample_id, study)
seurat.obj@meta.data$cell_type[seurat.obj@meta.data$cell_type == "microglial cell"] <- "Microglia"

print(unique(seurat.obj@meta.data$cell_type))
print(unique(seurat.obj@meta.data$condition))
print(unique(seurat.obj@meta.data$gender))
print(unique(seurat.obj@meta.data$sample_id))

Idents(seurat.obj)<-seurat.obj@meta.data$cell_type
seurat.obj<-subset(seurat.obj, idents = "Microglia")
saveRDS(seurat.obj, file = paste0("microglia_markers/Atlas/", study, ".rds"))

files<-list.files('microglia_markers/Atlas')

a <- readRDS(paste0('microglia_markers/Atlas/', files[1]))
a@meta.data$study<-files[1]
b <- readRDS(paste0('microglia_markers/Atlas/', files[2]))
b@meta.data$study<-files[2]
c <- readRDS(paste0('microglia_markers/Atlas/', files[3]))
c@meta.data$study<-files[3]
d <- readRDS(paste0('microglia_markers/Atlas/', files[4]))
d@meta.data$study<-files[4]
e <- readRDS(paste0('microglia_markers/Atlas/', files[5]))
e@meta.data$study<-files[5]
f <- readRDS(paste0('microglia_markers/Atlas/', files[7]))
f@meta.data$study<-files[7]
g <- readRDS(paste0('microglia_markers/Atlas/', files[8]))
g@meta.data$study<-files[8]
h <- readRDS(paste0('microglia_markers/Atlas/', files[9]))
h@meta.data$study<-files[9]
i <- readRDS(paste0('microglia_markers/Atlas/', files[10]))
i@meta.data$study<-files[10]
j <- readRDS(paste0('microglia_markers/Atlas/', files[11]))
j@meta.data$study<-files[11]
k <- readRDS(paste0('microglia_markers/Atlas/', files[12]))
k@meta.data$study<-files[12]

merged1<-merge(a,b)
merged1<-merge(merged1, c)
merged1<-merge(merged1, d)
merged1<-merge(merged1, e)
merged1<-merge(merged1, f)
merged1<-merge(merged1, g)
merged1<-merge(merged1, h)
merged1<-merge(merged1, i)
merged1<-merge(merged1, j)
merged1<-merge(merged1, k)

merged1@meta.data$condition[merged1@meta.data$condition=="Control"] <- 'control'
merged1@meta.data$condition[merged1@meta.data$condition=="PD"] <- "Parkinson's disease"
merged1@meta.data$gender[merged1@meta.data$gender=="M"] <- 'Male'
merged1@meta.data$gender[merged1@meta.data$gender=="male"] <- 'Male'
merged1@meta.data$gender[merged1@meta.data$gender=="F"] <- 'Female'
merged1@meta.data$gender[merged1@meta.data$gender=="female"] <- 'Female'

print(unique(merged1@meta.data$cell_type))
print(unique(merged1@meta.data$condition))
print(unique(merged1@meta.data$gender))
print(unique(merged1@meta.data$sample_id))
print(unique(merged1@meta.data$study))

merged1<- NormalizeData(object = merged1)
merged1 <- FindVariableFeatures(object = merged1)
merged1 <- ScaleData(object = merged1)
merged1<- RunPCA(object = merged1)
merged1 <- FindNeighbors(object = merged1, dims = 1:15)
merged1 <- FindClusters(object = merged1, resolution = 0.1)
merged1 <- RunUMAP(object = merged1, dims = 1:15)

p1 <- DimPlot(merged1, reduction = 'umap', group.by = "study", raster = FALSE)
p1

saveRDS(merged1, "microglia_markers/Atlas/merged1.rds")

#############################################################################
merged1<-readRDS("microglia_markers/Atlas/merged1.rds")

obj.list <- SplitObject(merged1, split.by = 'study')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features) #takes abut 24 hours and 256gb ram
merged1 <- IntegrateData(anchorset = anchors)
merged1 <- ScaleData(object = merged1)
merged1 <- RunPCA(object = merged1)
merged1 <- RunUMAP(object = merged1, dims = 1:15)
saveRDS(merged1, file = "microglia_markers/Atlas/merged1-int.rds")
##########################################################################################
merged1<-readRDS("microglia_markers/Atlas/merged1-int.rds")

p2 <- DimPlot(merged1, reduction = 'umap', group.by = 'study', raster = FALSE)
p3 <- DimPlot(merged1, reduction = 'umap', group.by = 'condition', raster = FALSE)

library(gridExtra)
grid.arrange(p2, p3, ncol = 2)