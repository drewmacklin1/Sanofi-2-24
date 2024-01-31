#GEO to seurat
#https://www.youtube.com/watch?v=5e_8wr5Xx_Q
#https://www.youtube.com/watch?v=43Z13DS_emQ

cd cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/DMD_datasets/temp/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213925/suppl/GSE213925_RAW.tar
tar xvf GSE213925_RAW.tar
rm  GSE213925_RAW.tar
y

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(tidyverse)
library(Seurat)

sample_names<- list.files("DMD_datasets/temp", pattern = "*matrix.mtx.gz")
sample_names<- str_replace(sample_names, "(GSM[0-9]+_.+)_matrix.mtx.gz", "\\1")
for (x in 1:length(sample_names)){
  dir.create(paste0("DMD_datasets/temp/", sample_names[x]))
  for(y in 1:length(list.files("DMD_datasets/temp"))){
    if(grepl(sample_names[x], list.files("DMD_datasets/temp")[y])& grepl(".gz", list.files("DMD_datasets/temp")[y])){
      file.rename(from = paste0("DMD_datasets/temp/", list.files("DMD_datasets/temp")[y]),  to = paste0("DMD_datasets/temp/", sample_names[x], "/",list.files("DMD_datasets/temp")[y]))
    }
  }
}

sample_names<-list.files("DMD_datasets/temp")[c(1,3,5,7,9,11)]
for (x in 1:length(sample_names)){
  dir.create(paste0("DMD_datasets/temp/", sample_names[x]))
  for(y in 1:length(list.files("DMD_datasets/temp"))){
    if(grepl(sample_names[x], list.files("DMD_datasets/temp")[y])& grepl(".gz", list.files("DMD_datasets/temp")[y])){
      file.rename(from = paste0("DMD_datasets/temp/", list.files("DMD_datasets/temp")[y]),  to = paste0("DMD_datasets/temp/", sample_names[x], "/",list.files("DMD_datasets/temp")[y]))
    }
  }
}

sample_names<-list.files("DMD_datasets/temp")
for(x in 1:length(sample_names)){
  for(y in 1:3){
    if(grepl("barcodes", paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]))){
      file.rename(from = paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]),  to = paste0("DMD_datasets/temp/", sample_names[x], "/barcodes.tsv.gz"))
    }
    if(grepl("features", paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]))){
      file.rename(from = paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]),  to = paste0("DMD_datasets/temp/", sample_names[x], "/features.tsv.gz"))
    }
    if(grepl("matrix", paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]))){
      file.rename(from = paste0("DMD_datasets/temp/", sample_names[x], "/", list.files(paste0("DMD_datasets/temp/", sample_names[x]))[y]),  to = paste0("DMD_datasets/temp/", sample_names[x], "/matrix.mtx.gz"))
    }
  }
}

files<-list.files("DMD_datasets/temp")
print(length(files))
GSM6596509_wtNSG01<-Read10X(paste0("DMD_datasets/temp/", files[1]))
GSM6596510_wtNSG02<-Read10X(paste0("DMD_datasets/temp/", files[2]))
GSM6596511_mdxNSG01<-Read10X(paste0("DMD_datasets/temp/", files[3]))
GSM6596512_mdxNSG02<-Read10X(paste0("DMD_datasets/temp/", files[4]))
GSM6596513_mdxD2NSG01<-Read10X(paste0("DMD_datasets/temp/", files[5]))
GSM6596514_mdxD2NSG02<-Read10X(paste0("DMD_datasets/temp/", files[6]))

GSM6596509_wtNSG01<-CreateSeuratObject(counts = GSM6596509_wtNSG01, project = "GSM6596509_wtNSG01", min.cells = 3, min.features = 200)
GSM6596510_wtNSG02<-CreateSeuratObject(counts = GSM6596510_wtNSG02, project = "GSM6596510_wtNSG02", min.cells = 3, min.features = 200)
GSM6596511_mdxNSG01<-CreateSeuratObject(counts = GSM6596511_mdxNSG01, project = "GSM6596511_mdxNSG01", min.cells = 3, min.features = 200)
GSM6596512_mdxNSG02<-CreateSeuratObject(counts = GSM6596512_mdxNSG02, project = "GSM6596512_mdxNSG02", min.cells = 3, min.features = 200)
GSM6596513_mdxD2NSG01<-CreateSeuratObject(counts = GSM6596513_mdxD2NSG01, project = "GSM6596513_mdxD2NSG01", min.cells = 3, min.features = 200)
GSM6596514_mdxD2NSG02<-CreateSeuratObject(counts = GSM6596514_mdxD2NSG02, project = "GSM6596514_mdxD2NSG02", min.cells = 3, min.features = 200)

merged<-merge(GSM6596509_wtNSG01, c(GSM6596510_wtNSG02,GSM6596511_mdxNSG01,GSM6596512_mdxNSG02,GSM6596513_mdxD2NSG01,GSM6596514_mdxD2NSG02))

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

p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = "orig.ident", raster = FALSE)
p1

merged_seurat_filtered@meta.data$condition[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596509_wtNSG01", "GSM6596510_wtNSG02")]<-"control"
merged_seurat_filtered@meta.data$condition[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596511_mdxNSG01", "GSM6596512_mdxNSG02")]<-"dystrophic"
merged_seurat_filtered@meta.data$condition[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596513_mdxD2NSG01", "GSM6596514_mdxD2NSG02")]<-"severely dystrophic"
merged_seurat_filtered@meta.data$genotype[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596509_wtNSG01", "GSM6596510_wtNSG02")]<-"Scid , Il2rgnull"
merged_seurat_filtered@meta.data$genotype[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596511_mdxNSG01", "GSM6596512_mdxNSG02")]<-"Mdx, Scid, Il2rgnull"
merged_seurat_filtered@meta.data$genotype[merged_seurat_filtered@meta.data$orig.ident %in% c("GSM6596513_mdxD2NSG01", "GSM6596514_mdxD2NSG02")]<-"Mdx, Scid, Il2rgnull, Anxa6, Ltbp4"
merged_seurat_filtered@meta.data$age<-"8 weeks"
merged_seurat_filtered@meta.data$tissue<-"Right Hindlimb Gastrocnemius"
merged_seurat_filtered@meta.data$treatment<-"not treated"
View(merged_seurat_filtered@meta.data)

saveRDS(merged_seurat_filtered, file = "DMD_datasets/temp/seurat.rds")

#upload to BioTuring for celltype ID if no batch correction needed