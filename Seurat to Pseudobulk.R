#convert seurat objects into pseudobulk counts and meta data files for every celltype present
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
library(tidyverse)
library(Seurat)

study<-list.files("datasets")[1]
seurat.obj_name <- list.files(paste0("datasets/", study, "/seurat_object"))
seurat.obj <- readRDS(paste0('datasets/', study, "/seurat_object/", seurat.obj_name))

#standardize metadata, Define cell_type, sample_id, condition, and gender column names and define control (do this manually)
metadata2 <- as.data.frame(seurat.obj@meta.data)
View(metadata2)

metadata2$cell_type <- seurat.obj@meta.data$
  metadata2$sample_id <- seurat.obj@meta.data$
  metadata2$gender <- seurat.obj@meta.data$
  metadata2$condition <- seurat.obj@meta.data$
  metadata2$cell_type[metadata2$cell_type == ''] <- ""
metadata2$condition[metadata2$condition == ''] <- "control"

print(unique(metadata2$cell_type))
print(unique(metadata2$condition))
print(unique(metadata2$gender))
print(unique(metadata2$sample_id))

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
    write.csv(output_counts, paste0('datasets/', study, "/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_data.csv"))
    write.csv(output_metadata, paste0('datasets/', study, "/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_metadata.csv"))
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
write.csv(output_counts, paste0('datasets/', study, "/pseudobulk_counts_matrices/", "All_Cells_pseudobulk_data.csv"))
write.csv(output_metadata, paste0('datasets/', study, "/pseudobulk_counts_matrices/", "All_Cells_pseudobulk_metadata.csv"))