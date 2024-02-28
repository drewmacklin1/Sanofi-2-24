setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage")

library(Seurat)
library(tidyverse)

seurat.obj1 <- readRDS("PMCB/SCB/hOligo_MS_data/MS_merged_Jakel_Absinthe_Schimmer_Proto_oligo_annotated_seurat.RDS")
p1<-DimPlot(seurat.obj1, reduction = 'umap', group.by = 'orig.ident', raster = FALSE)
p2<-DimPlot(seurat.obj1, reduction = 'umap', raster = FALSE)

#remove opcs
Idents(seurat.obj1)<-seurat.obj1@meta.data$sargent_cellstates
seurat.obj<-subset(seurat.obj1, idents = "OLIGODENDROCYTE")
p3<-DimPlot(seurat.obj, reduction = 'umap', group.by= "seurat_clusters", raster = FALSE)

#compare bottom and top blob #1 *Do this 5 different ways for Anna...
#umap_coords<- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
#bottom<-rownames(umap_coords[(umap_coords$umap_1+1)^2 + (umap_coords$umap_2+6)^2 < 15 ,])
#top<-rownames(umap_coords[(umap_coords$umap_1+1)^2 + (umap_coords$umap_2-6)^2 < 15,])
#bottom<-rownames(umap_coords[umap_coords$umap_2 < -0,])
#top<-rownames(umap_coords[umap_coords$umap_2 > 0,])

seurat.obj@meta.data$top_bottom<-"ignore"
seurat.obj@meta.data$top_bottom[seurat.obj@meta.data$seurat_clusters %in% c(2,3,4,5,8)]<-"top"
seurat.obj@meta.data$top_bottom[seurat.obj@meta.data$seurat_clusters %in% c(0,1)]<-"bottom"

#seurat.obj@meta.data$top_bottom[rownames(seurat.obj@meta.data) %in% top] <-"top"
#seurat.obj@meta.data$top_bottom[rownames(seurat.obj@meta.data) %in% bottom] <-"bottom"
Idents(seurat.obj)<-seurat.obj@meta.data$top_bottom
p4<-DimPlot(seurat.obj, reduction = 'umap', raster = FALSE)

seurat.obj@meta.data$top_bottom_condition<-paste(seurat.obj@meta.data$top_bottom, seurat.obj@meta.data$disease)
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="ignore MS"] <- "ignore"
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="ignore Control"] <- "ignore"
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="top Control"] <- "Control"
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="bottom Control"] <- "Control"
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="top MS"] <- "MS"
seurat.obj@meta.data$top_bottom_condition[seurat.obj@meta.data$top_bottom_condition=="bottom MS"] <- "MS"
unique(seurat.obj@meta.data$top_bottom_condition)
Idents(seurat.obj)<-seurat.obj@meta.data$top_bottom_condition

p5<-DimPlot(seurat.obj, reduction = 'umap', group.by = 'top_bottom_condition', raster = FALSE)
Idents(seurat.obj)<-seurat.obj@meta.data$top_bottom

#MS to healthy ratio
top<-subset(seurat.obj, idents = "top")
bottom<-subset(seurat.obj, idents = "bottom")
TopMS<-sum(top@meta.data$disease == "MS")
TopControl<-sum(top@meta.data$disease == "Control")
BottomMS<-sum(bottom@meta.data$disease == "MS")
BottomControl<-sum(bottom@meta.data$disease == "Control")
TopMS/TopControl
BottomMS/BottomControl


#pseudobulking
#seurat.obj[["SCT"]]@data ##normalized data
#seurat.obj[["SCT"]]@scale.data ## supposedly stabilized data but takes too long to run
#compare disease vs healthy in top vs bottom
metadata2 <- as.data.frame(seurat.obj@meta.data)
metadata2$cell_type <- seurat.obj@meta.data$top_bottom
metadata2$sample_id <- seurat.obj@meta.data$sample_id
#metadata2$gender <- seurat.obj@meta.data$sex
metadata2$condition <- seurat.obj@meta.data$disease
metadata2$condition[metadata2$condition == 'Control'] <- "control"

print(unique(metadata2$cell_type))
print(unique(metadata2$condition))
#print(unique(metadata2$gender))
print(unique(metadata2$sample_id))

celltypes_list<-unique(metadata2$cell_type)

#split by cell type
for(x in 1:length(celltypes_list)){
  celltype_rownames <- rownames(subset(metadata2, metadata2$cell_type==celltypes_list[x]))
  if(length(celltype_rownames)>200){
    celltype_counts <- seurat.obj[["SCT"]]@counts[,colnames(seurat.obj[["SCT"]]@counts) %in% celltype_rownames]
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
    write.csv(output_counts, paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/', celltypes_list[x], "_pseudobulk_data.csv"))
    write.csv(output_metadata, paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/', celltypes_list[x], "_pseudobulk_metadata.csv"))
  }
}
###########################################################################################################
#directly compare top vs bottom regardless of disease state
metadata2 <- as.data.frame(seurat.obj@meta.data)
metadata2$sample_id <- seurat.obj@meta.data$sample_id
metadata2$condition <- seurat.obj@meta.data$top_bottom

print(unique(metadata2$condition))
print(unique(metadata2$sample_id))

sample_ids <- unique(metadata2$sample_id)

#split by sample
output_counts <- data.frame()
output_metadata <- data.frame()
for(y in 1:length(sample_ids)){
  sample_ids_rownames <- rownames(subset(metadata2, metadata2$sample_id==sample_ids[y]))
  grouped_sample_counts <- as.matrix(seurat.obj[["SCT"]]@counts[,colnames(seurat.obj[["SCT"]]@counts) %in% sample_ids_rownames])
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
write.csv(output_counts, paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/All_Cells_pseudobulk_data.csv'))
write.csv(output_metadata, paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/All_Cells_pseudobulk_metadata.csv'))

#Diffex
library(DESeq2)
files<-list.files('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices')
study<-"OligoAnalysis"
for (j in c(1)){
  data_counts<-as.data.frame(read.csv(paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/', files[j])))
  rownames(data_counts)<-data_counts$X
  data_counts<-data_counts[-1]
  data_counts <- round(data_counts)
  metadata<-as.data.frame(read.csv(paste0('Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/', files[j+1])))
  rownames(metadata) <- colnames(data_counts)
  metadata$sample <- rownames(metadata)
  
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
  
  result <- results(dds, contrast = c('condition', 'top', 'bottom'))
  
  #Convert Results to dataframe and add celltype identifier and comparison identifier
  final_results <- data.frame(result)
  final_results$celltype <- gsub(".csv", "", gsub("_pseudobulk_data", "", files[j]))
  final_results$gene_name <- rownames(final_results)
  #final_results$comparison <- gsub(" ", "_", paste(disease, "vs.",'Control'))
  final_results$comparison <- "Top vs Bottom"
  final_results$dataset.origin <- study
  
  #Export results as csv to folder
  write.csv(final_results, paste0("Drew_Macklin/OligoAnalysis/differential_expression_results/diffex_results_", final_results$celltype[1], final_results$comparison[1], ".csv"))
}

#DEG List
significant <- subset(filter(final_results, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
significant <- subset(filter(significant, baseMean>150))
write.csv(significant, "Drew_Macklin/OligoAnalysis/differential_expression_results/DEG_list.csv")

#top 10 DEGs by L2FC
significant$abs<-abs(significant$log2FoldChange)
significant<-significant[order(-significant$abs), ]
top10<-rownames(significant)[1:10]

#Feature plot
FeaturePlot(object = seurat.obj, features = top10)

#Normalized counts
counts <- as.data.frame(counts(dds, normalized = TRUE))
counts$gene_name <- rownames(counts)
write.csv(counts, "Drew_Macklin/OligoAnalysis/pseudobulk_counts_matrices/All_Cells_pseudobulk_metadata(normalized).csv")
counts$mean<-rowMeans(counts[1:21],na.rm = TRUE)
morethan500<-rownames(counts[counts$mean>500,])

#Split plot
DimPlot(seurat.obj, reduction = 'umap', split.by= "disease", raster = FALSE)

#Cluster contributions of each dataset
Idents(seurat.obj)<-seurat.obj@meta.data$seurat_clusters
seurat.obj@meta.data$seurat_clusters<-as.character(seurat.obj@meta.data$seurat_clusters)
seurat.obj@meta.data$project_name<-as.character(seurat.obj@meta.data$project_name)
dataset_fraction<-data.frame(matrix(0, nrow =length(unique(seurat.obj@meta.data$seurat_clusters)) , ncol = length(unique(seurat.obj@meta.data$project_name))))
colnames(dataset_fraction)<-unique(seurat.obj@meta.data$project_name)
rownames(dataset_fraction)<-unique(seurat.obj@meta.data$seurat_clusters)
for(i in 1:dim(dataset_fraction)[1]){ #for every seurat cluster
  for(j in 1:length(unique(seurat.obj@meta.data$project_name))){ #for every dataset
    dataset_fraction[unique(seurat.obj@meta.data$seurat_clusters)[i], unique(seurat.obj@meta.data$project_name)[j]]<-sum(seurat.obj@meta.data$seurat_clusters %in% unique(seurat.obj@meta.data$seurat_clusters)[i] & seurat.obj@meta.data$project_name %in% unique(seurat.obj@meta.data$project_name)[j])
  }
  dataset_fraction$total[i]<-  print(sum(seurat.obj@meta.data$seurat_clusters %in% unique(seurat.obj@meta.data$seurat_clusters)[i]))
}
dataset_fraction[1]<-dataset_fraction[1]/dataset_fraction[5]
dataset_fraction[2]<-dataset_fraction[2]/dataset_fraction[5]
dataset_fraction[3]<-dataset_fraction[3]/dataset_fraction[5]
dataset_fraction[4]<-dataset_fraction[4]/dataset_fraction[5]
write.csv(dataset_fraction, "Drew_Macklin/OligoAnalysis/dataset_cluster_contributions.csv")

#top/bottom contributions of each dataset
seurat.obj@meta.data$top_bottom<-as.character(seurat.obj@meta.data$top_bottom)
seurat.obj@meta.data$project_name<-as.character(seurat.obj@meta.data$project_name)
dataset_fraction2<-data.frame(matrix(0, nrow =length(unique(seurat.obj@meta.data$top_bottom)) , ncol = length(unique(seurat.obj@meta.data$project_name))))
colnames(dataset_fraction2)<-unique(seurat.obj@meta.data$project_name)
rownames(dataset_fraction2)<-unique(seurat.obj@meta.data$top_bottom)
for(i in 1:dim(dataset_fraction2)[1]){ #for every top/bottom
  for(j in 1:length(unique(seurat.obj@meta.data$project_name))){ #for every dataset
    dataset_fraction2[unique(seurat.obj@meta.data$top_bottom)[i], unique(seurat.obj@meta.data$project_name)[j]]<-sum(seurat.obj@meta.data$top_bottom %in% unique(seurat.obj@meta.data$top_bottom)[i] & seurat.obj@meta.data$project_name %in% unique(seurat.obj@meta.data$project_name)[j])
  }
  dataset_fraction2$total[i]<-  print(sum(seurat.obj@meta.data$top_bottom %in% unique(seurat.obj@meta.data$top_bottom)[i]))
}
dataset_fraction2[1]<-dataset_fraction2[1]/dataset_fraction2[5]
dataset_fraction2[2]<-dataset_fraction2[2]/dataset_fraction2[5]
dataset_fraction2[3]<-dataset_fraction2[3]/dataset_fraction2[5]
dataset_fraction2[4]<-dataset_fraction2[4]/dataset_fraction2[5]
write.csv(dataset_fraction2, "Drew_Macklin/OligoAnalysis/dataset_TopBottom_contributions.csv")

#Pseudobulk barplots
library(ggprism)
counts <- as.data.frame(counts(dds, normalized = TRUE))
counts$gene_name <- rownames(counts)
df_counts <- gather(as.data.frame(counts), key = 'sample', value = 'counts', -gene_name)
df_counts <- left_join(df_counts, select(metadata, sample, condition), by = 'sample')
program_gene_list<-top10
for(z in 1:length(program_gene_list)){
  gene_plot <- program_gene_list[z]
  df_counts_plot <- filter(df_counts, gene_name == gene_plot)
  df_counts_plot_summary <- df_counts_plot %>%
    group_by(condition) %>%
    summarise(mean = mean(counts),
              sd = sd(counts))
  diseases <- as.character(unique(metadata[metadata$condition!='control', 'condition']))                                   ###NEW
  df_counts_plot_summary <- df_counts_plot_summary %>% mutate(condition = fct_relevel(condition, "control", diseases))    ###NEW
  bar_plot <- ggplot(df_counts_plot_summary, aes(x = condition, y = mean, fill = condition, colour = condition)) +
    geom_bar(stat = 'identity', width = 0.5, alpha = .5) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
    geom_jitter(data = df_counts_plot, aes(x = condition, y = counts, shape = condition), size = 1, width = 0.25) +
    theme_prism(base_size=10) +
    scale_fill_prism(palette='shades_of_gray')+
    scale_colour_prism(palette='shades_of_gray')+
    scale_shape_prism(palette='filled')+
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 10)) +
    ggtitle(gene_plot) +
    ylab('normalized counts')
  ggsave(filename = paste0(gene_plot, "_barplot.png"),
         plot = bar_plot,
         device = 'png',
         path = paste0("Drew_Macklin/OligoAnalysis/bar_plots"),
         width = 4,
         height = 4,
         units = "in")
}

#recluster new integrated object
seurat.obj<-readRDS("OligoAnalysis/seurat_object/Oligo_integrated(1).rds")
seurat.obj2<-readRDS("OligoAnalysis/seurat_object/Oligo_merged_seurat_filtered.rds")
print(paste(dim(seurat.obj@meta.data),dim(seurat.obj2@meta.data)))

seurat.obj <- FindNeighbors(object = seurat.obj, dims = 1:15)
seurat.obj <- FindClusters(object = seurat.obj, resolution = 0.8)
seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:15)
View(seurat.obj@meta.data)
DimPlot(seurat.obj, reduction = 'umap', group.by = "seurat_clusters", raster = FALSE)
#seurat.obj@meta.data$seurat_clusters<-seurat.obj@meta.data$integrated_snn_res.0.05
for (x in 1:length(unique(seurat.obj@meta.data$seurat_clusters))){
  MSfraction<-sum(seurat.obj@meta.data$seurat_clusters == unique(seurat.obj@meta.data$seurat_clusters)[x]&seurat.obj@meta.data$condition == "MS")/sum(seurat.obj@meta.data$seurat_clusters == unique(seurat.obj@meta.data$seurat_clusters)[x])
  Ctrlfraction<-sum(seurat.obj@meta.data$seurat_clusters == unique(seurat.obj@meta.data$seurat_clusters)[x]&seurat.obj@meta.data$condition == "control")/sum(seurat.obj@meta.data$seurat_clusters == unique(seurat.obj@meta.data$seurat_clusters)[x])
  print(paste("ms fraction", MSfraction, x-1))
  print(paste("control fraction", Ctrlfraction, x-1))
}
FeaturePlot(object = seurat.obj, features = 'KLK6')

#compare diffex of GSE180759 oligos to integrated all_cells of integrated oligo object (oligos)
seurat.obj<-readRDS("OligoAnalysis/seurat_object/OligoGSE180759.rds")
metadata2 <- as.data.frame(seurat.obj@meta.data)
metadata2$sample_id <- seurat.obj@meta.data$sample_id
metadata2$condition <- "MS"
metadata2$condition[metadata2$Condition == "Normal control white matter"] <- "control"

print(unique(metadata2$condition))
print(unique(metadata2$sample_id))

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
write.csv(output_counts, paste0('OligoAnalysis/pseudobulk_counts_matrices/All_Cells(GSE180759)_pseudobulk_data.csv'))
write.csv(output_metadata, paste0('OligoAnalysis/pseudobulk_counts_matrices/All_Cells(GSE180759)_pseudobulk_metadata.csv'))

#Diffex
library(DESeq2)
files<-list.files('OligoAnalysis/pseudobulk_counts_matrices')
study<-"OligoAnalysis"
for (j in c(4)){
  data_counts<-as.data.frame(read.csv(paste0('OligoAnalysis/pseudobulk_counts_matrices/', files[j])))
  rownames(data_counts)<-data_counts$X
  data_counts<-data_counts[-1]
  data_counts <- round(data_counts)
  metadata<-as.data.frame(read.csv(paste0('OligoAnalysis/pseudobulk_counts_matrices/', files[j+1])))
  rownames(metadata) <- colnames(data_counts)
  metadata$sample <- rownames(metadata)
  
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
  
  result <- results(dds, contrast = c('condition', 'control', 'MS'))
  
  #Convert Results to dataframe and add celltype identifier and comparison identifier
  final_results <- data.frame(result)
  final_results$celltype <- gsub(".csv", "", gsub("_pseudobulk_data", "", files[j]))
  final_results$gene_name <- rownames(final_results)
  #final_results$comparison <- gsub(" ", "_", paste(disease, "vs.",'Control'))
  final_results$comparison <- "MS vs Control"
  final_results$dataset.origin <- study
  
  #Export results as csv to folder
  write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_", final_results$celltype[1], final_results$comparison[1], "(GSE180759).csv"))
}

#make metadata chart
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
seurat.obj3<-readRDS("OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj3<-readRDS("OligoAnalysis/seurat_object/GSE118257.rds")
seurat.obj3<-readRDS("OligoAnalysis/seurat_object/PMID31316211.rds")

y<-21
x<-unique(seurat.obj3@meta.data$Sample.ID)[y]
x
unique(seurat.obj3@meta.data$Condition[seurat.obj3@meta.data$Sample.ID == x])
unique(seurat.obj3@meta.data$Tissue[seurat.obj3@meta.data$Sample.ID == x])
unique(seurat.obj3@meta.data$Gender[seurat.obj3@meta.data$Sample.ID == x])
unique(seurat.obj3@meta.data$Age..years.[seurat.obj3@meta.data$Sample.ID == x])
sum(seurat.obj3@meta.data$cell.type...subgroup..standardized.[seurat.obj3@meta.data$Sample.ID == x]=="oligodendrocyte")

#recreate GSE180579's paper's extended data figure 3
seurat.obj<-readRDS("OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj@meta.data$keep<-FALSE
seurat.obj@meta.data$keep[seurat.obj@meta.data$Author.s.cell.type %in% c("Oligodendrocytes", "OPC - oligodendrocyte progenitor cells")]<-TRUE
Idents(seurat.obj)<-seurat.obj@meta.data$keep
seurat.obj<-subset(seurat.obj, idents = TRUE)
seurat.obj@meta.data$cell_subtype<- seurat.obj@meta.data$Author.s.cell.subtype..Oligodendrocyte.
seurat.obj@meta.data$cell_subtype[seurat.obj@meta.data$cell_subtype == "Unassigned"]<- seurat.obj@meta.data$Author.s.cell.type[seurat.obj@meta.data$cell_subtype == "Unassigned"]
#seurat.obj <- FindVariableFeatures(object = seurat.obj)
#seurat.obj <- ScaleData(object = seurat.obj)
#seurat.obj <- RunPCA(object = seurat.obj)
#seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:5)
#DimPlot(seurat.obj, reduction = 'umap', group.by = "cell_subtype", raster = FALSE)
#VlnPlot(seurat.obj,  c("SOX6", "OLIG2","PDGFRA","CNP", "ASPA", "MAG", "PLP1","MOG","HSP90AA1", "LINGO1", "KLK6","OPALIN"), cols = NULL, pt.size = 0, idents = NULL, sort = FALSE, assay = NULL, group.by = "cell_subtype", split.by = NULL, adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE, ncol = NULL, split.plot = FALSE, stack = TRUE, combine = TRUE, fill.by = "feature", flip = FALSE, add.noise = TRUE, raster = NULL)

seurat.obj@meta.data$cell_type<-seurat.obj@meta.data$cell_subtype
seurat.obj@meta.data$keep<-FALSE
seurat.obj@meta.data$keep[seurat.obj@meta.data$cell_subtype %in% c("stressed oligodendrocyte", "normal oligodendrocyte")]<-TRUE
Idents(seurat.obj)<-seurat.obj@meta.data$keep
seurat.obj<-subset(seurat.obj, idents = TRUE)
seurat.obj@meta.data$sample_id<-seurat.obj@meta.data$Sample.ID
seurat.obj@meta.data$condition<-seurat.obj@meta.data$cell_subtype
seurat.obj@meta.data$gender<-seurat.obj@meta.data$Gender
seurat.obj$patient_info<-paste(seurat.obj@meta.data$Age..years., seurat.obj@meta.data$gender..standardized., seurat.obj@meta.data$Clinical.notes)

metadata2 <- as.data.frame(seurat.obj@meta.data)
sample_ids <- unique(metadata2$sample_id)
output_counts <- data.frame()
output_metadata <- data.frame()
for(y in 1:length(sample_ids)){
  sample_ids_rownames <- rownames(subset(metadata2, metadata2$sample_id==sample_ids[y]))
  grouped_sample_counts <- as.matrix(seurat.obj$RNA@counts[,colnames(seurat.obj$RNA@counts) %in% sample_ids_rownames])
  grouped_sample_metadata <- metadata2[rownames(metadata2) %in% sample_ids_rownames,]
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

library(DESeq2)
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
#significant <- subset(filter(final_results, abs(log2FoldChange)>.585))
#significant <- subset(filter(significant, pvalue<.05))
#significant <- subset(filter(significant, baseMean>50))
write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_(stressed.vs.normal oligos gse180759 - patient controlled).csv"))





#Normal vs stressed oligos with tissue split and patient id control
seurat.obj@meta.data$cell_type<-seurat.obj@meta.data$Condition #celltype is to run tissue split
metadata2 <- as.data.frame(seurat.obj@meta.data)
celltypes_list<-unique(metadata2$cell_type)

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
    #significant <- subset(filter(final_results, abs(log2FoldChange)>.585))
    #significant <- subset(filter(significant, pvalue<.05))
    #significant <- subset(filter(significant, baseMean>50))
    write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_(stressed.vs.normal oligos gse180759 - patient controlled)", final_results$tissue[1], ".csv"))
     #write.csv(output_counts, paste0('datasets/', study, "/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_data.csv"))
    #write.csv(output_metadata, paste0('datasets/', study, "/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_metadata.csv"))
  }
}

#cuprizone and EAE KLK6
library(Seurat)
library(tidyverse)
library(DESeq2)
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/PMCB/SCB/hOligo_MS_data")
cup<-readRDS("EBPi_snRNAseq_full_seurat.RDS")
# cup <- FindVariableFeatures(object = cup)
# cup <- ScaleData(object = cup)
# cup <- RunPCA(object = cup)
# cup <- RunUMAP(object = cup, dims = 1:15)
# FeaturePlot(object = cup, features = "Klk6")
# DimPlot(cup, reduction = 'umap', group.by = "sargent_cellstates", raster = FALSE)

View(cup@meta.data)
unique(cup@meta.data$disease)
unique(cup@meta.data$cuprizone)
unique(cup@meta.data$treatment)
unique(cup@meta.data$time..wks.)
sum(is.na(cup@meta.data$cuprizone))
cup@meta.data$cell_type<-cup@meta.data$sargent_cellstates
#cup@meta.data$cell_type<-cup@meta.data$signacx_cellstates
cup@meta.data$sample_id
cup@meta.data$condition<-"control"
cup@meta.data$condition[cup@meta.data$cuprizone == "Cuprizone"]<-"MS"
cup@meta.data$condition[cup@meta.data$treatment == "RA742"]<-"treated"
unique(cup@meta.data$cell_type)
unique(cup@meta.data$condition)

metadata2<-cup@meta.data
seurat.obj<-cup

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

celltypes_list<-unique(metadata2$cell_type)

#split by cell type
for(x in 15:length(celltypes_list)){
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
    write.csv(output_counts, paste0("OligoAnalysis/Cuprizone/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_data.csv"))
    write.csv(output_metadata, paste0("OligoAnalysis/Cuprizone/pseudobulk_counts_matrices/", celltypes_list[x], "_pseudobulk_metadata.csv"))
    data_counts<-output_counts
    data_counts <- round(data_counts)
    metadata<-output_metadata
    metadata$sample <- rownames(metadata)
    dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
    dds <- estimateSizeFactors(dds)
    keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
    dds <- dds[keep,]
    dds <- DESeq(dds)
    result <- results(dds, contrast = c('condition', 'MS', 'control'))
    final_results <- data.frame(result)
    final_results$gene_name <- rownames(final_results)
    final_results$comparison <- "Cuprizone vs control"
    final_results$cell_type<-celltypes_list[x]
    write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_Cuprizone vs Control", final_results$cell_type[1], ".csv"))
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
write.csv(output_counts, paste0("OligoAnalysis/Cuprizone/pseudobulk_counts_matrices/All_Cells_pseudobulk_data.csv"))
write.csv(output_metadata, paste0("OligoAnalysis/Cuprizone/pseudobulk_counts_matrices/All_Cells_pseudobulk_metadata.csv"))
data_counts<-output_counts
data_counts <- round(data_counts)
metadata<-output_metadata
metadata$sample <- rownames(metadata)
dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
dds <- dds[keep,]
dds <- DESeq(dds)
result <- results(dds, contrast = c('condition', 'MS', 'control'))
final_results <- data.frame(result)
final_results$gene_name <- rownames(final_results)
final_results$comparison <- "Cuprizone vs control"
final_results$cell_type<-celltypes_list[x]
write.csv(final_results, paste0("OligoAnalysis/differential_expression_results/diffex_results_Cuprizone vs Control All_Cells.csv"))
counts <- as.data.frame(counts(dds, normalized = TRUE))
counts$gene_name <- rownames(counts)
write.csv(counts, paste0("OligoAnalysis/Cuprizone/pseudobulk_counts_matrices/All_Cells_Normalized pseudobulk_data.csv"))

library(ggprism)

df_counts <- gather(as.data.frame(counts), key = 'sample', value = 'counts', -gene_name)
df_counts <- left_join(df_counts, select(metadata, sample, condition), by = 'sample')

#create barplot for EVERY gene and celltype and indication in datasets
all_genes <- "Klk6"
celltype<- "All_cells"

for(x in 1:length(all_genes)){
  df_counts_plot <- filter(df_counts, gene_name == all_genes[x])
  df_counts_plot_summary <- df_counts_plot %>%
    group_by(condition) %>%
    summarise(mean = mean(counts),
              sd = sd(counts))
  diseases <- as.character(unique(metadata[metadata$condition!='control', 'condition']))
  df_counts_plot_summary <- df_counts_plot_summary %>% mutate(condition = fct_relevel(condition, "control", diseases))
  bar_plot <- ggplot(df_counts_plot_summary, aes(x = condition, y = mean, fill = condition, colour = condition)) +
    geom_bar(stat = 'identity', width = 0.5, alpha = .5) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
    geom_jitter(data = df_counts_plot, aes(x = condition, y = counts, shape = condition), size = 1, width = 0.25) +
    theme_prism(base_size=10) +
    scale_fill_prism(palette='shades_of_gray')+
    scale_colour_prism(palette='shades_of_gray')+
    scale_shape_prism(palette='filled')+
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 10)) +
    ggtitle(paste(celltype, all_genes[x])) +
    ylab('normalized counts')
}

klk6_cup<-list.files("OligoAnalysis/differential_expression_results/")[grepl("uprizone", list.files("OligoAnalysis/differential_expression_results/"))]
find_klk6p<-function(x){
  b<-read.csv(paste0("OligoAnalysis/differential_expression_results/",x))
  print(b$padj[b$X=="Klk6"])
  print(x)
}
for(a in 1:length(klk6_cup)){
  find_klk6p(klk6_cup[a])
}


EAE <- readRDS("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/mouse_datasets/EAE_Internal_Hammond/seurat_object/CSF1r inhibitor mouse EAE.rds")
EAE <- FindVariableFeatures(object = EAE)
EAE <- ScaleData(object = EAE)
EAE <- RunPCA(object = EAE)
EAE <- RunUMAP(object = EAE, dims = 1:15)
FeaturePlot(object = EAE, reduction = 'umap', features = "KLK6")
DimPlot(EAE, reduction = 'umap', group.by = "cell_type", raster = FALSE)

#KLK6 gene correlation
all_cells_counts<-output_counts
all_cells_counts<-read.csv("OligoAnalysis/pseudobulk_counts_matrices/All_Cells(new_object)_pseudobulk_data.csv")
rownames(all_cells_counts)<-all_cells_counts$X
all_cells_counts<-all_cells_counts[-1]
all_cells_counts$baseMean<-rowMeans(all_cells_counts)
keep <- rowSums(all_cells_counts == 0)<21
all_cells_counts <- all_cells_counts[keep,]
all_cells_counts<-all_cells_counts/all_cells_counts$baseMean
all_cells_counts<-all_cells_counts[,1:20]
klk6<-all_cells_counts["KLK6",]

sampleDists <- dist(all_cells_counts)
sampleDistMatrix <- as.matrix(sampleDists)
klk6_df<-as.data.frame(sampleDistMatrix["KLK6",])
colnames(klk6_df)<-"KLK6_correlation"
write.csv(klk6_df, 'OligoAnalysis/KLK6_correlation')

mat<-cor(t(output_counts), method = "pearson", use = "complete.obs")
mat_df<-as.data.frame(mat["KLK6",])
colnames(mat_df)<-"KLK6_Pearson_correlation (stressed and normal oligos from 180759)"

library("Hmisc")
mat2 <- rcorr(as.matrix(t(output_counts)))
mat2r<-as.data.frame(mat2$r["KLK6",])
colnames(mat2r)<-"KLK6 correlation"
mat2P<-as.data.frame(mat2$P["KLK6",])
dim(mat2P)==dim(mat2r)
mat2r$pvalue<-mat2P[,1]
View(mat2r)

keep <- rowSums(all_cells_counts == 0)
for(x in 1:45){
  print(sum(keep==x))
}
