setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage")

library(Seurat)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggprism)

seurat.obj1<-readRDS("Drew_Macklin/OligoAnalysis/seurat_object/GSE180759.rds")
seurat.obj2<-readRDS("Drew_Macklin/OligoAnalysis/seurat_object/GSE118257.rds")
seurat.obj3<-readRDS("Drew_Macklin/OligoAnalysis/seurat_object/PMID31316211.rds")
seurat.obj4<-readRDS("PMCB/SCB/hOligo_MS_data/EBPi_snRNAseq_full_seurat.RDS")
seurat.obj5<-readRDS("Drew_Macklin/mouse_datasets/EAE_Internal_Hammond/seurat_object/CSF1r inhibitor mouse EAE.rds")

###############  filter to just oligos ###################
oligo_only<-function(object, column, name){
  Idents(object)<-object@meta.data[,colnames(object@meta.data)==column]
  object<-subset(object, ident = name)
  return(object)
}

seurat.obj1<-oligo_only(seurat.obj1, "cell.type...subgroup..standardized.", "oligodendrocyte")
seurat.obj2<-oligo_only(seurat.obj2, "cell.type...subgroup..standardized.", "oligodendrocyte")
seurat.obj3<-oligo_only(seurat.obj3, "cell.type...subgroup..standardized.", "oligodendrocyte")
seurat.obj4<-oligo_only(seurat.obj4, "sargent_cellstates", "OLIGODENDROCYTES")
seurat.obj5<-oligo_only(seurat.obj5, "cell_type", "Oligodendrocyte")

###############   #set condition, cell_type, sample_id, etc #############################
View(seurat.obj1@meta.data)
View(seurat.obj2@meta.data)
View(seurat.obj3@meta.data)
View(seurat.obj4@meta.data)
View(seurat.obj5@meta.data)

check<-function(object){
  print(unique(object@meta.data$sample_id))
  print(unique(object@meta.data$cell_type))
  print(unique(object@meta.data$condition))
  print(unique(object@meta.data$gender))
}

#1
seurat.obj1@meta.data$cell_type<-seurat.obj1@meta.data$cell.type...subgroup..standardized.
seurat.obj1@meta.data$sample_id<-seurat.obj1@meta.data$Sample.ID
seurat.obj1@meta.data$condition<-"control"
seurat.obj1@meta.data$condition[seurat.obj1@meta.data$Condition != "Normal control white matter"] <- "MS"
seurat.obj1@meta.data$gender<-seurat.obj1@meta.data$gender..standardized.

#2
seurat.obj2@meta.data$cell_type<-seurat.obj2@meta.data$cell.type...subgroup..standardized.
seurat.obj2@meta.data$sample_id<-seurat.obj2@meta.data$Sample.ID
seurat.obj2@meta.data$condition<-"control"
seurat.obj2@meta.data$condition[seurat.obj2@meta.data$Condition == "Multiple sclerosis (MS)"] <- "MS"
seurat.obj2@meta.data$gender<-seurat.obj2@meta.data$gender..standardized.

#3
seurat.obj3@meta.data$cell_type<-seurat.obj3@meta.data$cell.type...subgroup..standardized.
seurat.obj3@meta.data$sample_id<-seurat.obj3@meta.data$Sample.ID
seurat.obj3@meta.data$condition<-"control"
seurat.obj3@meta.data$condition[seurat.obj3@meta.data$Condition == "Multiple sclerosis (MS)"] <- "MS"
seurat.obj3@meta.data$gender<-seurat.obj3@meta.data$gender..standardized.

#4
seurat.obj4@meta.data$cell_type<-seurat.obj4@meta.data$sargent_cellstates
seurat.obj4@meta.data$condition<-"control"
seurat.obj4@meta.data$condition[seurat.obj4@meta.data$cuprizone == "Cuprizone"]<-"MS"
seurat.obj4@meta.data$condition[seurat.obj4@meta.data$treatment == "RA742"]<-"treated"

#5
seurat.obj5@meta.data$sample_id<-seurat.obj5@meta.data$orig.ident
seurat.obj5@meta.data$condition[seurat.obj5@meta.data$condition == "Naive"]<-"control"
seurat.obj5@meta.data$condition[seurat.obj5@meta.data$condition == "EAE+vehicle"]<-"MS"
seurat.obj5@meta.data$condition[seurat.obj5@meta.data$condition == "EAE+CSF1Ri"]<-"treated"

#pseudobulk and diffex and filter
genelist<-c("ABCC1", "ABLIM3", "ACSBG1", "ACTN2", "ADAMTS12", "ARHGAP26", "BACE2", "C2CD2", "C2CD6", "COL18A1", "CTDSPL", "DHCR24", "ELOVL2", "ENPP6", "FBXL2", "FSTL5", "GAS2", 
            "GTF2IRD1", "HAPLN2", "HEG1", "HIVEP3", "INPP5F", "JPH1", "KANK4", "LIMK2", "LINC02343", "LRCH1", "MACROD2", "MSRA", "MTHFD1L", "NCK2", "PLCH2", "PLEKHA6", "PLEKHG1",
            "PPM1E", "PRTG", "RASGRF1", "RBFOX1", "SAMD4A", "SETBP1", "SGSM1", "SLC25A29", "SLC35F3", "SLC5A11", "TENM4", "TLE2", "TLK1", "TMEM178A", "TNFRSF21", "TRIM36", "USH1C", 
            "VOPP1", "ZSCAN31", "ABCG1", "AC092691", "AL354809", "APCDD1", "ASTN2", "FCHSD2", "FRY", "HHIP", "HIVEP2", "KCNIP4", "LAMA2", "LINC00461", "LINC01792", "LRIG3", "NTRK2", 
            "PLXDC2", "PPFIA2", "PPP3CA", "PRKD1", "RUNX1T1", "SLC6A15", "TMEM132C")

pseudo_and_diffex<-function(seurat.obj, dataset_name){
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
      dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = metadata, design = ~ condition)
      dds <- estimateSizeFactors(dds)
      keep <- rowSums(dds@colData$condition == 'control' & counts(dds, normalized=TRUE)==0) <= .25*sum(dds@colData$condition == 'control') & rowSums(dds@colData$condition != 'control' & counts(dds, normalized=TRUE)==0) <=.25*sum(dds@colData$condition != 'control')
      dds <- dds[keep,]
      dds <- DESeq(dds)
      result <- results(dds, contrast = c('condition', 'MS', 'control'))
      final_results <- data.frame(result)
      final_results$gene_name <- toupper(rownames(final_results))
      final_results$comparison <- "MS vs Control"
      final_results$cell_type<-celltypes_list[x]
      final_results<-final_results[final_results$gene_name %in% genelist, ]
      write.csv(final_results, paste0("Drew_Macklin/OligoAnalysis/differential_expression_results/(75 genes)_MS_vs_control_", dataset_name, ".csv"))
    
      #heatmap
      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      heatmap_info <- subset(final_results, select = c('log2FoldChange','padj', 'pvalue'))
      heatmap_sig <- heatmap_info$pvalue
      heatmap_sig <- ifelse(heatmap_sig <= .05, '*', '')
      heatmap_sig <- as.matrix(heatmap_sig)
      rownames(heatmap_sig) <- rownames(heatmap_info)
      heatmap_info[is.na(heatmap_info)] <- as.double(0)
      if(!(all(heatmap_info==0))){
          myBreaks <- c(seq(min(heatmap_info$log2FoldChange), 0, length.out=ceiling(paletteLength/2) + 1), 
                        seq(max(heatmap_info$log2FoldChange)/paletteLength, max(heatmap_info$log2FoldChange), length.out=floor(paletteLength/2)))
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
          ph <- pheatmap(heatmap_info$log2FoldChange,
                         color=myColor, 
                         breaks=myBreaks,
                         cluster_cols = FALSE,
                         cellwidth=15,
                         cellheight = 15,
                         display_numbers = heatmap_sig,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         labels_row = rownames(heatmap_info),
                         main = paste(dataset_name, ' Log2FC'))
          ggsave(filename = paste0("(75 genes)_MS_vs_control_", dataset_name,"_heatmap.png"),
                 plot = ph,
                 device = 'png',
                 path = paste0("Drew_Macklin/OligoAnalysis/Heatmaps"),
                 width = 8,
                 height = 20,
                 units = "in")
      }
    }
  }
}

pseudo_and_diffex(seurat.obj1, "GSE180759")
pseudo_and_diffex(seurat.obj2, "GSE118257")
pseudo_and_diffex(seurat.obj3, "PMID31316211")
pseudo_and_diffex(seurat.obj4, "Cuprizone")
pseudo_and_diffex(seurat.obj5, "EAE")



#Table
files<-paste0("OligoAnalysis/differential_expression_results/",list.files("OligoAnalysis/differential_expression_results")[1:5])
cuprizone<-read.csv(files[1])
EAE<-read.csv(files[2])
GSE118257<-read.csv(files[3])
GSE180759<-read.csv(files[4])
PMID31316211<-read.csv(files[5])
cuprizone$X<-toupper(cuprizone$X)

diffexfiles<-c("cuprizone", "EAE", "GSE118257", "GSE180759", "PMID31316211")
allDEGs<-unique(c(cuprizone$X, EAE$X, GSE118257$X, GSE180759$X, PMID31316211$X))
colnames<-c("gene", "cuprizone_padj","cuprizone_l2fc", "EAE_padj", "EAE_l2fc", "GSE118257_padj", "GSE118257_l2fc", "GSE180759_padj", "GSE180759_l2fc", "PMID31316211_padj", "PMID31316211_l2fc")
df <- data.frame(matrix(data = NA, ncol = length(colnames), nrow = length(allDEGs)))
rownames(df)<-allDEGs
colnames(df)<-colnames

df$gene[1:length(allDEGs)]<-allDEGs
for(j in 1:length(allDEGs)){
  if(allDEGs[j] %in% cuprizone$X){
    df[j, "cuprizone_l2fc"]<-cuprizone$log2FoldChange[cuprizone$X == allDEGs[j]]
    df[j, "cuprizone_padj"]<-cuprizone$padj[cuprizone$X == allDEGs[j]]
  }
  if(allDEGs[j] %in% EAE$X){
    df[j, "EAE_l2fc"]<-EAE$log2FoldChange[EAE$X == allDEGs[j]]
    df[j, "EAE_padj"]<-EAE$padj[EAE$X == allDEGs[j]]
  }
  if(allDEGs[j] %in% GSE118257$X){
    df[j, "GSE118257_l2fc"]<-GSE118257$log2FoldChange[GSE118257$X == allDEGs[j]]
    df[j, "GSE118257_padj"]<-GSE118257$padj[GSE118257$X == allDEGs[j]]
  }
  if(allDEGs[j] %in% GSE180759$X){
    df[j, "GSE180759_l2fc"]<-GSE180759$log2FoldChange[GSE180759$X == allDEGs[j]]
    df[j, "GSE180759_padj"]<-GSE180759$padj[GSE180759$X == allDEGs[j]]
  }
  if(allDEGs[j] %in% PMID31316211$X){
    df[j, "PMID31316211_l2fc"]<-PMID31316211$log2FoldChange[PMID31316211$X == allDEGs[j]]
    df[j, "PMID31316211_padj"]<-PMID31316211$padj[PMID31316211$X == allDEGs[j]]
  }
}

write.csv(df, paste0("OligoAnalysis/differential_expression_results/(75 genes) Expression Table.csv"))