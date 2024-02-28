#1 Human up regulated UpSetR
#2 Human down regulated UpSetR
#3 human datasets heatmap
#4 human sc heatmap
#5 Mouse up regulated UpSetR
#6 Mouse down regulated UpSetR
#7 mouse datasets heatmap
#8 mouse sc datasets heatmap
#9 all datasets heatmap
#10 all datasets (minus MicroArray) heatmap
#11 all datasets pulmonary fibrosis genes heatmap
#12 mouse TGFB1 genes
#13 human TGFB1 genes
#14 mouse TNF genes
#15 human TNF genes
#16 mouse upstream regulators heatmap
#17 human upstream regulators heatmap
#18 all upstream regulators heatmap
#19 Mouse pathway heatmap
#20 Human pathway heatmap

#note: if a gene does not show up - that does not mean there is no L2FC, it means that the associated padj was >.05 so it is not included as it
#would visually lend more confidence to genes where the confidence is undue
#note2: the 4 UpSetR plots must be saved manually

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")
library(tidyverse)
library(pheatmap)
library(UpSetR)

make_heatmap<-function(df,name, height, width){
  anno_df<-annotation_df[rownames(annotation_df) %in% colnames(df),]
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(min(df, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(df, na.rm = TRUE)/paletteLength, max(df, na.rm = TRUE), length.out=floor(paletteLength/2)))
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
  ph <- pheatmap(df,
                 annotation_col = anno_df,
                 color=myColor,
                 breaks=myBreaks,
                 na_col = "black",
                 main = name)
  ggsave(filename = paste0(name, "_heatmap.png"),
         plot = ph,
         device = 'png',
         path = paste0("DMD_datasets/X temp/Heatmaps"),
         width = width,
         height = height,
         units = "in",
         limitsize = FALSE)
}

moderate109<-read.csv("DMD_datasets/X temp/MicroArray/moderate109_clean_diffex.csv")
moderate136<-read.csv("DMD_datasets/X temp/MicroArray/moderate136_clean_diffex.csv")
moderate384<-read.csv("DMD_datasets/X temp/MicroArray/moderate384_clean_diffex.csv")
moderate601<-read.csv("DMD_datasets/X temp/MicroArray/moderate601_clean_diffex.csv")
severe109<-read.csv("DMD_datasets/X temp/MicroArray/severe109_clean_diffex.csv")
severe136<-read.csv("DMD_datasets/X temp/MicroArray/severe136_clean_diffex.csv")
severe384<-read.csv("DMD_datasets/X temp/MicroArray/severe384_clean_diffex.csv")
moderate213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate604<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA604175_Bulk_20240110/differential_expression_results/PRJNA604175_All_results.csv")
moderate658<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA658118_Bulk_20240122/differential_expression_results/PRJNA658118_All_results.csv")
moderate771<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA771932_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
moderate772<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Moderate)_vs._Control.csv")
severe213<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE213925_scRNAseq_20240125/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe218<-read.csv("DMD_datasets/RareMuscle_DMD_External_GSE218201_scRNAseq_20240112/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
severe772<-read.csv("DMD_datasets/RareMuscle_DMD_External_PRJNA772047_scRNAseq_20240105/differential_expression_results/diffex_results_All_CellsDMD_(Severe)_vs._Control.csv")
moderate213$X<-toupper(moderate213$X)
moderate213<-subset(filter(moderate213, padj<.05))
moderate658$X<-toupper(moderate658$X)
moderate658<-subset(filter(moderate658, padj<.05))
moderate771$X<-toupper(moderate771$X)
moderate771<-subset(filter(moderate771, padj<.05))
moderate772$X<-toupper(moderate772$X)
moderate772<-subset(filter(moderate772, padj<.05))
moderate604$X<-toupper(moderate604$X)
moderate604<-subset(filter(moderate604, padj<.05))
severe213$X<-toupper(severe213$X)
severe213<-subset(filter(severe213, padj<.05))
severe218$X<-toupper(severe218$X)
severe218<-subset(filter(severe218, padj<.05))
severe772$X<-toupper(severe772$X)
severe772<-subset(filter(severe772, padj<.05))

diffexfiles<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658", "severe772", "moderate772", "moderate109", "moderate136", "moderate384", "moderate601", "severe109", "severe136", "severe384")
allDEGs<-unique(c(severe213$X, severe218$X, moderate213$X, moderate771$X, moderate658$X, moderate604$X, severe772$X, moderate772$X, moderate109$X, moderate136$X, moderate384$X, moderate601$X, severe109$X, severe136$X, severe384$X))
overlap <- data.frame(matrix(data = NA, ncol = length(diffexfiles), nrow = length(allDEGs)))
colnames(overlap)<-diffexfiles
rownames(overlap)<-allDEGs

for (a in 1:length(allDEGs)){
  if(allDEGs[a] %in% severe213$X){
    overlap$severe213[a]<-severe213$log2FoldChange[severe213$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe218$X){
    overlap$severe218[a]<-severe218$log2FoldChange[severe218$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate213$X){
    overlap$moderate213[a]<-moderate213$log2FoldChange[moderate213$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate771$X){
    overlap$moderate771[a]<-moderate771$log2FoldChange[moderate771$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate604$X){
    overlap$moderate604[a]<-moderate604$log2FoldChange[moderate604$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate658$X){
    overlap$moderate658[a]<-moderate658$log2FoldChange[moderate658$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate772$X){
    overlap$moderate772[a]<-moderate772$log2FoldChange[moderate772$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe772$X){
    overlap$severe772[a]<-severe772$log2FoldChange[severe772$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate109$X){
    overlap$moderate109[a]<-moderate109$logFC[moderate109$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate136$X){
    overlap$moderate136[a]<-moderate136$logFC[moderate136$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate384$X){
    overlap$moderate384[a]<-moderate384$logFC[moderate384$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate601$X){
    overlap$moderate601[a]<-moderate601$logFC[moderate601$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe109$X){
    overlap$severe109[a]<-severe109$logFC[severe109$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe136$X){
    overlap$severe136[a]<-severe136$logFC[severe136$X == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe384$X){
    overlap$severe384[a]<-severe384$logFC[severe384$X == allDEGs[a]]
  }
}

colnames<-c("species","severity", "type", "confidence")
annotation_df<- data.frame(matrix(data = NA, ncol = length(colnames), nrow = length(diffexfiles)))
colnames(annotation_df)<-colnames
rownames(annotation_df)<-diffexfiles
annotation_df[1,]<-c("mouse", "severe", "sc", "med", "213")
annotation_df[2,]<-c("mouse", "severe", "sc", "med", "218")
annotation_df[3,]<-c("mouse", "moderate", "sc", "med", "213")
annotation_df[4,]<-c("mouse", "moderate", "sc", "high", "771")
annotation_df[5,]<-c("mouse", "moderate", "bulk", "med", "604")
annotation_df[6,]<-c("mouse", "moderate", "bulk", "high", "658")
annotation_df[7,]<-c("human", "severe", "sc", "med", "772")
annotation_df[8,]<-c("human", "moderate", "sc", "low", "772")
annotation_df[9,]<-c("human", "moderate", "MicroArray", "high", "109")
annotation_df[10,]<-c("human", "moderate", "MicroArray", "low", "136")
annotation_df[11,]<-c("human", "moderate", "MicroArray", "high", "384")
annotation_df[12,]<-c("human", "moderate", "MicroArray", "high", "601")
annotation_df[13,]<-c("human", "severe", "MicroArray", "low", "109")
annotation_df[14,]<-c("human", "severe", "MicroArray", "med", "136")
annotation_df[15,]<-c("human", "severe", "MicroArray", "low", "384")

#1
human_overlap<-overlap[,c(7:15)]
human_overlap_upset_up<-human_overlap
human_overlap_upset_up[human_overlap_upset_up>=.585]<-1
human_overlap_upset_up[human_overlap_upset_up<.585]<-0
human_overlap_upset_up[is.na(human_overlap_upset_up)]<-0
human_upset_up <- upset(human_overlap_upset_up, sets = colnames(human_overlap_upset_up), order.by = "freq")

#2
human_overlap_upset_down<-human_overlap
human_overlap_upset_down[human_overlap_upset_down>= -.585]<-0
human_overlap_upset_down[human_overlap_upset_down < -.585]<-1
human_overlap_upset_down[is.na(human_overlap_upset_down)]<-0
human_upset_down <- upset(human_overlap_upset_down, sets = colnames(human_overlap_upset_down), order.by = "freq")

#3
human_overlap<-human_overlap[rowSums(is.na(human_overlap))<=dim(human_overlap)[2]/3,]
human_overlap$median<-apply(human_overlap, 1, median, na.rm=T)
human_overlap<-human_overlap[abs(human_overlap$median)>2.25,]
human_overlap<-human_overlap[,-c(length(human_overlap))]
human_all<-make_heatmap(human_overlap, "human_all_expr", height = 20, width = 8)

#4
human_overlap<-overlap[,c(7:8)]
human_overlap<-human_overlap[rowSums(is.na(human_overlap))<=dim(human_overlap)[2]/3,]
human_overlap$median<-apply(human_overlap, 1, median, na.rm=T)
human_overlap<-human_overlap[abs(human_overlap$median)>2.25,]
human_overlap<-human_overlap[,-c(length(human_overlap))]
human_sc<-make_heatmap(human_overlap, "human_sc_expr", height = 15, width = 8)

#5
mouse_overlap<-overlap[,c(1:6)]
mouse_overlap_upset_up<-mouse_overlap
mouse_overlap_upset_up[mouse_overlap_upset_up>=.585]<-1
mouse_overlap_upset_up[mouse_overlap_upset_up<.585]<-0
mouse_overlap_upset_up[is.na(mouse_overlap_upset_up)]<-0
mouse_upset_up <- upset(mouse_overlap_upset_up, sets = colnames(mouse_overlap_upset_up), order.by = "freq")

#6
mouse_overlap_upset_down<-mouse_overlap
mouse_overlap_upset_down[mouse_overlap_upset_down>= -.585]<-0
mouse_overlap_upset_down[mouse_overlap_upset_down < -.585]<-1
mouse_overlap_upset_down[is.na(mouse_overlap_upset_down)]<-0
mouse_upset_down <- upset(mouse_overlap_upset_down, sets = colnames(mouse_overlap_upset_down), order.by = "freq")

#7
mouse_overlap<-mouse_overlap[rowSums(is.na(mouse_overlap))<=dim(mouse_overlap)[2]/3,]
mouse_overlap$median<-apply(mouse_overlap, 1, median, na.rm=T)
mouse_overlap<-mouse_overlap[abs(mouse_overlap$median)>2.25,]
mouse_overlap<-mouse_overlap[,-c(length(mouse_overlap))]
mouse_all<-make_heatmap(mouse_overlap, "mouse_all_expr", height = 60, width = 8)

#8
mouse_overlap<-overlap[,c(1:4)]
mouse_overlap<-mouse_overlap[rowSums(is.na(mouse_overlap))<=dim(mouse_overlap)[2]/3,]
mouse_overlap$median<-apply(mouse_overlap, 1, median, na.rm=T)
mouse_overlap<-mouse_overlap[abs(mouse_overlap$median)>2.25,]
mouse_overlap<-mouse_overlap[,-c(length(mouse_overlap))]
mouse_sc<-make_heatmap(mouse_overlap, "mouse_sc_expr", height = 100, width = 8)

#9
all_overlap<-overlap
all_overlap<-all_overlap[rowSums(is.na(all_overlap))<=dim(all_overlap)[2]/3,]
all_overlap$median<-apply(all_overlap, 1, median, na.rm=T)
all_overlap<-all_overlap[abs(all_overlap$median)>2.25,]
all_overlap<-all_overlap[,-c(length(all_overlap))]
all_both<-make_heatmap(all_overlap, "all_all_expr", height = 15, width = 8)

#10
all_NMA_overlap<-overlap[,c(1:8)]
all_NMA_overlap<-all_NMA_overlap[rowSums(is.na(all_NMA_overlap))<=dim(all_NMA_overlap)[2]/3,]
all_NMA_overlap$median<-apply(all_NMA_overlap, 1, median, na.rm=T)
all_NMA_overlap<-all_NMA_overlap[abs(all_NMA_overlap$median)>2.25,]
all_NMA_overlap<-all_NMA_overlap[,-c(length(all_NMA_overlap))]
all_NMA_overlap<-make_heatmap(all_NMA_overlap[,1:7], "all_No-MicroArray_expr", height = 15, width = 8) #too many NAs in column 8 that dendrogram cant form so throws an error

#11
pulmonaryfibrosis<-read.csv("DMD_datasets/X temp/Pulmonary Fibrosis Pathway.csv")
pf_overlap<-overlap[rownames(overlap) %in% pulmonaryfibrosis$Symbol,]
pf_overlap_plot<-make_heatmap(pf_overlap[,c(1:2,4:15)], "pulmonary fibrosis genes expr", height = 15, width = 8) #too many NAs in column 3 that dendrogram cant form so throws an error

make_heatmap2<-function(df,name, height, width, gene_col){
  anno_df<-annotation_df[rownames(annotation_df) %in% colnames(df),]
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(min(df, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(df, na.rm = TRUE)/paletteLength, max(df, na.rm = TRUE), length.out=floor(paletteLength/2)))
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
  ph <- pheatmap(df,
                 annotation_col = anno_df,
                 annotation_row = gene_col,
                 color=myColor,
                 breaks=myBreaks,
                 na_col = "black",
                 main = name)
  ggsave(filename = paste0(name, "_heatmap.png"),
         plot = ph,
         device = 'png',
         path = paste0("DMD_datasets/X temp/Heatmaps"),
         width = width,
         height = height,
         units = "in",
         limitsize = FALSE)
}

#12
TGFB1<-read.csv("DMD_datasets/X temp/TGFB1 genes.csv")
mouse_TGFB1_overlap<-overlap[rownames(overlap) %in% TGFB1$ID,c(1:2,4:6)]
mouse_TGFB1_overlap<-mouse_TGFB1_overlap[rowSums(is.na(mouse_TGFB1_overlap))<=dim(mouse_TGFB1_overlap)[2]/2,]
mouse_gene_col<-TGFB1[TGFB1$ID %in% rownames(mouse_TGFB1_overlap),c(1,3)]
rownames_mouse_gene_col<-mouse_gene_col$ID
mouse_gene_col<-as.data.frame(mouse_gene_col[,c(-1)])
rownames(mouse_gene_col)<-rownames_mouse_gene_col
colnames(mouse_gene_col)<-"Activation in Pathway"
mouse_TGFB1_overlap_plot<-make_heatmap2(mouse_TGFB1_overlap, "Mouse TGFB1 genes", height = 15, width = 8, gene_col = mouse_gene_col)

#13
human_TGFB1_overlap<-overlap[rownames(overlap) %in% TGFB1$ID,c(7:9,11,13,15)]
human_TGFB1_overlap<-human_TGFB1_overlap[rowSums(is.na(human_TGFB1_overlap))<=dim(human_TGFB1_overlap)[2]/2,]
human_gene_col<-TGFB1[TGFB1$ID %in% rownames(human_TGFB1_overlap),c(1,3)]
rownames_human_gene_col<-human_gene_col$ID
human_gene_col<-as.data.frame(human_gene_col[,c(-1)])
rownames(human_gene_col)<-rownames_human_gene_col
colnames(human_gene_col)<-"Activation in Pathway"
human_TGFB1_overlap_plot<-make_heatmap2(human_TGFB1_overlap, "Human TGFB1 genes", height = 15, width = 8, gene_col = human_gene_col)

#14
TNF<-read.csv("DMD_datasets/X temp/TNF genes.csv")
mouse_TNF_overlap<-overlap[rownames(overlap) %in% TNF$ID,c(1:6)]
mouse_TNF_overlap<-mouse_TNF_overlap[rowSums(is.na(mouse_TNF_overlap))<=dim(mouse_TNF_overlap)[2]/2,]
mouse_gene_col<-TNF[TNF$ID %in% rownames(mouse_TNF_overlap),c(1,3)]
rownames_mouse_gene_col<-mouse_gene_col$ID
mouse_gene_col<-as.data.frame(mouse_gene_col[,c(-1)])
rownames(mouse_gene_col)<-rownames_mouse_gene_col
colnames(mouse_gene_col)<-"Activation in Pathway"
mouse_TNF_overlap_plot<-make_heatmap2(mouse_TNF_overlap, "Mouse TNF genes", height = 15, width = 8, gene_col = mouse_gene_col)

#15
human_TNF_overlap<-overlap[rownames(overlap) %in% TNF$ID,c(7:15)]
human_TNF_overlap<-human_TNF_overlap[rowSums(is.na(human_TNF_overlap))<=dim(human_TNF_overlap)[2]/2,]
human_gene_col<-TNF[TNF$ID %in% rownames(human_TNF_overlap),c(1,3)]
rownames_human_gene_col<-human_gene_col$ID
human_gene_col<-as.data.frame(human_gene_col[,c(-1)])
rownames(human_gene_col)<-rownames_human_gene_col
colnames(human_gene_col)<-"Activation in Pathway"
human_TNF_overlap_plot<-make_heatmap2(human_TNF_overlap, "Human TNF genes", height = 15, width = 8, gene_col = human_gene_col)

################ Upstream regulators ##################
significant <- function (y){
  x<-read.csv(y)
  significant <- subset(filter(x, p.value.of.overlap<.05))
  return(significant)
}

files<-paste0("DMD_datasets/X temp/Upstream_regulators/", list.files("DMD_datasets/X temp/Upstream_regulators"))
moderate109<-significant(files[2])
moderate136<-significant(files[3])
moderate213<-significant(files[4])
moderate384<-significant(files[5])
moderate601<-significant(files[6])
moderate604<-significant(files[7])
moderate658<-significant(files[8])
moderate771<-significant(files[9])
moderate772<-significant(files[10])
severe109<-significant(files[11])
severe136<-significant(files[12])
severe213<-significant(files[13])
severe218<-significant(files[14])
severe384<-significant(files[15])
severe772<-significant(files[16])

diffexfiles<-c("severe213", "severe218", "moderate213", "moderate771", "moderate604", "moderate658", "severe772", "moderate772", "moderate109", "moderate136", "moderate384", "moderate601", "severe109", "severe136", "severe384")
allDEGs<-unique(c(severe213$Upstream.Regulator, severe218$Upstream.Regulator, moderate213$Upstream.Regulator, moderate771$Upstream.Regulator, moderate658$Upstream.Regulator, moderate604$Upstream.Regulator, severe772$Upstream.Regulator, 
                  moderate772$Upstream.Regulator, moderate109$Upstream.Regulator, moderate136$Upstream.Regulator, moderate384$Upstream.Regulator, moderate601$Upstream.Regulator, severe109$Upstream.Regulator, severe136$Upstream.Regulator, severe384$Upstream.Regulator))
overlap <- data.frame(matrix(data = NA, ncol = length(diffexfiles), nrow = length(allDEGs)))
rownames(overlap)<-allDEGs
colnames(overlap)<-diffexfiles
for (a in 1:length(allDEGs)){
  if(allDEGs[a] %in% severe213$Upstream.Regulator){
    overlap$severe213[a]<-severe213$Activation.z.score[severe213$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe218$Upstream.Regulator){
    overlap$severe218[a]<-severe218$Activation.z.score[severe218$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate213$Upstream.Regulator){
    overlap$moderate213[a]<-moderate213$Activation.z.score[moderate213$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate771$Upstream.Regulator){
    overlap$moderate771[a]<-moderate771$Activation.z.score[moderate771$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate604$Upstream.Regulator){
    overlap$moderate604[a]<-moderate604$Activation.z.score[moderate604$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate658$Upstream.Regulator){
    overlap$moderate658[a]<-moderate658$Activation.z.score[moderate658$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe772$Upstream.Regulator){
    overlap$severe772[a]<-severe772$Activation.z.score[severe772$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate772$Upstream.Regulator){
    overlap$moderate772[a]<-moderate772$Activation.z.score[moderate772$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate109$Upstream.Regulator){
    overlap$moderate109[a]<-moderate109$Activation.z.score[moderate109$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate136$Upstream.Regulator){
    overlap$moderate136[a]<-moderate136$Activation.z.score[moderate136$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate384$Upstream.Regulator){
    overlap$moderate384[a]<-moderate384$Activation.z.score[moderate384$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% moderate601$Upstream.Regulator){
    overlap$moderate601[a]<-moderate601$Activation.z.score[moderate601$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe109$Upstream.Regulator){
    overlap$severe109[a]<-severe109$Activation.z.score[severe109$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe136$Upstream.Regulator){
    overlap$severe136[a]<-severe136$Activation.z.score[severe136$Upstream.Regulator == allDEGs[a]]
  }
  if(allDEGs[a] %in% severe384$Upstream.Regulator){
    overlap$severe384[a]<-severe384$Activation.z.score[severe384$Upstream.Regulator == allDEGs[a]]
  }
}

#16
mouse_overlap<-overlap[,1:6]
mouse_overlap<-mouse_overlap[rowSums(is.na(mouse_overlap))<=dim(mouse_overlap)[2]/3,]
mouse_overlap$median<-apply(mouse_overlap, 1, median, na.rm=T)
mouse_overlap<-mouse_overlap[abs(mouse_overlap$median)>2.25,]
mouse_overlap<-mouse_overlap[,-c(length(mouse_overlap))]
mouse_USR<-make_heatmap(mouse_overlap, "mouse_USR", height = 60, width = 8)

#17
human_overlap<-overlap[,7:15]
human_overlap<-human_overlap[rowSums(is.na(human_overlap))<=dim(human_overlap)[2]/3,]
human_overlap$median<-apply(human_overlap, 1, median, na.rm=T)
human_overlap<-human_overlap[abs(human_overlap$median)>2.25,]
human_overlap<-human_overlap[,-c(length(human_overlap))]
human_USR<-make_heatmap(human_overlap, "human_USR", height = 60, width = 8)

#18
all_overlap<-overlap
all_overlap<-all_overlap[rowSums(is.na(all_overlap))<=dim(all_overlap)[2]/3,]
all_overlap$median<-apply(all_overlap, 1, median, na.rm=T)
all_overlap<-all_overlap[abs(all_overlap$median)>2.25,]
all_overlap<-all_overlap[,-c(length(all_overlap))]
all_USR<-make_heatmap(all_overlap, "all_USR", height = 45, width = 8)

#19
mouse_pathway<-read.csv("DMD_datasets/X temp/Mouse heatmap info (pathways).csv")
rownames(mouse_pathway)<-mouse_pathway$Canonical.Pathways
mouse_pathway<-mouse_pathway[,c(-1)]
colnames(mouse_pathway)<-c("severe213", "severe218", "moderate658", "moderate213", "moderate604", "moderate771")
for(x in 1:dim(mouse_pathway)[2]){
  mouse_pathway[,x]<-as.numeric(mouse_pathway[,x])  
}
mouse_pathway<-mouse_pathway[rowSums(is.na(mouse_pathway))<=dim(mouse_pathway)[2]/3,]
mouse_pathway$median<-apply(mouse_pathway, 1, median, na.rm=T)
mouse_pathway<-mouse_pathway[abs(mouse_pathway$median)>2.25,]
mouse_pathway<-mouse_pathway[,-c(length(mouse_pathway))]
mouse_pathway_plot<-make_heatmap(mouse_pathway, "mouse_pathways", height = 15, width = 15)

#20
human_pathway<-read.csv("DMD_datasets/X temp/Human heatmap info (pathways).csv")
rownames(human_pathway)<-human_pathway$Canonical.Pathways
human_pathway<-human_pathway[,c(-1)]
colnames(human_pathway)<-c("moderate109", "moderate384", "severe136", "severe384", "moderate136", "severe772", "severe109", "moderate601", "moderate772")
for(x in 1:dim(human_pathway)[2]){
  human_pathway[,x]<-as.numeric(human_pathway[,x])  
}
human_pathway<-human_pathway[rowSums(is.na(human_pathway))<=dim(human_pathway)[2]/3,]
human_pathway$median<-apply(human_pathway, 1, median, na.rm=T)
human_pathway<-human_pathway[abs(human_pathway$median)>2.25,]
human_pathway<-human_pathway[,-c(length(human_pathway))]
human_pathway_plot<-make_heatmap(human_pathway, "human_pathways", height = 15, width = 15)

#21 
all_pathway<-read.csv("All heatmap info (pathways).csv")
rownames(all_pathway)<-all_pathway$Canonical.Pathways
all_pathway<-all_pathway[,c(-1)]
colnames(all_pathway)<-c("moderate109", "moderate658", "moderate384", "moderate771", "severe136", "moderate604", "severe384", "severe218", "moderate136", "severe772", "severe109", "moderate601", "moderate213", "severe213", "moderate772")
for(x in 1:dim(all_pathway)[2]){
  all_pathway[,x]<-as.numeric(all_pathway[,x])  
}
all_pathway<-all_pathway[rowSums(is.na(all_pathway))<=dim(all_pathway)[2]/3,]
all_pathway$median<-apply(all_pathway, 1, median, na.rm=T)
all_pathway<-all_pathway[abs(all_pathway$median)>2.25,]
all_pathway<-all_pathway[,-c(length(all_pathway))]
all_pathway_plot<-make_heatmap(all_pathway, "All_pathways", height = 15, width = 15)
