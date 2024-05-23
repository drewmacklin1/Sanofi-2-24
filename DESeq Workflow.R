#this code outputs l2fc, pvalue, padj csvs for every dataset for every celltype comparing disease vs control
#this code generates a barplot for every gene for every celltype for every dataset
#this code generates a heatmap and summary file for the l2fc and padj for "program" (short gene lists that are important to a certain program)

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggprism)

targets <- as.data.frame(as.matrix(read.csv('results_target_credentialing/genelist_targets.csv'), row.names = 1))
program_list <- unique(targets$program)[1:7] #ignore PMP
studylist<-c(list.files("datasets")[1:10], list.files("datasets")[13:16]) #ignore the overlap folders
studylist<-studylist[c(3,4,5,7,8,10,11,12,13,14)]

for(s in 1:length(studylist)){
  study <- studylist[s]
  files <- list.files(path= paste0('datasets/', study, '/pseudobulk_counts_matrices'))
  for (y in 1:length(program_list)){
    program1 <- program_list[y]
    dir.create(paste0("results_target_credentialing_121223/", program1))
    dir.create(paste0("results_target_credentialing_121223/", program1, "/", study))
  }

  #for every cell type in pseudobulk folder 
  odds <- function(x) x[ x %% 2 == 1 ]
  modifier <- odds(1:length(files))
  for (j in modifier){
    data_counts<-as.data.frame(read.csv(paste0('datasets/', study, '/pseudobulk_counts_matrices/', files[j])))
    rownames(data_counts)<-data_counts$X
    data_counts<-data_counts[-1]
    data_counts <- round(data_counts)
    metadata<-as.data.frame(read.csv(paste0('datasets/', study, '/pseudobulk_counts_matrices/', files[j+1])))
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
    
    #define disease/treatment states and compare each one to control
    for (i in 1:length(unique(metadata$condition))){
      disease <- as.character(unique(metadata$condition)[i])
      if(disease!='control'){
        result <- results(dds, contrast = c('condition', disease, 'control')) #this does not work for tissue-split datasets unless the objects are re-run through pseudobulk not split by tissue
        
        #Convert Results to dataframe and add celltype identifier and comparison identifier
        final_results <- data.frame(result)
        final_results$celltype <- gsub(".csv", "", gsub("_pseudobulk_data", "", files[j]))
        final_results$gene_name <- rownames(final_results)
        final_results$comparison <- gsub(" ", "_", paste(disease, "vs.",'Control'))
        final_results$dataset.origin <- study
        
        #Export results as csv to folder
        write.csv(final_results, paste0("datasets/", study, "/differential_expression_results/diffex_results_", final_results$celltype[1], final_results$comparison[1], ".csv"))
        
        #Heat Maps and summary file
        paletteLength <- 50
        myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
        for (x in 1:length(program_list)){
          program1 <- program_list[x]
          program_gene_list <- filter(targets, program %in% program_list[x])$gene_name
          heatmap_info <- subset(filter(final_results, gene_name %in% program_gene_list), select = c('log2FoldChange','padj', 'pvalue'))
          not_in <- program_gene_list[!(program_gene_list %in% final_results$gene_name)]
          if(length(not_in)!=0){
            start <- dim(heatmap_info)[1]
            end <- as.integer(length(not_in))
            for (y in 1:end){
              rownames(heatmap_info[start+y,]) <- not_in[y]
            }
          }
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
                           main = paste(program1, ' Log2FC'))
            ggsave(filename = paste0(final_results$comparison[1],"_heatmap.png"),
                   plot = ph,
                   device = 'png',
                   path = paste0("results_target_credentialing_121223/", program1, "/", study, "/", final_results$celltype[1]),
                   width = 8,
                   height = 8,
                   units = "in")
          }
          program_gene_list <- filter(targets, program %in% program_list[x])$gene_name
          summary_info <- subset(filter(final_results, gene_name %in% program_gene_list), select = c('log2FoldChange','padj', 'pvalue', 'celltype', 'comparison', 'dataset.origin'))
          if(!is.na(summary_info[1,1])){
            summary_info$program<-program1
            write.csv(summary_info, paste0("results_target_credentialing_121223/", program1, "/", study, "/", final_results$celltype[1], "/", final_results$comparison[1], "_summary.csv"))
          }
        }
      }
    }
    
    #Prep for graphs
    # counts <- as.data.frame(counts(dds, normalized = TRUE))
    # counts$gene_name <- rownames(counts)
    # write.csv(counts, paste0("datasets/", study, "/normalized_pseudobulk_counts_matrices/", final_results$celltype[1], ".csv"))
    # df_counts <- gather(as.data.frame(counts), key = 'sample', value = 'counts', -gene_name)
    # df_counts <- left_join(df_counts, select(metadata, sample, condition), by = 'sample')
    # 
    #create barplot
    # for (y in 1:length(program_list)){
    #   program1 <- program_list[y]
    #   program_gene_list <- filter(targets, program %in% program1)$gene_name
    #   for(z in 1:length(program_gene_list)){
    #     gene_plot <- program_gene_list[z]
    #     df_counts_plot <- filter(df_counts, gene_name == gene_plot)
    #     df_counts_plot_summary <- df_counts_plot %>%
    #       group_by(condition) %>%
    #       summarise(mean = mean(counts),
    #                 sd = sd(counts))
    #     diseases <- as.character(unique(metadata[metadata$condition!='control', 'condition']))                                   ###NEW
    #     df_counts_plot_summary <- df_counts_plot_summary %>% mutate(condition = fct_relevel(condition, "control", diseases))    ###NEW
    #     bar_plot <- ggplot(df_counts_plot_summary, aes(x = condition, y = mean, fill = condition, colour = condition)) +
    #       geom_bar(stat = 'identity', width = 0.5, alpha = .5) +
    #       geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
    #       geom_jitter(data = df_counts_plot, aes(x = condition, y = counts, shape = condition), size = 1, width = 0.25) +
    #       theme_prism(base_size=10) +
    #       scale_fill_prism(palette='shades_of_gray')+
    #       scale_colour_prism(palette='shades_of_gray')+
    #       scale_shape_prism(palette='filled')+
    #       theme(text = element_text(size = 10),
    #             axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    #             axis.title.x = element_blank(),
    #             legend.position = 'none',
    #             plot.title = element_text(size = 10)) +
    #       ggtitle(gene_plot) +
    #       ylab('normalized counts')
    #     
    #     ggsave(filename = paste0(final_results$celltype[1],"_", gene_plot, "_barplot.png"),
    #            plot = bar_plot,
    #            device = 'png',
    #            path = paste0("results_target_credentialing_121223/", program1, "/", study, "/", final_results$celltype[1]),
    #            width = 4,
    #            height = 4,
    #            units = "in")
    #   }
    # }
  }
}

#database-wide and program summary
library(fs)

total_summary<- NULL
summary_locs<-NULL
programs<- list.files("results_target_credentialing_121223")
programs<- programs[!is_file(paste0("results_target_credentialing_121223/", programs))]
for(a in 1:length(programs)){
  program_summary_locs<-NULL
  studies<-list.files(paste0("results_target_credentialing_121223/", programs[a]))
  studies<- studies[!is_file(paste0("results_target_credentialing_121223/", programs[a], "/", studies))]
  for(b in 1:length(studies)){
    celltypes<-list.files(paste0("results_target_credentialing_121223/", programs[a],"/", studies[b]))
    for(c in 1:length(celltypes)){
      sumsandplots<-list.files(paste0("results_target_credentialing_121223/", programs[a],"/", studies[b], "/", celltypes[c]))
      for(d in 1:length(sumsandplots)){
        if(grepl("summary", sumsandplots[d])){
          summary_locs<-c(summary_locs, paste0("results_target_credentialing_121223/", programs[a],"/", studies[b], "/", celltypes[c], "/", sumsandplots[d]))
          program_summary_locs<-c(program_summary_locs, paste0("results_target_credentialing_121223/", programs[a],"/", studies[b], "/", celltypes[c], "/", sumsandplots[d]))
        }
      }
    }
  }
  if(length(program_summary_locs)>0){
    program_total_summary<-read.csv(program_summary_locs[1])
    if(length(program_summary_locs)>1){
      for(e in 2:length(program_summary_locs)){
        program_temp<-read.csv(program_summary_locs[e])
        program_total_summary<-rbind(program_total_summary,program_temp)
      }
    }
  }
  write.csv(program_total_summary, paste0("results_target_credentialing_121223/", programs[a], "/", programs[a], "_total_summary.csv"))
  program_significant <- subset(filter(program_total_summary, abs(log2FoldChange)>.585))
  program_significant <- subset(filter(program_significant, padj<.05))
  write.csv(program_significant, paste0("results_target_credentialing_121223/", programs[a], "/", programs[a], "_total_summary(significant).csv"))
  print(programs[a])
  total_summary<-rbind(total_summary,program_total_summary)
}
write.csv(total_summary, "results_target_credentialing_121223/total_summary.csv")
significant <- subset(filter(total_summary, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
write.csv(significant, "results_target_credentialing_121223/total_summary(significant).csv")
