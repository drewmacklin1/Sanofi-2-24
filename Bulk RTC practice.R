setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")
  
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggprism)
  
targets <- as.data.frame(as.matrix(read.csv('results_target_credentialing/genelist_targets.csv'), row.names = 1))
program_list <- unique(targets$program)[1:7]
  
studylist<-c("Neuroinflammation_ALOX15_iMG_bulk_20230928_practice")

for(s in 1:length(studylist)){
  study <- studylist[s]
  for (y in 1:length(program_list)){
    program1 <- program_list[y]
    dir.create(paste0("bulk_RTC_practice/", program1))
    dir.create(paste0("bulk_RTC_practice/", program1, "/", study))
  }
  
  final_results_complete<-read.csv(paste0(study, "/", list.files(study)[grepl("differential_expression_complete.csv", list.files(study))]))
  z<-dim(final_results_complete)[2]/3-1
  for(i in 1:z){
    final_results<-final_results_complete[c(3,3*i+1, 3*i+2,3*i+3)]
    final_results$comparison<-gsub("_log2FC", "", colnames(final_results)[2])
    colnames(final_results)<-c("gene_name", "log2FoldChange", "padj", "pvalue", "comparison")
  
    #Heat Maps
    paletteLength <- 50
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    for (x in 1:length(program_list)){
      program1 <- program_list[x]
      program_gene_list <- filter(targets, program %in% program_list[x])$gene_name
      heatmap_info <- subset(filter(final_results, gene_name %in% program_gene_list), select = c('gene_name', 'log2FoldChange','padj', 'pvalue'))
      rownames(heatmap_info)<-heatmap_info$gene_name
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
               path = paste0("bulk_RTC_practice/", program1, "/", study),
               width = 8,
               height = 8,
               units = "in")
      }
      program_gene_list <- filter(targets, program %in% program_list[x])$gene_name
      summary_info <- subset(filter(final_results, gene_name %in% program_gene_list), select = c('gene_name', 'log2FoldChange','padj', 'pvalue', 'comparison'))
      if(!is.na(summary_info[1,1])){
        summary_info$program<-program1
        write.csv(summary_info, paste0("bulk_RTC_practice/", program1, "/", study, "/", final_results$comparison[1], "_summary.csv"))
      }
    }
  }
 
  #Prep for graphs
  counts1<-read.csv(paste0(study, "/", list.files(study)[grepl("counts_normalized.csv", list.files(study))]))
  counts1<-counts1[, order(colnames(counts1))]
  counts1 <- counts1 %>% select("X", "gene_ID", "gene_name", everything())
  counts1<- counts1[-c(1,2)]
  df_counts1 <- gather(as.data.frame(counts1), key = 'sample', value = 'counts', -gene_name)
  df_counts1$condition = substr(df_counts1$sample,1,nchar(df_counts1$sample)-2)
  
  #create barplot
  for (y in 1:length(program_list)){
    program1 <- program_list[y]
    program_gene_list <- filter(targets, program %in% program1)$gene_name
    for(z in 1:length(program_gene_list)){
      gene_plot <- program_gene_list[z]
      df_counts_plot <- filter(df_counts1, gene_name == gene_plot)
      df_counts_plot_summary <- df_counts_plot %>%
        group_by(condition) %>%
        summarise(mean = mean(counts),
                  sd = sd(counts))
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
      ggsave(filename = paste0(study,"_", gene_plot, "_barplot.png"),
             plot = bar_plot,
             device = 'png',
             path = paste0("bulk_RTC_practice/", program1, "/", study),
             width = 4,
             height = 4,
             units = "in")
    }
  }
}

#database-wide and program summary
library(fs)
total_summary<- NULL
summary_locs<-NULL
programs<- list.files("bulk_RTC_practice")
programs<- programs[!is_file(paste0("bulk_RTC_practice/", programs))]
for(a in 1:length(programs)){
  program_summary_locs<-NULL
  studies<-list.files(paste0("bulk_RTC_practice/", programs[a]))
  studies<- studies[!is_file(paste0("bulk_RTC_practice/", programs[a], "/", studies))]
  for(b in 1:length(studies)){
    sumsandplots<-list.files(paste0("bulk_RTC_practice/", programs[a],"/", studies[b]))
    for(d in 1:length(sumsandplots)){
      if(grepl("summary", sumsandplots[d])){
        summary_locs<-c(summary_locs, paste0("bulk_RTC_practice/", programs[a],"/", studies[b], "/", sumsandplots[d]))
        program_summary_locs<-c(program_summary_locs, paste0("bulk_RTC_practice/", programs[a],"/", studies[b], "/", sumsandplots[d]))
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
  write.csv(program_total_summary, paste0("bulk_RTC_practice/", programs[a], "/", programs[a], "_total_summary.csv"))
  program_significant <- subset(filter(program_total_summary, abs(log2FoldChange)>.585))
  program_significant <- subset(filter(program_significant, padj<.05))
  write.csv(program_significant, paste0("bulk_RTC_practice/", programs[a], "/", programs[a], "_total_summary(significant).csv"))
  print(programs[a])
  total_summary<-rbind(total_summary,program_total_summary)
}
write.csv(total_summary, "bulk_RTC_practice/total_summary.csv")
significant <- subset(filter(total_summary, abs(log2FoldChange)>.585))
significant <- subset(filter(significant, padj<.05))
write.csv(significant, "bulk_RTC_practice/total_summary(significant).csv")
