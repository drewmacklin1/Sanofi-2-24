#L2FC>.25, rawp<.05, >.66 datasets, ranked by rowSums(abs(L2FC))
library(tidyverse)
#library(UpSetR)
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")

all_genes <- read.csv("results_target_credentialing/all_genes.csv", row.names = 1)
all_genes<-all_genes$x
studylist<-c("AD_EC_SFG_PMID33432193", "ALS_Spinalcord_Internal_Zelic", "ALS_FTLD_PrimMotorCortex_GSE174332", "FTD_FrontTempCortex_PMID35879464", "MS_Cortex_Internal_Proto", "MS_CorticalSubcorticalLesions_PRJNA544731", "MS_Lesions_PMID34497421", "PD_Midbrain_PMID34919646", "PD_PreFCortex_PPR453641", "PD_Putamen_Internal_Unknown", "PD_SubNigra_PMID35513515", "AD_OC_OTC_PMID33609158")
all_celltypes<- c("Astrocyte", "Neuron", "Macrophage", "Endothelial", "Ependymal", "Fibroblast", "OPC", "Microglia", "Oligodendrocyte", "Pericyte", "Lymphocyte")
combos<-c('AD', 'AD_ALS', 'AD_ALS_FTD', 'AD_ALS_FTD_FTLD', 'AD_ALS_FTD_FTLD_LBD', 'AD_ALS_FTD_FTLD_MS', 'AD_ALS_FTD_FTLD_MS_LBD', 'AD_ALS_FTD_FTLD_MS_PD', 'AD_ALS_FTD_FTLD_MS_PD_LBD', 'AD_ALS_FTD_FTLD_PD', 'AD_ALS_FTD_FTLD_PD_LBD', 'AD_ALS_FTD_LBD', 'AD_ALS_FTD_MS', 'AD_ALS_FTD_MS_LBD','AD_ALS_FTD_MS_PD', 'AD_ALS_FTD_MS_PD_LBD', 'AD_ALS_FTD_PD', 'AD_ALS_FTD_PD_LBD', 'AD_ALS_FTLD', 'AD_ALS_FTLD_LBD', 'AD_ALS_FTLD_MS', 'AD_ALS_FTLD_MS_LBD', 'AD_ALS_FTLD_MS_PD', 'AD_ALS_FTLD_MS_PD_LBD', 'AD_ALS_FTLD_PD', 'AD_ALS_FTLD_PD_LBD', 'AD_ALS_LBD', 'AD_ALS_MS', 'AD_ALS_MS_LBD', 'AD_ALS_MS_PD', 'AD_ALS_MS_PD_LBD', 'AD_ALS_PD', 'AD_ALS_PD_LBD', 'AD_FTD', 'AD_FTD_FTLD', 'AD_FTD_FTLD_LBD', 'AD_FTD_FTLD_MS', 'AD_FTD_FTLD_MS_LBD', 'AD_FTD_FTLD_MS_PD', 'AD_FTD_FTLD_MS_PD_LBD', 'AD_FTD_FTLD_PD', 'AD_FTD_FTLD_PD_LBD', 'AD_FTD_LBD', 'AD_FTD_MS', 'AD_FTD_MS_LBD', 'AD_FTD_MS_PD', 'AD_FTD_MS_PD_LBD', 'AD_FTD_PD', 'AD_FTD_PD_LBD', 'AD_FTLD', 'AD_FTLD_LBD', 'AD_FTLD_MS', 'AD_FTLD_MS_LBD', 'AD_FTLD_MS_PD', 'AD_FTLD_MS_PD_LBD', 'AD_FTLD_PD', 'AD_FTLD_PD_LBD', 'AD_LBD', 'AD_MS', 'AD_MS_LBD', 'AD_MS_PD', 'AD_MS_PD_LBD', 'AD_PD', 'AD_PD_LBD', 'ALS', 'ALS_FTD', 'ALS_FTD_FTLD', 'ALS_FTD_FTLD_LBD', 'ALS_FTD_FTLD_MS', 'ALS_FTD_FTLD_MS_LBD', 'ALS_FTD_FTLD_MS_PD', 'ALS_FTD_FTLD_MS_PD_LBD', 'ALS_FTD_FTLD_PD', 'ALS_FTD_FTLD_PD_LBD', 'ALS_FTD_LBD', 'ALS_FTD_MS', 'ALS_FTD_MS_LBD', 'ALS_FTD_MS_PD', 'ALS_FTD_MS_PD_LBD', 'ALS_FTD_PD', 'ALS_FTD_PD_LBD', 'ALS_FTLD', 'ALS_FTLD_LBD', 'ALS_FTLD_MS', 'ALS_FTLD_MS_LBD', 'ALS_FTLD_MS_PD', 'ALS_FTLD_MS_PD_LBD', 'ALS_FTLD_PD', 'ALS_FTLD_PD_LBD', 'ALS_LBD', 'ALS_MS', 'ALS_MS_LBD', 'ALS_MS_PD', 'ALS_MS_PD_LBD', 'ALS_PD', 'ALS_PD_LBD', 'FTD', 'FTD_FTLD', 'FTD_FTLD_LBD', 'FTD_FTLD_MS', 'FTD_FTLD_MS_LBD', 'FTD_FTLD_MS_PD', 'FTD_FTLD_MS_PD_LBD','FTD_FTLD_PD', 'FTD_FTLD_PD_LBD', 'FTD_LBD', 'FTD_MS', 'FTD_MS_LBD', 'FTD_MS_PD', 'FTD_MS_PD_LBD', 'FTD_PD', 'FTD_PD_LBD', 'FTLD', 'FTLD_LBD', 'FTLD_MS', 'FTLD_MS_LBD', 'FTLD_MS_PD', 'FTLD_MS_PD_LBD', 'FTLD_PD', 'FTLD_PD_LBD', 'LBD', 'MS', 'MS_LBD', 'MS_PD', 'MS_PD_LBD', 'PD', 'PD_LBD')
celltype_df <- read.csv("results_target_credentialing/celltype_df.csv")
celltype<- lapply(celltype_df, function(x) x[x!=""])

overlap_table_by_celltype<- data.frame(matrix(ncol = length(all_celltypes), nrow = 127))
colnames(overlap_table_by_celltype) <- all_celltypes
rownames(overlap_table_by_celltype) <- combos

#for all celltypes
for (q in 1:(length(all_celltypes)-1)){ #remove -1 to add lymphocytes
  
  #find all celltype diffex files
  choice<- all_celltypes[q]
  print(choice) #delete this
  diffex_locations <- NULL
  for(i in 1:length(studylist)){
    files <- list.files(path=gsub(" ", '', paste('datasets/', studylist[i], '/differential_expression_results')))
    for(j in 1:length(files)){
      for(k in 1:length(celltype[[choice]])){
        if(grepl(celltype[[choice]][k], files[j])){
          diffex_locations<-c(diffex_locations, paste0('datasets/', studylist[i], '/differential_expression_results/', files[j]))
        }
      }
    }
  }
  diffex_locations<-unique(diffex_locations)

  #add all significant genes to empty dataframe
  data <- data.frame(matrix(ncol = length(diffex_locations), nrow = length(all_genes)))
  colnames(data) <- diffex_locations
  rownames(data) <- all_genes
  for (x in 1:length(diffex_locations)){
    datax<-read.csv(diffex_locations[x])
    significant <- subset(filter(datax, abs(log2FoldChange)>.25))
    significant <- subset(filter(significant, pvalue<.05))
    significant_genes <- significant$X
    if(length(significant_genes)!=0){
      for(z in 1:length(significant_genes)){
        data[significant_genes[z],diffex_locations[x]] <-significant$log2FoldChange[z]
      }
    }
  }
  dataClean<-data[rowSums(is.na(data)) != ncol(data),]
  colnames(dataClean)<-gsub("differential_expression_results/diffex_results_", "", gsub("datasets/","", colnames(dataClean)))
  
  #in 2/3 of disease & group by disease
  in_33_disease<-data.frame(matrix(ncol = 7, nrow = length(rownames(dataClean))))
  colnames(in_33_disease) <- c("AD", "ALS", "FTD", "FTLD", "MS", "PD", "LBD")
  rownames(in_33_disease) <- rownames(dataClean)
  for(i in 1:length(rownames(dataClean))){
    if(sum(is.na(dataClean[i,grepl("AD", colnames(dataClean))]))/length(dataClean[i,grepl("AD", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("AD", colnames(dataClean))])!=0){
      in_33_disease[i,"AD"] <-mean(as.numeric(dataClean[i,grepl("AD", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("ALSvs", colnames(dataClean))]))/length(dataClean[i,grepl("ALSvs", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("ALSvs", colnames(dataClean))])!=0){
      in_33_disease[i,"ALS"] <-mean(as.numeric(dataClean[i,grepl("ALSvs", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("FTLDvs", colnames(dataClean))]))/length(dataClean[i,grepl("FTLDvs", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("FTLDvs", colnames(dataClean))])!=0){
      in_33_disease[i,"FTLD"] <-mean(as.numeric(dataClean[i,grepl("FTLDvs", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("FTD", colnames(dataClean))]))/length(dataClean[i,grepl("FTD", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("FTD", colnames(dataClean))])!=0){
      in_33_disease[i,"FTD"] <-mean(as.numeric(dataClean[i,grepl("FTD", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("MS", colnames(dataClean))]))/length(dataClean[i,grepl("MS", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("MS", colnames(dataClean))])!=0){
      in_33_disease[i,"MS"] <-mean(as.numeric(dataClean[i,grepl("MS", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("PD", colnames(dataClean))]))/length(dataClean[i,grepl("PD", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("PD", colnames(dataClean))])!=0){
      in_33_disease[i,"PD"] <-mean(as.numeric(dataClean[i,grepl("PD", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("LBD", colnames(dataClean))]))/length(dataClean[i,grepl("LBD", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("LBD", colnames(dataClean))])!=0){
      in_33_disease[i,"LBD"] <-mean(as.numeric(dataClean[i,grepl("LBD", colnames(dataClean))]), na.rm=TRUE)
    }
  }
  
  #sort by significance
  in33Clean<-in_33_disease[rowSums(is.na(in_33_disease)) != ncol(in_33_disease),] 
  #lmpyrows<-rownames(in33Clean)### delete after lymphocytes
  keep<-colSums(!is.na(in33Clean))>0
  in33Clean<-in33Clean[,keep]
  #in33Clean<-as.data.frame(in33Clean)### delete after lymphocytes
  #rownames(in33Clean)<-lmpyrows ### delete after lymphocytes
  #colnames(in33Clean)<- "FTD" ### delete after lymphocytes
  in33Clean[,"Significance"] <-rowSums(abs(in33Clean), na.rm = TRUE)
  in33Clean<-in33Clean[order(-in33Clean$Significance), ]
  write.csv(in33Clean, paste0("datasets/overlap/11-14_Avg_disease_Exp_Dir_", choice, ".csv"))
  
  #upset plot
  # astro_upset<-in33Clean[,1:dim(in33Clean)[2]-1] #check the celltype first
  # astro_upset[!is.na(astro_upset[])]<-1
  # astro_upset[is.na(astro_upset[])]<-0
  # plot <- upset(astro_upset, sets = colnames(astro_upset), order.by = "freq")
  # ggsave(filename = paste0(choice, "_UpsetPlot.png"), plot = print(plot), path = "datasets/overlap" , width = 10, height = 8, units ="in", device ="png", dpi=700)
  
  #look for overlaps between disease states
  genecombos<- rep(NA, 10000)
  assoc_combo_names<-NULL
  for(a in 1:length(combos)){
    disease_states<-strsplit(combos[a],"_")
    if(all(disease_states[[1]] %in% colnames(in33Clean))){
      
      if(length(disease_states[[1]])==1){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==2){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==3){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & !is.na(in33Clean[c,disease_states[[1]][3]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!=disease_states[[1]][3] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==4){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & !is.na(in33Clean[c,disease_states[[1]][3]]) & !is.na(in33Clean[c,disease_states[[1]][4]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!=disease_states[[1]][3] & colnames(in33Clean)!=disease_states[[1]][4] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==5){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & !is.na(in33Clean[c,disease_states[[1]][3]]) & !is.na(in33Clean[c,disease_states[[1]][4]]) & !is.na(in33Clean[c,disease_states[[1]][5]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!=disease_states[[1]][3] & colnames(in33Clean)!=disease_states[[1]][4] & colnames(in33Clean)!=disease_states[[1]][5] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==6){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & !is.na(in33Clean[c,disease_states[[1]][3]]) & !is.na(in33Clean[c,disease_states[[1]][4]]) & !is.na(in33Clean[c,disease_states[[1]][5]]) & !is.na(in33Clean[c,disease_states[[1]][6]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!=disease_states[[1]][3] & colnames(in33Clean)!=disease_states[[1]][4] & colnames(in33Clean)!=disease_states[[1]][5] & colnames(in33Clean)!=disease_states[[1]][6] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(disease_states[[1]])==7){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,disease_states[[1]][1]]) & !is.na(in33Clean[c,disease_states[[1]][2]]) & !is.na(in33Clean[c,disease_states[[1]][3]]) & !is.na(in33Clean[c,disease_states[[1]][4]]) & !is.na(in33Clean[c,disease_states[[1]][5]]) & !is.na(in33Clean[c,disease_states[[1]][6]]) & !is.na(in33Clean[c,disease_states[[1]][7]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=disease_states[[1]][1] & colnames(in33Clean)!=disease_states[[1]][2] & colnames(in33Clean)!=disease_states[[1]][3] & colnames(in33Clean)!=disease_states[[1]][4] & colnames(in33Clean)!=disease_states[[1]][5] & colnames(in33Clean)!=disease_states[[1]][6] & colnames(in33Clean)!=disease_states[[1]][7] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
    }
  }
  keep<-colnames(genecombos)== "combo_holder"
  genecombos<-genecombos[,keep]
  colnames(genecombos)[colnames(genecombos)== "combo_holder"]<-assoc_combo_names
  write.csv(genecombos, paste0("datasets/overlap/11-14_genecombos_", choice, ".csv"))
}
combofiles<-list.files("datasets/overlap")[12:22]
for(h in 1:length(combofiles)){
  combofile<-read.csv(paste0("datasets/overlap/", combofiles[h]))
  cell<-gsub("genecombos_", "", gsub( ".csv", "" , combofiles[h]))
  for(k in 2:length(colnames(combofile))){
    overlap_table_by_celltype[colnames(combofile)[k],cell]<-length(unique(combofile[,k]))
  }
}
write.csv(overlap_table_by_celltype, "datasets/overlap/11-14_overlap_table_by_celltype.csv")