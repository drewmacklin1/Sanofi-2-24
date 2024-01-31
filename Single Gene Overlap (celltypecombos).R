#L2FC>.25, rawp<.05, >.66 datasets, ranked by rowSums(abs(L2FC))
library(tidyverse)
library(UpSetR)
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")

all_genes <- read.csv("results_target_credentialing/all_genes.csv", row.names = 1)
all_genes<-all_genes$x
studylist<-c("AD_EC_SFG_PMID33432193", "ALS_Spinalcord_Internal_Zelic", "ALS_FTLD_PrimMotorCortex_GSE174332", "FTD_FrontTempCortex_PMID35879464", "MS_Cortex_Internal_Proto", "MS_CorticalSubcorticalLesions_PRJNA544731", "MS_Lesions_PMID34497421", "PD_Midbrain_PMID34919646", "PD_PreFCortex_PPR453641", "PD_Putamen_Internal_Unknown", "PD_SubNigra_PMID35513515")
all_celltypes<- c("Astrocyte", "Ependymal", "Neuron", "Microglia", "Oligodendrocyte", "OPC")
combos<-c('Astro', 'Astro_Epen', 'Epen', 'Astro_Neuro', 'Astro_Epen_Neuro', 'Epen_Neuro', 'Neuro', 'Astro_Micro', 'Astro_Epen_Micro', 'Epen_Micro', 'Astro_Neuro_Micro', 'Astro_Epen_Neuro_Micro', 'Epen_Neuro_Micro', 'Neuro_Micro', 'Micro', 'Astro_Oligo', 'Astro_Epen_Oligo', 'Epen_Oligo', 'Astro_Neuro_Oligo', 'Astro_Epen_Neuro_Oligo', 'Epen_Neuro_Oligo', 'Neuro_Oligo', 'Astro_Micro_Oligo', 'Astro_Epen_Micro_Oligo', 'Epen_Micro_Oligo', 'Astro_Neuro_Micro_Oligo', 'Astro_Epen_Neuro_Micro_Oligo', 'Epen_Neuro_Micro_Oligo', 'Neuro_Micro_Oligo', 'Micro_Oligo', 'Oligo', 'Astro_OPC', 'Astro_Epen_OPC', 'Epen_OPC', 'Astro_Neuro_OPC', 'Astro_Epen_Neuro_OPC', 'Epen_Neuro_OPC', 'Neuro_OPC', 'Astro_Micro_OPC', 'Astro_Epen_Micro_OPC', 'Epen_Micro_OPC', 'Astro_Neuro_Micro_OPC', 'Astro_Epen_Neuro_Micro_OPC', 'Epen_Neuro_Micro_OPC', 'Neuro_Micro_OPC', 'Micro_OPC', 'Astro_Oligo_OPC', 'Astro_Epen_Oligo_OPC', 'Epen_Oligo_OPC', 'Astro_Neuro_Oligo_OPC', 'Astro_Epen_Neuro_Oligo_OPC', 'Epen_Neuro_Oligo_OPC', 'Neuro_Oligo_OPC', 'Astro_Micro_Oligo_OPC', 'Astro_Epen_Micro_Oligo_OPC', 'Epen_Micro_Oligo_OPC', 'Astro_Neuro_Micro_Oligo_OPC', 'Astro_Epen_Neuro_Micro_Oligo_OPC', 'Epen_Neuro_Micro_Oligo_OPC', 'Neuro_Micro_Oligo_OPC', 'Micro_Oligo_OPC', 'Oligo_OPC', 'OPC')
celltype_df <- read.csv("results_target_credentialing/celltype_df.csv")
celltype<- lapply(celltype_df, function(x) x[x!=""])
diseasestate_df <- read.csv("results_target_credentialing/disease_state.csv")
diseasestate<- lapply(diseasestate_df, function(x) x[x!=""])
  
overlap_table_by_diseasestate<- data.frame(matrix(ncol = length(colnames(diseasestate_df)), nrow = length(combos)))
colnames(overlap_table_by_diseasestate) <- colnames(diseasestate_df)
rownames(overlap_table_by_diseasestate) <- combos

#for all disease states
for (q in 1:(length(diseasestate))){ 
  
  #find all celltype diffex files
  choice<- colnames(diseasestate_df)[q]
  diffex_locations <- NULL
  for(i in 1:length(studylist)){
    files <- list.files(path=gsub(" ", '', paste('datasets/', studylist[i], '/differential_expression_results')))
    for(j in 1:length(files)){
      for(k in 1:length(diseasestate[[q]])){
        if(grepl(diseasestate[[q]][k], files[j])){
          diffex_locations<-c(diffex_locations, paste0('datasets/', studylist[i], '/differential_expression_results/', files[j]))
        }
      }
    }
  }
  
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
  in_33_disease<-data.frame(matrix(ncol = 6, nrow = length(rownames(dataClean))))
  colnames(in_33_disease) <- c("Astro", "Epen", "Neuro", "Micro", "Oligo", "OPC")
  rownames(in_33_disease) <- rownames(dataClean)
  for(i in 1:length(rownames(dataClean))){
    if(sum(is.na(dataClean[i,grepl("stro", colnames(dataClean))]))/length(dataClean[i,grepl("stro", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("stro", colnames(dataClean))])!=0){
      in_33_disease[i,"Astro"] <-mean(as.numeric(dataClean[i,grepl("stro", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("pend", colnames(dataClean))]))/length(dataClean[i,grepl("pend", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("pend", colnames(dataClean))])!=0){
      in_33_disease[i,"Epend"] <-mean(as.numeric(dataClean[i,grepl("pend", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i, grepl(paste(c("euro", "xcitatory", "nhibitory", "GABA"), collapse = "|"), colnames(dataClean))]))/length(dataClean[i, grepl(paste(c("euro", "xcitatory", "nhibitory", "GABA"), collapse = "|"), colnames(dataClean))])<=.34 & length(dataClean[i, grepl(paste(c("euro", "xcitatory", "nhibitory", "GABA"), collapse = "|"), colnames(dataClean))])!=0){
      in_33_disease[i,"Neuro"] <-mean(as.numeric(dataClean[i, grepl(paste(c("euro", "xcitatory", "nhibitory", "GABA"), collapse = "|"), colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("icro", colnames(dataClean))]))/length(dataClean[i,grepl("icro", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("icro", colnames(dataClean))])!=0){
      in_33_disease[i,"Micro"] <-mean(as.numeric(dataClean[i,grepl("icro", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl("ligo", colnames(dataClean))]))/length(dataClean[i,grepl("ligo", colnames(dataClean))])<=.34 & length(dataClean[i,grepl("ligo", colnames(dataClean))])!=0){
      in_33_disease[i,"Oligo"] <-mean(as.numeric(dataClean[i,grepl("ligo", colnames(dataClean))]), na.rm=TRUE)
    }
    if(sum(is.na(dataClean[i,grepl(paste(c("OPC", "opc"), collapse = "|"), colnames(dataClean))]))/length(dataClean[i,grepl(paste(c("OPC", "opc"), collapse = "|"), colnames(dataClean))])<=.34 & length(dataClean[i,grepl(paste(c("OPC", "opc"), collapse = "|"), colnames(dataClean))])!=0){
      in_33_disease[i,"OPC"] <-mean(as.numeric(dataClean[i,grepl(paste(c("OPC", "opc"), collapse = "|"), colnames(dataClean))]), na.rm=TRUE)
    }
  }
  
  #sort by significance
  in33Clean<-in_33_disease[rowSums(is.na(in_33_disease)) != ncol(in_33_disease),] 
  keep<-colSums(!is.na(in33Clean))>0
  in33Clean<-in33Clean[,keep]
  in33Clean[,"Significance"] <-rowSums(abs(in33Clean), na.rm = TRUE)
  in33Clean<-in33Clean[order(-in33Clean$Significance), ]
  write.csv(in33Clean, paste0("datasets/overlap_celltypecombos/Avg_disease_Exp_Dir_", choice, ".csv"))
  
  # # upset plot
  # in33Clean<-read.csv("datasets/overlap_celltypecombos/Avg_disease_Exp_Dir_PD.csv")
  # astro_upset<-in33Clean[,1:dim(in33Clean)[2]-1][-1]
  # astro_upset[!is.na(astro_upset[])]<-1
  # astro_upset[is.na(astro_upset[])]<-0
  # upset(astro_upset, sets = colnames(astro_upset), order.by = "freq")
  # #ggsave(filename = paste0(choice, "_UpsetPlot.png"), plot = print(plot), path = "datasets/overlap" , width = 10, height = 8, units ="in", device ="png", dpi=700)
  
  #look for overlaps between disease states
  genecombos<- rep(NA, 10000)
  assoc_combo_names<-NULL
  for(a in 1:length(combos)){
    cell_states<-strsplit(combos[a],"_")
    if(all(cell_states[[1]] %in% colnames(in33Clean))){
      
      if(length(cell_states[[1]])==1){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==2){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==3){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & !is.na(in33Clean[c,cell_states[[1]][3]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!=cell_states[[1]][3] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==4){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & !is.na(in33Clean[c,cell_states[[1]][3]]) & !is.na(in33Clean[c,cell_states[[1]][4]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!=cell_states[[1]][3] & colnames(in33Clean)!=cell_states[[1]][4] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==5){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & !is.na(in33Clean[c,cell_states[[1]][3]]) & !is.na(in33Clean[c,cell_states[[1]][4]]) & !is.na(in33Clean[c,cell_states[[1]][5]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!=cell_states[[1]][3] & colnames(in33Clean)!=cell_states[[1]][4] & colnames(in33Clean)!=cell_states[[1]][5] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==6){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & !is.na(in33Clean[c,cell_states[[1]][3]]) & !is.na(in33Clean[c,cell_states[[1]][4]]) & !is.na(in33Clean[c,cell_states[[1]][5]]) & !is.na(in33Clean[c,cell_states[[1]][6]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!=cell_states[[1]][3] & colnames(in33Clean)!=cell_states[[1]][4] & colnames(in33Clean)!=cell_states[[1]][5] & colnames(in33Clean)!=cell_states[[1]][6] & colnames(in33Clean)!="Significance"])]))){
            combo_holder <-c(combo_holder, rownames(in33Clean)[c])
          }
        }
        genecombos<-cbind(genecombos,combo_holder)
        if(length(combo_holder)>0){
          assoc_combo_names<-c(assoc_combo_names, combos[a])
        }
      }
      
      if(length(cell_states[[1]])==7){
        combo_holder<-NULL
        for(c in 1:length(rownames(in33Clean))){
          if(!is.na(in33Clean[c,cell_states[[1]][1]]) & !is.na(in33Clean[c,cell_states[[1]][2]]) & !is.na(in33Clean[c,cell_states[[1]][3]]) & !is.na(in33Clean[c,cell_states[[1]][4]]) & !is.na(in33Clean[c,cell_states[[1]][5]]) & !is.na(in33Clean[c,cell_states[[1]][6]]) & !is.na(in33Clean[c,cell_states[[1]][7]]) & all(is.na(in33Clean[c,colnames(in33Clean[colnames(in33Clean)!=cell_states[[1]][1] & colnames(in33Clean)!=cell_states[[1]][2] & colnames(in33Clean)!=cell_states[[1]][3] & colnames(in33Clean)!=cell_states[[1]][4] & colnames(in33Clean)!=cell_states[[1]][5] & colnames(in33Clean)!=cell_states[[1]][6] & colnames(in33Clean)!=cell_states[[1]][7] & colnames(in33Clean)!="Significance"])]))){
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
  write.csv(genecombos, paste0("datasets/overlap_celltypecombos/genecombos_", choice, ".csv"))
}
combofiles<-list.files("datasets/overlap_celltypecombos")[8:14]
for(h in 1:length(combofiles)){
  combofile<-read.csv(paste0("datasets/overlap_celltypecombos/", combofiles[h]))
  cell<-gsub("genecombos_", "", gsub( ".csv", "" ,combofiles[h]))
  for(k in 2:length(colnames(combofile))){
    overlap_table_by_diseasestate[colnames(combofile)[k],cell]<-length(unique(combofile[,k]))
  }
}
write.csv(overlap_table_by_diseasestate, "datasets/overlap_celltypecombos/overlap_table_by_celltype.csv")