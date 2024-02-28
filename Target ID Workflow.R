library(biomaRt)
library(tidyverse)
library(httr)
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")

files <- list.files(path = 'datasets/overlap')
ensembl = useEnsembl(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     mirror = "useast")

for(a in c(22, 27, 28, 29)){
  dis_exp_dir <- read.csv(paste0('datasets/overlap/', files[a]))
  significant <- subset(filter(dis_exp_dir, Significance > 2))
  significant$drugName <- NA
  significant$isApproved <- NA
  significant$ensemblID <- NA
  significant$phenotype <-NA
  significant$evidence <-NA
  
  #get ensembl ID
  for (b in 1:length(rownames(significant))) {
    id <- getBM(
      attributes = 'ensembl_gene_id',
      filters = 'hgnc_symbol',
      values = significant$X[b],
      mart = ensembl
    )[[1]]
    if (!is.logical(id)) {
      significant$ensemblID[b] <- paste(id[[1]])
      print(paste(a, "of 11", b, "of", length(rownames(significant))))
    }
  }
  
  for (b in 1:dim(significant)[1]) {
    gene_id <- significant$ensemblID[[b]][1] #only checks against first ensembl id
    
    # Build query string to get info about known drugs
    query_string = "
    query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
     approvedSymbol
     approvedName
     knownDrugs{
      rows{
        drug{
          id
          name
          isApproved
          }
        }
      }
    }
    }
    "
    # Set base URL of GraphQL API endpoint, # Set variables object of arguments to be passed to endpoint,  # Construct POST request body object with query string and variables,  # Perform POST request
    base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
    variables <- list("ensemblId" = gene_id)
    post_body <- list(query = query_string, variables = variables)
    r <- POST(url = base_url, body = post_body, encode = 'json')
    
    if (length(content(r)$data$target$knownDrugs$rows) > 0) {
      drugs <- NULL
      status <- NULL
      for (c in 1:length(content(r)$data$target$knownDrugs[[1]])) {
        drugs <-paste(drugs,content(r)$data$target$knownDrugs$rows[[c]]$drug$name,",")
        status <- paste(status,content(r)$data$target$knownDrugs$rows[[c]]$drug$isApproved,",")
      }
      significant$drugName[b] <- drugs
      significant$isApproved[b] <- status
    }
    
    # get information about GiTIPs
    query_string = "query InternalGitip($ensemblId: String!) {
    internalGitip(ensemblId: $ensemblId) {
      data {
        phenotype
        causal_evidence
        }
      }
    }
    "
    base_url <- "https://opentargets.sanofi.com/api/v4/graphql" #GiTIPs is on the private Sanofi server (genetics informed target indication pairs)
    variables <- list("ensemblId" = gene_id)
    post_body <- list(query = query_string, variables = variables)
    httr::set_config(config(ssl_verifypeer = 0))                    #required to get around Sanofi Proxy
    r <- POST(url = base_url, body = post_body, encode = 'json')
    
    if (length(content(r)$data$internalGitip$data) > 0) {
      phenos <- NULL
      evidence <- NULL
      for (c in 1:length(content(r)$data$internalGitip$data)) {
        if(grepl("lzheimer",content(r)$data$internalGitip$data[[c]]$phenotype) | grepl("myotrophic",content(r)$data$internalGitip$data[[c]]$phenotype) | grepl("ALS",content(r)$data$internalGitip$data[[c]]$phenotype) | grepl("arkinson",content(r)$data$internalGitip$data[[c]]$phenotype) && !(content(r)$data$internalGitip$data[[c]]$causal_evidence %in% c("NO EVIDENCE", "VERY WEAK", "WEAK"))){
          phenos <-paste(phenos,content(r)$data$internalGitip$data[[c]]$phenotype,",")
          evidence <- paste(evidence,content(r)$data$internalGitip$data[[c]]$causal_evidence,",")
        }
      }
      if(!is.null(phenos)){
        significant$phenotype[b] <- phenos
        significant$evidence[b] <- evidence
      }
    }
  }
  significant$phenotype2<-NA
  significant$evidence2<-NA
  for(d in 1:dim(significant)[1]){
    if(grepl("lzheimer", significant$phenotype[d])){
      significant$phenotype2[d] <-"Alzheimer's Disease"
    }
    if(grepl("myotrophic", significant$phenotype[d])){
      significant$phenotype2[d] <-"ALS"
    }
    if(grepl("arkinson", significant$phenotype[d])){
      significant$phenotype2[d] <-"Parkinson's Disease"
    }
    if(grepl("EVIDENCE", significant$evidence[d])){
      significant$evidence2[d] <-1
    }
    if(grepl("WEAK", significant$evidence[d])){ #captures weak and very weak so do first
      significant$evidence2[d] <-3
    }
    if(grepl("VERY WEAK", significant$evidence[d])){
      significant$evidence2[d] <-2
    }
    if(grepl("MODERATE", significant$evidence[d])){
      significant$evidence2[d] <-4
    }
    if(grepl("HIGH", significant$evidence[d])){ #captures high and very high so do first
      significant$evidence2[d] <-5
    }
    if(grepl("VERY HIGH", significant$evidence[d])){
      significant$evidence2[d] <-6
    }
  }
  write.csv(significant, paste0("opentarget/",gsub("Avg_disease_Exp_Dir_", "", files[a])))
}

########################################################################################################################################################
#query TNE for info on total score and druggability
setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")
library(readr)
celltype_files<- c("11-14_Astrocyte", "11-14_Microglia", "11-14_Neuron", "11-14_Oligodendrocyte")

#do everything below for each celltype
celltype<-celltype_files[4] #change for each celltype
csv <- read_csv(paste0("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/opentarget/", celltype, ".csv"))
list<-NULL
for(a in 1:length(csv$X)){
  list<-c(list, paste0("'",csv$X[a], "',"))
}
copy<-c('`Summary_View_General`.`hgnc_symbol` in (',list, ")")
cat(copy, file = paste0("opentarget/", celltype, ".txt"))

#open the txt file in Drew_Macklin/opentarget. Delete the comma and space after the last gene of interest and copy the file.
# go to https://pmcb.sanofi.com/TNETEST/query. Open filter icon. paste the results into the textbox, send query. Click export.
#rename file to "*Celltype_TNEresults" and upload to magellan opentarget folder

TNEresults<-read_csv(paste0("opentarget/", celltype, "_TNEresults.csv"))
TNEresults_filter<-as.data.frame(cbind(TNEresults$`Ensembl Gene ID`, TNEresults$`HGNC Symbol`, TNEresults$`PD PD Total Score`, TNEresults$`AZD AZD Total Score`, TNEresults$`ALS ALS Total Score`, TNEresults$`Druggability - Biologics`, TNEresults$`Druggability - Small Molecule`, TNEresults$`Sanofi Drugs`))
colnames(TNEresults_filter)<-c("`Ensembl Gene ID`", "`HGNC Symbol`", "`PD PD Total Score`", "`AZD AZD Total Score`", "`ALS ALS Total Score`", "`Druggability - Biologics`", "`Druggability - Small Molecule`", "`Sanofi Drugs`")
csv2<-merge(csv, TNEresults_filter, by.x='X', by.y='`HGNC Symbol`')
csv2<-csv2[-2]
csv2$Target_PD<-NA
csv2$Target_AD<-NA
csv2$Target_ALS<-NA
csv2$Celltype<-celltype

#add information about triculture expression
tri_diffex<-list.files("opentarget/Triculture/differential_expression_results")
astro<-read_csv(paste0("opentarget/Triculture/differential_expression_results/",tri_diffex[4]))
astro<-astro[,c(1,2)]
colnames(astro)<-c("HGNC", "Astro BaseMean")
csv2<-merge(csv2, astro, by.x='X', by.y='HGNC')
micro<-read_csv(paste0("opentarget/Triculture/differential_expression_results/",tri_diffex[7]))
micro<-micro[,c(1,2)]
colnames(micro)<-c("HGNC", "Micro BaseMean")
csv2<-merge(csv2, micro, by.x='X', by.y='HGNC')
neuro<-read_csv(paste0("opentarget/Triculture/differential_expression_results/",tri_diffex[10]))
neuro<-neuro[,c(1,2)]
colnames(neuro)<-c("HGNC", "Neuro BaseMean")
csv2<-merge(csv2, neuro, by.x='X', by.y='HGNC')
if(!("PD" %in% colnames(csv2))){
  csv2$PD<-NA
}
if(!("MS" %in% colnames(csv2))){
  csv2$MS<-NA
}
if(!("LBD" %in% colnames(csv2))){
  csv2$LBD<-NA
}
if(!("AD" %in% colnames(csv2))){
  csv2$AD<-NA
}
if(!("ALS" %in% colnames(csv2))){
  csv2$ALS<-NA
}
if(!("FTD" %in% colnames(csv2))){
  csv2$FTD<-NA
}
if(!("FTLD" %in% colnames(csv2))){
  csv2$FTLD<-NA
}

#calculate druggability score
csv2$DruggabilityScore<-1
for(a in 1:length(rownames(csv2))){
  if(grepl("FALSE", csv2$isApproved[a]) | csv2[a,"`Druggability - Small Molecule`"] == "Y" | csv2[a,"`Druggability - Biologics`"] == "Y"){
    csv2$DruggabilityScore[a]<-2
  } 
  if(grepl("TRUE", csv2$isApproved[a])){
    csv2$DruggabilityScore[a]<-3
  } 
  
  #factor in triculture expression
  #if(csv2[a,"Micro BaseMean"]>40){   #change for each celltype
      csv2$Target_PD[a]<-abs(csv2$PD[a])*as.numeric(csv2[a, "`ALS ALS Total Score`"])*csv2$DruggabilityScore[a]
      csv2$Target_AD[a]<-abs(csv2$AD[a])*as.numeric(csv2[a, "`ALS ALS Total Score`"])*csv2$DruggabilityScore[a]
      csv2$Target_ALS[a]<-abs(csv2$ALS[a])*as.numeric(csv2[a, "`ALS ALS Total Score`"])*csv2$DruggabilityScore[a]
 #   }
  
  #calculate evidence score
  if(is.na(csv2$evidence2[a])){
    csv2$evidence2[a]<-1
  }
  if(grepl("ALS", csv2$phenotype2[a])){
    csv2$Target_ALS[a]<-csv2$Target_ALS[a]*(csv2$evidence2[a])
  }
  if(grepl("Alzheimer's Disease", csv2$phenotype2[a])){
    csv2$Target_AD[a]<-csv2$Target_AD[a]*(csv2$evidence2[a])
  }
  if(grepl("Parkinson's Disease", csv2$phenotype2[a])){
    csv2$Target_PD[a]<-csv2$Target_PD[a]*(csv2$evidence2[a])
  }
}
reorder<-c("X", "phenotype2", "evidence2", "DruggabilityScore", "ensemblID","`Ensembl Gene ID`", "Celltype", "Target_PD", "Target_AD", "Target_ALS", "PD", "`PD PD Total Score`", "AD", "`AZD AZD Total Score`", "ALS", "`ALS ALS Total Score`", "drugName", "isApproved","Astro BaseMean","Micro BaseMean", "Neuro BaseMean", "`Druggability - Small Molecule`","`Druggability - Biologics`", "`Sanofi Drugs`", "FTD", "FTLD", "MS", "LBD", "Significance", "phenotype", "evidence")
csv2<-csv2[,reorder]
write.csv(csv2,"opentarget/11-14_OligodendrocyteDeliverable.csv") #change for each celltype

############################################################################################################################################################
#combine all csvs into one
final_m<-read_csv("opentarget/11-14_MicrogliaDeliverable.csv")
final_n<-read_csv("opentarget/11-14_NeuronDeliverable.csv")
final_a<-read_csv("opentarget/11-14_AstrocyteDeliverable.csv")
final_o<-read_csv("opentarget/11-14_OligodendrocyteDeliverable.csv")
final<-rbind(final_m, final_n, final_a, final_o)
final<-final[,-c(1)]
final$Significance2<-NA
for(a in 1:length(rownames(final))){
  final$Significance2[a]<-sum(final$Target_PD[a], final$Target_AD[a], final$Target_ALS[a], na.rm = TRUE)
}
final<-final[order(-final$Significance2), ]
topids<-head(final, 500)
topids<-topids[,c("X", "Celltype", "Target_AD", "Target_ALS", "Target_PD", "phenotype2", "evidence2", "DruggabilityScore", "AD", "`AZD AZD Total Score`", "ALS", "`ALS ALS Total Score`", "PD", "`PD PD Total Score`", "drugName", "isApproved","Astro BaseMean","Micro BaseMean", "Neuro BaseMean", "`Druggability - Small Molecule`","`Druggability - Biologics`", "`Sanofi Drugs`", "ensemblID", "FTD", "FTLD", "MS", "LBD", "phenotype", "evidence", "Significance2")]
topids<-topids[!duplicated(topids), ]
write.csv(topids,"opentarget/11-14_TargetIds.csv")
