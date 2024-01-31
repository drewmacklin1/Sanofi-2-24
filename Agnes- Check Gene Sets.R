setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

library(tidyverse)

data<-read.csv("Agnes/Agnes Bulk Genes Data.csv")

# data<-data[order(data$ALP_VEH.vs.WT_ALP_VEH.RawP, decreasing = FALSE), ]
# ALP<-head(data$GeneName, 200)
# data<-data[order(data$JCK_VEH.vs.WT_JCK_VEH.RawP, decreasing = FALSE), ]
# JCK<-head(data$GeneName, 200)
# data<-data[order(data$POD_VEH.vs.WT_POD_VEH.RawP, decreasing = FALSE), ]
# POD<-head(data$GeneName, 200)
# all<-unique(c(ALP,JCK,POD)) #517, includes top 200 from each dataset
# all2<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.05 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.05 & data$POD_VEH.vs.WT_POD_VEH.RawP<.05 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >1.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >1.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >1.5]
all2e<-data$geneid[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.05 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.05 & data$POD_VEH.vs.WT_POD_VEH.RawP<.05 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >1.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >1.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >1.5]
# all3<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH <.05 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.05 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.05 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >1.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >1.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >1.5]
# data$total_FDR<-data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH + data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH + data$POD_VEH.vs.WT_POD_VEH.FDR_BH
# data<-data[order(data$total_FDR, decreasing = FALSE), ]
# all4<-head(data$GeneName, 500)

praw.0001<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.0001 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.0001 & data$POD_VEH.vs.WT_POD_VEH.RawP<.0001]
padj.0001<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH<.0001 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.0001 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.0001]
praw.000001<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.000001 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.000001 & data$POD_VEH.vs.WT_POD_VEH.RawP<.000001]
padj.000001<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH<.000001 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.000001 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.000001]
praw.01andl2fc2.5<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.01 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.01 & data$POD_VEH.vs.WT_POD_VEH.RawP<.01 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]
padj.01andl2fc2.5<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH<.01 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.01 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.01 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]
praw.001andl2fc2.5<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.001 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.001 & data$POD_VEH.vs.WT_POD_VEH.RawP<.001 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]
padj.001andl2fc2.5<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH<.001 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.001 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.001 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]
praw.01andl2fc2.25<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.01 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.01 & data$POD_VEH.vs.WT_POD_VEH.RawP<.01 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.25 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.25 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.25]
padj.01andl2fc2.25<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.FDR_BH<.01 & data$JCK_VEH.vs.WT_JCK_VEH.FDR_BH<.01 & data$POD_VEH.vs.WT_POD_VEH.FDR_BH<.01 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.25 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.25 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.25]
praw.0001andl2fc2.5<-data$GeneName[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.0001 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.0001 & data$POD_VEH.vs.WT_POD_VEH.RawP<.0001 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]

praw.01andl2fc2.5_e<-data$geneid[data$ALP_VEH.vs.WT_ALP_VEH.RawP<.01 & data$JCK_VEH.vs.WT_JCK_VEH.RawP<.01 & data$POD_VEH.vs.WT_POD_VEH.RawP<.01 & abs(data$ALP_VEH.vs.WT_ALP_VEH.FoldChange) >2.5 & abs(data$JCK_VEH.vs.WT_JCK_VEH.FoldChange) >2.5 & abs(data$POD_VEH.vs.WT_POD_VEH.FoldChange) >2.5]

verify<-read.csv("Agnes/Agnes Gene list.csv")
FALSE %in% (verify$Ensembl.ID %in% all2e)
sum((all2e  %in% verify$Ensembl.ID) == FALSE)
all2e[(all2e  %in% verify$Ensembl.ID) == FALSE]

FALSE %in% (verify$Ensembl.ID %in% praw.01andl2fc2.5_e)
sum((praw.01andl2fc2.5_e  %in% verify$Ensembl.ID) == FALSE)
praw.01andl2fc2.5_e[(praw.01andl2fc2.5_e %in% verify$Ensembl.ID) == FALSE]

write.csv(praw.01andl2fc2.5, "Agnes/Geneset(praw.01andl2fc2.5).csv")
