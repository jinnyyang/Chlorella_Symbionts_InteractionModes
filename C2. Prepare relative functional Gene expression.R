library(readxl)
library(reshape2)
library(dplyr)
library(vegan)
############### Prepare rarefied table contain only DE genes, more expressed under Treatment F and I seperatively
count_rarefied=read.csv("Rarefired_Cur_counts.csv",row.names = 1)
####
DEseq2_Cur=read.csv("Deseq2_Curvibacter.csv")
colnames(DEseq2_Cur)[1]="GeneID"
DEseq2_Cur_DEgenes_GeneID_F=DEseq2_Cur[which(DEseq2_Cur$log2FoldChange>=1 & DEseq2_Cur$padj<=0.05),]$GeneID
DEseq2_Cur_DEgenes_GeneID_I=DEseq2_Cur[which(DEseq2_Cur$log2FoldChange<=-1 & DEseq2_Cur$padj<=0.05),]$GeneID
####
DEgene_count_rarefied_F=count_rarefied[count_rarefied$GeneID %in% DEseq2_Cur_DEgenes_GeneID_F,]
DEgene_count_rarefied_I=count_rarefied[count_rarefied$GeneID %in% DEseq2_Cur_DEgenes_GeneID_I,]
########

############
Cur_COG=read.csv("Cur_DOM_GeneID_Matched.csv",row.names = 1)

Cur_DOM_F=merge(DEgene_count_rarefied_F,Cur_COG,by="GeneID")
Cur_F=melt(data.frame(Cur_DOM_F[2:4],Cur_DOM_F$L1))
colnames(Cur_F)=c("Functions","Treatments","Ratio")

Cur_DOM_I=merge(DEgene_count_rarefied_I,Cur_COG,by="GeneID")
Cur_I=melt(data.frame(Cur_DOM_I[5:7],Cur_DOM_I$L1))
colnames(Cur_I)=c("Functions","Treatments","Ratio")

Cur=rbind(Cur_F,Cur_I)
Cur <- Cur %>%
  group_by(Functions, Treatments) %>%
  summarise(
    DEgene_counts = n(),           # Count the number of rows in each group
    Ratio = sum(Ratio)        # Sum the Ratio values in each group
  )


Replicates=sapply(strsplit(as.character(Cur$Treatments),"_", fixed = TRUE), "[", 2) 
Cur=data.frame(Cur,Replicates)
Cur$Ratio=Cur$Ratio/2040648
write.csv(Cur,"Cur_DOM_Ratio_DEgenes_ForPlot.csv")

####################################
### KEGG #####
Annotation_Cur=read.table("Data/Curvibacter_KAAS.txt",header=F,sep="\t")
colnames(Annotation_Cur)=c("GeneID","KO")
Cur_KEGG_F=merge(Annotation_Cur,DEgene_count_rarefied_F,by="GeneID")
Cur_KEGG_I=merge(Annotation_Cur,DEgene_count_rarefied_I,by="GeneID")

Cur_Forplots_F=melt(Cur_KEGG_F[,2:5])
Cur_Forplots_I=melt(data.frame(Cur_KEGG_I[,2],Cur_KEGG_I[,6:8]))
colnames(Cur_Forplots_F)=c("Functions","Treatments","Ratio")
colnames(Cur_Forplots_I)=c("Functions","Treatments","Ratio")

Cur_Forplots=rbind(Cur_Forplots_F,Cur_Forplots_I)

Cur = Cur_Forplots %>%
  group_by(Functions, Treatments) %>%
  summarise(
    Ratio = sum(Ratio)
  )
Replicates=sapply(strsplit(as.character(Cur$Treatments),"_", fixed = TRUE), "[", 2) 
Cur=data.frame(Cur,Replicates)
Cur$Ratio=Cur$Ratio/2040648
write.csv(Cur,"Cur_KEGG_Ratio_DEgenes_ForPlot.csv")

#########
count_rarefied=read.csv("Rarefired_Ros_counts.csv",row.names = 1)
####
DEseq2_Ros=read.csv("Deseq2_Roseomonas.csv")
colnames(DEseq2_Ros)[1]="GeneID"
DEseq2_Ros_DEgenes_GeneID_F=DEseq2_Ros[which(DEseq2_Ros$log2FoldChange>0 & DEseq2_Ros$padj<=0.05),]$GeneID
DEseq2_Ros_DEgenes_GeneID_I=DEseq2_Ros[which(DEseq2_Ros$log2FoldChange<0 & DEseq2_Ros$padj<=0.05),]$GeneID
####
DEgene_count_rarefied_F=count_rarefied[count_rarefied$GeneID %in% DEseq2_Ros_DEgenes_GeneID_F,]
DEgene_count_rarefied_I=count_rarefied[count_rarefied$GeneID %in% DEseq2_Ros_DEgenes_GeneID_I,]
########

#############
############
Ros_COG=read.csv("Ros_DOM_GeneID_Matched.csv",row.names = 1)

Ros_DOM_F=merge(DEgene_count_rarefied_F,Ros_COG,by="GeneID")
Ros_F=melt(data.frame(Ros_DOM_F[2:4],Ros_DOM_F$L1))
colnames(Ros_F)=c("Functions","Treatments","Ratio")


Ros_DOM_I=merge(DEgene_count_rarefied_I,Ros_COG,by="GeneID")
Ros_I=melt(data.frame(Ros_DOM_I[5:7],Ros_DOM_I$L1))
colnames(Ros_I)=c("Functions","Treatments","Ratio")

Ros=rbind(Ros_F,Ros_I)
Ros <- Ros %>%
  group_by(Functions, Treatments) %>%
  summarise(
    DEgene_counts = n(),           # Count the number of rows in each group
    Ratio = sum(Ratio)        # Sum the Ratio values in each group
  )


Replicates=sapply(strsplit(as.character(Ros$Treatments),"_", fixed = TRUE), "[", 2) 
Ros=data.frame(Ros,Replicates)
Ros$Ratio=Ros$Ratio/2040648
write.csv(Ros,"Ros_DOM_Ratio_DEgenes_ForPlot.csv")

####################################
### KEGG #####
Annotation_Ros=read.table("Data/Falsiroseomonas_KAAS.txt",header=F,sep="\t")
colnames(Annotation_Ros)=c("GeneID","KO")
Ros_KEGG_F=merge(Annotation_Ros,DEgene_count_rarefied_F,by="GeneID")
Ros_KEGG_I=merge(Annotation_Ros,DEgene_count_rarefied_I,by="GeneID")

Ros_Forplots_F=melt(Ros_KEGG_F[,2:5])
Ros_Forplots_I=melt(data.frame(Ros_KEGG_I[,2],Ros_KEGG_I[,6:8]))
colnames(Ros_Forplots_F)=c("Functions","Treatments","Ratio")
colnames(Ros_Forplots_I)=c("Functions","Treatments","Ratio")

Ros_Forplots=rbind(Ros_Forplots_F,Ros_Forplots_I)

Ros = Ros_Forplots %>%
  group_by(Functions, Treatments) %>%
  summarise(
    Ratio = sum(Ratio)
  )
Replicates=sapply(strsplit(as.character(Ros$Treatments),"_", fixed = TRUE), "[", 2) 
Ros=data.frame(Ros,Replicates)
Ros$Ratio=Ros$Ratio/2040648
write.csv(Ros,"Ros_KEGG_Ratio_DEgenes_ForPlot.csv")