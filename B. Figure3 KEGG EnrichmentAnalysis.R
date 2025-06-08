library(gdata)
library(ggplot2)
library(readxl)
library(DESeq2)
library(tidyverse)
library(viridis)
library(ragp)
library(topGO)
library(clusterProfiler)
library(readr)
library(ggh4x)

####
Ros=read.csv("Ros_KEGG_Ratio_DEgenes_ForPlot.csv")
Ros_F=Ros[which(grepl("F", Ros$Treatments)),]
Ros_I=Ros[which(grepl("I", Ros$Treatments)),]

ko_list_F=Ros_F$Functions
ko_list_I=Ros_I$Functions

ekegg_ros_F = enrichKEGG(gene= ko_list_F,organism= "ko",keyType= "kegg")
ekegg_ros_I = enrichKEGG(gene= ko_list_I,organism= "ko",keyType= "kegg")

RosF=dotplot(ekegg_ros_F, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))
RosI=dotplot(ekegg_ros_I, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))

####
Cur=read.csv("Cur_KEGG_Ratio_DEgenes_ForPlot.csv")
Cur_F=Ros[which(grepl("F", Cur$Treatments)),]
Cur_I=Ros[which(grepl("I", Cur$Treatments)),]

ko_list_F=Cur_F$Functions
ko_list_I=Cur_I$Functions

ekegg_ros_F = enrichKEGG(gene= ko_list_F,organism= "ko",keyType= "kegg")
ekegg_ros_I = enrichKEGG(gene= ko_list_I,organism= "ko",keyType= "kegg")

CurF=dotplot(ekegg_ros_F, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))
CurI=dotplot(ekegg_ros_I, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))

FigureS4 = ggarrange(CurI,CurF,common.legend = T,labels=c("Treatment I","Treatment F"))

########Adding count ratio
Chl=read.csv("Chlo_KO_Count.csv",row.names = 1)
Chl_Cur=Chl[which(grepl("C_", Chl$Treatments)),]
Chl_Ros=Chl[which(grepl("R_", Chl$Treatments)),]

ekegg_Chl_Cur = enrichKEGG(gene= Chl_Cur$KO,organism= "ko",keyType= "kegg",
                           pvalueCutoff = 0.05)
ekegg_Chl_Ros = enrichKEGG(gene= Chl_Ros$KO,organism= "ko",keyType= "kegg")

ekegg_Chl_Ros_edited = ekegg_Chl_Ros@result[!grepl("Diabetic cardiomyopathy", ekegg_Chl_Ros@result$Description), ]

ekegg_Chl_Ros@result=ekegg_Chl_Ros_edited

Chl_Cur=dotplot(ekegg_Chl_Cur, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))
Chl_Ros=dotplot(ekegg_Chl_Ros, showCategory = 10)+ theme(axis.text.y = element_text(size = 9))+ theme(axis.text.y = element_text(size = 9))

ggarrange(Chl_Cur,Chl_Ros,common.legend = T,labels=c("+Cur","+Ros"))

library(ggpubr)
ggarrange(RosI,RosF,Chl_Cur,Chl_Ros,common.legend = T,
          labels=c("(A) Falsiroseomonas sp. (Treatment I)",
                   "(B) Falsiroseomonas sp. (Treatment F)",
                   "(C) C. sorokiniana (+Curvibacter sp.)",
                   "(D) C. sorokiniana (+Falsiroseomonas sp.)"),
          font.label = list(size = 12, color = "black", face = "bold")) 

Figure3 = ggarrange(RosI,RosF,Chl_Ros,common.legend = T,labels=c("Ros I","Ros F","Chl_Ros"),nrow=1)

#######
Annotation=read.table("Data/Chlo_KAAS_GreenAlgae.txt",header=F,sep="\t")
colnames(Annotation)=c("GeneID","KO")

DEseq2_Ros=read.csv("Deseq2_Chlorella_FP_F.csv")
colnames(DEseq2_Ros)[1]="GeneID"
DEseq2_Ros_DEgenes_GeneID=DEseq2_Ros[which(DEseq2_Ros$log2FoldChange>0 & DEseq2_Ros$padj<=0.05),]$GeneID
DEseq2_Ax_DEgenes_GeneID=DEseq2_Ros[which(DEseq2_Ros$log2FoldChange<0 & DEseq2_Ros$padj<=0.05),]$GeneID
####
KO_Ros=Annotation[Annotation$GeneID %in% DEseq2_Ros_DEgenes_GeneID,]
KO_Ax=Annotation[Annotation$GeneID %in% DEseq2_Ax_DEgenes_GeneID,]

ekegg_Chl_Ros= enrichKEGG(gene= KO_Ros$KO,organism= "ko",keyType= "kegg",
                           pvalueCutoff = 0.05)
ekegg_Chl_Ax = enrichKEGG(gene= KO_Ax$KO,organism= "ko",keyType= "kegg",
                           pvalueCutoff = 0.05)

Chl_Rps=dotplot(ekegg_Chl_Ros, showCategory = 10)
Chl_Ax=dotplot(ekegg_Chl_Ax, showCategory = 10)

Figure5=ggarrange(Chl_Cur,Chl_Ros,common.legend = T,labels=c("+Ros","Ax"))

