############## Prepare rarefied gene expression table ######
library(readxl)
library(reshape2)
library(dplyr)  # or library(tidyverse)
library(vegan)
####
count_raw=read.table("Data/Falsiroseomonas_counts.txt", row.name=1)
count=count_raw[,3:8]
count_rarefied = t(rrarefy(t(count), sample=2040648))
colnames(count_rarefied)=c("F_1","F_2","F_3","I_1","I_2","I_3")
count_rarefied=data.frame(count_rarefied,rownames(count_rarefied))
colnames(count_rarefied)[7]="GeneID"
write.csv(count_rarefied,"Rarefired_Ros_counts.csv")

count_raw=read.table("Data/Curvibacter_counts.txt", row.name=1)
count=count_raw[,3:8]
count_rarefied = t(rrarefy(t(count), sample=2040648))
colnames(count_rarefied)=c("F_1","F_2","F_3","I_1","I_2","I_3")
count_rarefied=data.frame(count_rarefied,rownames(count_rarefied))
colnames(count_rarefied)[7]="GeneID"
write.csv(count_rarefied,"Rarefired_Cur_counts.csv")

count_raw=read.table("Data/Chlorella_counts.txt", row.name=1)
count=count_raw[,3:11]
count_rarefied = t(rrarefy(t(count), sample=min(colSums(count)))) 
colSums(count_rarefied)
colnames(count_rarefied)=c("C_1","C_2","C_3","R_1","R_2","R_3","A_1","A_2","A_3")
write.csv(count_rarefied,"Chl_count_rarefied.csv")

################## Prepare matched DOC function ##############
COG_DOC_list_raw=read_excel("COG_DOMtransporters/DOC_COG_list.xlsx",col_names = FALSE)
COG_DOC_list_raw=as.data.frame(COG_DOC_list_raw)

Name=numeric()
COG_DOC_list=list()
for (i in 1:9){
  Name[i]=sapply(strsplit(as.character(COG_DOC_list_raw[i,]),":", fixed = TRUE), "[", 1)
  x=sapply(strsplit(as.character(COG_DOC_list_raw[i,]),":", fixed = TRUE), "[", 2)
  COG_DOC=numeric()
  for (j in 1:length(strsplit(x,";", fixed = TRUE)[[1]])){
    COG_DOC[j]=sapply(strsplit(sapply(strsplit(strsplit(x,";", fixed = TRUE)[[1]][j],",", fixed = TRUE), "[", 1)," ", fixed = TRUE), "[", 2)
  }
  COG_DOC_list[[i]]=COG_DOC
}
names(COG_DOC_list)=Name
COG_DOC_list=melt(COG_DOC_list)
colnames(COG_DOC_list)[1]="COG"
write.csv(COG_DOC_list,"COG_DOC_Cat.csv")

##### DOM 
COG_DOC_list=read.csv("COG_DOC_Cat.csv")
## Input annotation data
Data_raw=read_excel("Data/Falsiroseomonas_eggNOG.xlsx",col_names = T)
Annotation_Ros=Data_raw[3:nrow(Data_raw),]
colnames(Annotation_Ros)=Data_raw[2,]
Ros_COG_COG=sapply(strsplit(Annotation_Ros$eggNOG_OGs,"@", fixed = TRUE), "[", 1) 
Ros_COG_Name_COG=data.frame(Annotation_Ros$query,Ros_COG_COG)
Ros_COG_Name_COG=Ros_COG_Name_COG[grep("COG", Ros_COG_Name_COG$Ros_COG_COG),]
colnames(Ros_COG_Name_COG)=c("GeneID","COG")
Ros_COG=merge(COG_DOC_list,Ros_COG_Name_COG,by="COG")

write.csv(Ros_COG,"Ros_DOM_GeneID_Matched.csv")

####

Data_raw=read_excel("Data/Curvibacter_eggNOG.xlsx",col_names = T)
Annotation_Cur=Data_raw[3:nrow(Data_raw),]
colnames(Annotation_Cur)=Data_raw[2,]
Cur_COG_COG=sapply(strsplit(Annotation_Cur$eggNOG_OGs,"@", fixed = TRUE), "[", 1) 
Cur_COG_Name_COG=data.frame(Annotation_Cur$query,Cur_COG_COG)
Cur_COG_Name_COG=Cur_COG_Name_COG[grep("COG", Cur_COG_Name_COG$Cur_COG_COG),]
colnames(Cur_COG_Name_COG)=c("GeneID","COG")
Cur_COG=merge(COG_DOC_list,Cur_COG_Name_COG,by="COG")

write.csv(Cur_COG,"Cur_DOM_GeneID_Matched.csv")

############### Prepare KEGG #########

Annotation_Cur=read.table("Data/Curvibacter_KAAS.txt",header=F,sep="\t")
colnames(Annotation_Cur)=c("GeneID","KO")
Annotation_Ros=read.table("Data/Falsiroseomonas_KAAS.txt",header=F,sep="\t")
colnames(Annotation_Ros)=c("GeneID","KO")

######## Prepare KO and pathway list #####
library(KEGGREST)
U_KO=unique(c(Annotation_Ros$KO))
P_list=list()

for (i in 119:length(U_KO)){
  P=keggGet(U_KO[i])[[1]]$PATHWAY
  x=class(P)=="NULL"
  if (x==TRUE){
    P_list[[i]]=data.frame(rep(U_KO[i],1),"NULL")
  }else{
    P_list[[i]]=data.frame(rep(U_KO[i],length(P)),P)
  }
  colnames(P_list[[i]])=c("KO","Functions")
}

library(dplyr)
KO_Pathway_All = bind_rows(P_list)
#write.csv(KO_Pathway_All,"Bact_RosOnly_KEGG_PathwayName.csv")
write.csv(KO_Pathway_All,"Bact_CurOnly_KEGG_PathwayName.csv")
