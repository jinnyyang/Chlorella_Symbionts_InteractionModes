#### Figure 2: PC analysis on bacteria and algae whole gene expression profile #####
library(DESeq2)
library(ggplot2)

###26D3
RawcountData=read.table("Data/Curvibacter_counts.txt", header = F, sep = "\t")
countData=RawcountData[,4:9]
countData=countData[rowSums(countData)!=0,]
COUNT=as.matrix(countData)
colnames(COUNT)=c("FW","FW","FW","IW","IW","IW")
condition = factor(c("FW","FW","FW","IW","IW","IW"))
dds_Cur = DESeqDataSetFromMatrix(COUNT,DataFrame(condition), ~condition) 
vsd_Cur = vst(dds_Cur, blind=FALSE)
res = results(DESeq(dds_Cur))
write.csv(results(DESeq(dds_Cur)),"Deseq2_Curvibacter.csv")

### Ros
RawcountData=read.table("Data/Roseomonas_counts.txt", header = F, sep = "\t")
countData=RawcountData[,4:9]
countData=countData[rowSums(countData)!=0,]
COUNT=as.matrix(countData)
colnames(COUNT)=c("FP","FP","FP","IP","IP","IP")
condition = factor(c("FP","FP","FP","IP","IP","IP"))
dds_Ros = DESeqDataSetFromMatrix(COUNT,DataFrame(condition), ~condition) 
vsd_Ros = vst(dds_Ros, blind=FALSE)
res_Cur = results(DESeq(dds_Cur))
write.csv(results(DESeq(dds_Ros)),"Deseq2_Roseomonas.csv")

### C. sorokiniana
RawcountData=read.table("Data/Chlorella_counts.txt", row.name=1,header = F, sep = "\t")
countData=RawcountData[,3:11] 
countData=countData[rowSums(countData)>0,]
colSums(countData)
colnames(countData)=c("FW1","FW2","FW3","FP1","FP2","FP3","F1","F2","F3")
condition = factor(c("FW","FW","FW","FP","FP","FP","F","F","F"))
dds_Algae = DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
vsd_Algae = vst(dds_Algae,blind=FALSE) #VST (Variance Stabilizing Transformation)

####
FP_F_count=data.frame(count[,4:9])
FP_F_condition = factor(c("FP","FP","FP","F","F","F"))
condition = relevel(FP_F_condition, ref = "F")  # Sets F as the reference
dds = DESeqDataSetFromMatrix(FP_F_count, DataFrame(FP_F_condition), ~ FP_F_condition)
res = results(DESeq(dds))
write.csv(res,"Deseq2_Chlorella_FP_F.csv")
#
FW_F_count=data.frame(count[,1:3],count[,7:9])
FW_F_condition = factor(c("FW","FW","FW","F","F","F"))
condition = relevel(FW_F_condition, ref = "F")   #F as reference
dds = DESeqDataSetFromMatrix(FW_F_count, DataFrame(FW_F_condition), ~ FW_F_condition)
res = results(DESeq(dds))
write.csv(res,"Deseq2_Chlorella_FW_F.csv")
#
FW_FP_count=data.frame(count[1:6])
FW_FP_condition = factor(c("FW","FW","FW","FP","FP","FP"))
condition = relevel(FW_FP_condition, ref = "FP")   #FP as reference
dds = DESeqDataSetFromMatrix(FW_FP_count, DataFrame(FW_FP_condition), ~ FW_FP_condition)
res = results(DESeq(dds))
write.csv(res,"Deseq2_Chlorella_FW_FP.csv")
####
##### Figure 2 #######
vsd = vst(dds, blind=FALSE)


Cur=plotPCA(vsd_Cur) + theme_bw() + 
  scale_color_manual(values=c("darkgreen", "orange"),
                     labels = c("Treatment F", "Treatment I"),
                     name ="")+
  geom_point(colour="black",pch=21, size=3)+ 
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text=element_text(size=5),
        legend.position = c(0.75, 0.87),
        legend.key.size = unit(0.8, "mm"),
        plot.title = element_text(size=9),
        text = element_text(size=8),
        aspect.ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ylim(-20,20)+
  ggtitle("(A) Curvibacter sp.")
  annotate("text", x=8, y=10, label= "i")

Ros=plotPCA(vsd_Ros) + theme_bw() + 
  scale_color_manual(values=c("darkgreen", "orange"))+
  geom_point(colour="black",pch=21, size=3)+ 
  theme(legend.position = "none",
        plot.title = element_text(size=9),
        text = element_text(size=8),
        aspect.ratio=1)+
  ylim(-20,20)+
  theme(legend.key.size = unit(0.4, "cm"))+
  ggtitle("(B) Falsiroseomonas sp.")

###C. sorokiniana
Algae=plotPCA(vsd_Algae, intgroup = "condition") + theme_bw() + 
  scale_color_manual(values=c("darkgreen", "pink", "black"),
                     labels = c("Axenic", 
                                "+ Falsiroseomonas sp.", "+ Curvibacter sp."),
                     name ="")+
  ylim(-20,20)+
  geom_point(colour="black",pch=21, size=3)+ 
  theme(legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.68, 0.85),
        plot.title = element_text(size=9),
        text = element_text(size=8),
        aspect.ratio=1,
        legend.text = element_text(face = "italic",size=5))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle("(C) C. sorokiniana")

library(ggpubr)

ggarrange(Cur,Ros,Algae,ncol = 3, nrow = 1)
