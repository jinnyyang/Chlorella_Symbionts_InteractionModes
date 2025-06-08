library(DESeq2)
library(ggplot2)
library("pheatmap")

###26D3
count=read.table("Data/Curvibacter_counts.txt", row.name=1)
count=count[,3:8]
condition = factor(c("FW","FW","FW","IW","IW","IW"))

dds = DESeqDataSetFromMatrix(count, DataFrame(condition), ~ condition)
vsd = vst(dds, blind=FALSE)
Dists_W = as.matrix(dist(t(assay(vsd))))


###XPA
count=read.table("Data/Roseomonas_counts.txt", row.name=1)
count=count[,3:8]
condition = factor(c("FP","FP","FP","IP","IP","IP"))
dds = DESeqDataSetFromMatrix(count, DataFrame(condition), ~ condition)
vsd = vst(dds, blind=FALSE)
Dists_P = as.matrix(dist(t(assay(vsd))))


IP_FP=c(Dists_P[4:6,1],Dists_P[4:6,2],Dists_P[4:6,3])
IW_FW=c(Dists_W[4:6,1],Dists_W[4:6,2],Dists_W[4:6,3])
Bact_dis=data.frame(c(IW_FW,IP_FP),c(rep("Curvibacter sp.",9),rep("Falsiroseomonas sp.",9)))
colnames(Bact_dis)=c("Distance","Comparison")

anova_result = aov(Distance ~ Comparison, data = Bact_dis)
summary(anova_result)
tukey_result_B = TukeyHSD(anova_result)

# Display the results
summary(tukey_result_B)

Bacteria=ggplot(Bact_dis,aes(x=Comparison,y=Distance))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ylab("Euclidean distance\n(Treatment I vs F)")+xlab("")+ theme_minimal()+
  annotate(geom="text", x=1, y=90, label="a",color="red")+
  annotate(geom="text", x=2, y=90, label="b",color="red")+
  ggtitle("(A) Bacterial symbionts")

###C. sorokiniana
count_raw=read.table("Data/Chlorella_counts.txt", row.name=1)
countData=count_raw[,3:11]
colnames(countData)=c("FW1","FW2","FW3","FP1","FP2","FP3","F1","F2","F3")
colSums(countData)
condition = factor(c("FW","FW","FW","FP","FP","FP","F","F","F"))
dds = DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
vsd = vst(dds,blind=FALSE) #VST (Variance Stabilizing Transformation)
sampleDists = as.matrix(dist(t(assay(vsd))))

Ros_Ax=c(sampleDists[4:6,1],sampleDists[4:6,2],sampleDists[4:6,3])
Cur_Ax=c(sampleDists[7:9,1],sampleDists[7:9,2],sampleDists[7:9,3])
Cur_Ros=c(sampleDists[7:9,4],sampleDists[7:9,5],sampleDists[7:9,6])

DIST=data.frame(c(Cur_Ax,Ros_Ax,Cur_Ros),
           c(rep("+Curvibacter sp.\nvs\nAxenic",9),rep("+Falsiroseomonas sp.\nvs\nAxenic",9),rep("+Curvibacter sp.\nvs\n+Falsiroseomonas sp.",9)))

colnames(DIST)=c("Distance","Comparison")


anova_result_A = aov(Distance ~ Comparison, data = DIST)
summary(anova_result_A)
tukey_result_A = TukeyHSD(anova_result_A)


Algae=ggplot(DIST,aes(x=Comparison,y=Distance))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ylab("Euclidean distance")+xlab("")+ theme_minimal()+
  annotate(geom="text", x=1, y=40, label="a",color="red")+
  annotate(geom="text", x=2, y=40, label="b",color="red")+
  annotate(geom="text", x=3, y=40, label="c",color="red")+
  ggtitle("(B) C. sorokiniana")
  

library(ggpubr)
ggarrange(Bacteria,Algae, ncol = 2, nrow = 1)


### Output Tukey test
raw_Tueky=rbind(tukey_result_B$Comparison,tukey_result_A$Comparison)
write.csv(raw_Tueky,"Tables/TukeyTest_Raw.csv")

