library(DESeq2)
##Curvibacter_counts.csv
count=read.table("Data/Curvibacter_counts.txt", row.name=1)
count=count[,3:8]
condition = factor(c("FW","FW","FW","IW","IW","IW"))  #IW as reference
condition = relevel(condition, ref = "IW")  # Sets IW as the reference
dds = DESeqDataSetFromMatrix(count, DataFrame(condition), ~ condition)
res = results(DESeq(dds))
write.csv(results(DESeq(dds)),"Deseq2_Curvibacter.csv")
DEgenes=res[which(res$padj<=0.05),]
length(DEgenes) #337

length(which(res$log2FoldChange>=1 & res$padj<=0.05)) #8
length(which(res$log2FoldChange<=-1 & res$padj<=0.05)) #295

nrow(res) #6985
####
count=read.table("Data/Roseomonas_counts.txt", row.name=1)
count=count[,3:8]
condition = factor(c("FP","FP","FP","IP","IP","IP"))
condition = relevel(condition, ref = "IP")  # Sets IP as the reference
dds = DESeqDataSetFromMatrix(count, DataFrame(condition), ~ condition)
res = results(DESeq(dds))
write.csv(results(DESeq(dds)),"Deseq2_Roseomonas.csv")
DEgenes=rownames(res)[which(res$padj<=0.05)]
length(DEgenes) #3947'

length(which(res$log2FoldChange>=1 & res$padj<=0.05)) #1424
length(which(res$log2FoldChange<=-1 & res$padj<=0.05)) #1571

nrow(res)
#
count_raw=read.table("Data/Chlorella_counts.txt", row.name=1)
count=count_raw[,3:11]
colnames(count)=c("FW1","FW2","FW3","FP1","FP2","FP3","F1","F2","F3")

FP_F_count=data.frame(count[,4:9])
FP_F_condition = factor(c("FP","FP","FP","F","F","F"))
condition = relevel(FP_F_condition, ref = "F")  # Sets F as the reference
dds = DESeqDataSetFromMatrix(FP_F_count, DataFrame(FP_F_condition), ~ FP_F_condition)
res = results(DESeq(dds))
nrow(res)
length(which(res$log2FoldChange>=1 & res$padj<=0.05)) #6520
length(which(res$log2FoldChange<=-1 & res$padj<=0.05)) #63


write.csv(res,"Deseq2_Chlorella_FP_F.csv")
#
FW_F_count=data.frame(count[,1:3],count[,7:9])
FW_F_condition = factor(c("FW","FW","FW","F","F","F"))
condition = relevel(FW_F_condition, ref = "F")   #F as reference

dds = DESeqDataSetFromMatrix(FW_F_count, DataFrame(FW_F_condition), ~ FW_F_condition)
res = results(DESeq(dds))
write.csv(res,"Deseq2_Chlorella_FW_F.csv")

length(which(res$log2FoldChange>=1 & res$padj<=0.05)) #8676
length(which(res$log2FoldChange<=-1 & res$padj<=0.05)) #47


#
FW_FP_count=data.frame(count[1:6])
FW_FP_condition = factor(c("FW","FW","FW","FP","FP","FP"))
condition = relevel(FW_FP_condition, ref = "FP")   #FP as reference
dds = DESeqDataSetFromMatrix(FW_FP_count, DataFrame(FW_FP_condition), ~ FW_FP_condition)
res = results(DESeq(dds))
write.csv(res,"Deseq2_Chlorella_FW_FP.csv")
DEgenes=rownames(res)[which(res$padj<=0.05)]
DEgenes_Cur=res[which(res$log2FoldChange>=1 & res$padj<=0.05),]
DEgenes_Ros=res[which(res$log2FoldChange<=-1 & res$padj<=0.05),]
nrow(DEgenes_Cur) #8592
nrow(DEgenes_Ros) #6546

