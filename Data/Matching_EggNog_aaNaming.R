library(readxl)
###

##### EggNOG ####
EggNOG_1_raw=read_excel("Chl_EggNog_1_Viridplantae.xlsx",col_names = T)
EggNOG_1=EggNOG_1_raw[3:nrow(EggNOG_1_raw),]
colnames(EggNOG_1)=EggNOG_1_raw[2,]
length(unique(EggNOG_1$query)) #10370
nrow(EggNOG_1)

EggNOG_2_raw=read_excel("Chl_EggNog_2_Viridplantae.xlsx",col_names = T)
EggNOG_2=EggNOG_2_raw[3:nrow(EggNOG_2_raw),]
colnames(EggNOG_2)=EggNOG_2_raw[2,]
length(unique(EggNOG_2$query)) #12102
nrow(EggNOG_2)

###
fasta_1=read.table("Script_Fasta/augustus_Split1/augustus_S1.codingseq")
fasta_1_name=fasta_1$V1[grep("^>",fasta_1$V1)]
fasta_1_name=sub("^>", "", fasta_1_name)

fasta_1_g=sapply(strsplit(fasta_1_name,".", fixed = TRUE), "[", 2)
EggNOG_1_g=sapply(strsplit(EggNOG_1$query,".", fixed = TRUE), "[", 1)

10370-3

new_name_1_EggNOG=numeric()
for (i in 1:(length(EggNOG_1_g)-3)){
  new_name_1_EggNOG[i]=fasta_1_name[which(fasta_1_g==EggNOG_1_g[i])]
}
EggNOG_1_named=data.frame(new_name_1_EggNOG,EggNOG_1[1:10367,])
colnames(EggNOG_1_named)[1]="Gene_ID"

fasta=read.table("Script_Fasta/augustus_Split2/augustus_S2.codingseq")
fasta_name=fasta$V1[grep("^>",fasta$V1)]
fasta_name=sub("^>", "", fasta_name)
length(fasta_name) #28903

fasta_g=sapply(strsplit(fasta_name,".", fixed = TRUE), "[", 2)
EggNOG_g=sapply(strsplit(EggNOG_2$query,".", fixed = TRUE), "[", 1)

new_name_2_EggNOG=numeric()
for (i in 1:(nrow(EggNOG_2)-3)){
  new_name_2_EggNOG[i]=fasta_name[which(fasta_g==EggNOG_g[i])]
}
EggNOG_AA2_named=data.frame(new_name_EggNOG[1:12009],EggNOG_2[1:12009,])
colnames(EggNOG_AA2_named)[1]="Gene_ID"

EggNOG_AA_final=rbind(EggNOG_1_named,EggNOG_AA2_named)

write.csv(EggNOG_AA_final,"Chl_EggNOG_renamed.csv")







###############
KASS_AA1=read.table("KAAS/GreenAlgae1SHB.txt")

fasta=read.table("Script_Fasta/augustus_Split1/augustus_S1.codingseq")
fasta_name=fasta$V1[grep("^>",fasta$V1)]
fasta_name=sub("^>", "", fasta_name)

fasta_g=sapply(strsplit(fasta_name,".", fixed = TRUE), "[", 2)
KASS_g=sapply(strsplit(KASS_AA1$V1,".", fixed = TRUE), "[", 1)

new_name_KAAS=numeric()
for (i in 1:length(KASS_g)){
  new_name_KAAS[i]=fasta_name[which(fasta_g==KASS_g[i])]
}
KAAS_AA1_named=data.frame(new_name_KAAS,KASS_AA1)
colnames(KAAS_AA1_named)=c("Gene_ID","Old_ID","Ko")
###
KASS_AA2=read.table("KAAS/GreenAlgae2SHB.txt")

fasta=read.table("Fasta/augustus_Split2/augustus_S2.codingseq")
fasta_name=fasta$V1[grep("^>",fasta$V1)]
fasta_name=sub("^>", "", fasta_name)

fasta_g=sapply(strsplit(fasta_name,".", fixed = TRUE), "[", 2)
KASS_g=sapply(strsplit(KASS_AA2$V1,".", fixed = TRUE), "[", 1)

new_name_KAAS=numeric()
for (i in 1:length(KASS_g)){
  new_name_KAAS[i]=fasta_name[which(fasta_g==KASS_g[i])]
}
KAAS_AA2_named=data.frame(new_name_KAAS,KASS_AA2)
colnames(KAAS_AA2_named)=c("Gene_ID","Old_ID","Ko")

####
KAAS_final=rbind(KAAS_AA1_named,KAAS_AA2_named)
write.csv(KAAS_final,"KASS/Chlorella_KASS_combined_renamed.csv")


