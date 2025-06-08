ForPlot_DOM=read.csv("Ros_DOM_Ratio_DEgenes_ForPlot.csv",row.names = 1)
ForPlot_KEGG=read.csv("Ros_KEGG_Ratio_DEgenes_ForPlot.csv",row.names = 1)

ForPlot_DOM_Cur=read.csv("Cur_DOM_Ratio_DEgenes_ForPlot.csv",row.names = 1)

Overview=read.csv("Bact_KEGG_RatioSum_FroPlot_Overview.csv",row.names = 1)
AA=read.csv("Bact_KEGG_RatioSum_FroPlot_AA.csv")
Sec=read.csv("Bact_KEGG_RatioSum_FroPlot_Sec.csv")
Cofactors=read.csv("Bact_KEGG_RatioSum_FroPlot_Cofactors.csv")

Chl_Overview=read.csv("Chlo_KEGG_RatioSum_FroPlot_Overview.csv",row.names = 1)
Chl_AA=read.csv("Chlo_KEGG_RatioSum_FroPlot_AA.csv")
Chl_Sec=read.csv("Chlo_KEGG_RatioSum_FroPlot_Sec.csv")
Chl_Cofactors=read.csv("Chlo_KEGG_RatioSum_FroPlot_Cofactors.csv")

library(dplyr)
library(ggplot2)
### FUNCTIONS ####
Plot=function(ForPlot,TITLE){
  # Grouped
  P=ggplot(ForPlot, aes(fill=Treatments, y=Ratio*100, x=Functions)) + 
   geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=c("darkgreen","darkgreen","darkgreen","orange","orange","orange"))+
    theme_minimal() +
    labs(title = TITLE) +
    guides(color = guide_legend(override.aes = list(fill = "white")))+
    theme(legend.key.size = unit(0.4, "cm"),  # Smaller legend boxes
          legend.text = element_text(size = 7),  # Smaller text
          legend.title = element_text(size = 8),
          legend.position = "none",
          plot.title = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 80, hjust = 1),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 7)) +ylab("Relative abundance (%)")
    
  return(P)
}
Plot_Cur=function(ForPlot,TITLE){
  # Grouped
  P=ggplot(ForPlot, aes(fill=Treatments, y=Ratio*100, x=Functions)) + 
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=c("orange","orange","orange"))+
    theme_minimal() +
    labs(title = TITLE) +
    guides(color = guide_legend(override.aes = list(fill = "white")))+
    theme(legend.key.size = unit(0.4, "cm"),  # Smaller legend boxes
          legend.text = element_text(size = 7),  # Smaller text
          legend.title = element_text(size = 8),
          legend.position = "none",
          plot.title = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 80, hjust = 1),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 7)) +ylab("Relative abundance (%)")
  
  return(P)
}
Plot_Chl=function(ForPlot,TITLE){
  P=ggplot(ForPlot, aes(fill=Treatments, y=Ratio*100, x=Functions)) + 
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=c("black","black","black","red","red","red"))+
    theme_minimal() +
    labs(title = TITLE) +
    guides(color = guide_legend(override.aes = list(fill = "white")))+
    theme(legend.key.size = unit(0.4, "cm"),  # Smaller legend boxes
          legend.text = element_text(size = 7),  # Smaller text
          legend.title = element_text(size = 8),
          #legend.position = c(0.85, 0.7),
          legend.position = "none",
          plot.title = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 80, hjust = 1),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 7)) +ylab("Relative abundance (%)")
  
  return(P)
}
adjust_x=function(x){# Initialize an empty vector
  adjust_x <- numeric(x)
  adjust_x[1] <- 0.78
    for(i in 2:x) {
    if (i %% 2 == 0) {
      adjust_x[i] <- adjust_x[i - 1] + 0.45
    } else {
      adjust_x[i] <- adjust_x[i - 1] + 0.55
    }
  }
  return(adjust_x)
  }
WilcoxTest=function(ForPlot_DOM){
  ForPlot_SP=split(ForPlot_DOM,ForPlot_DOM$Functions)
  p=numeric()
  for (j in 1:length(names(ForPlot_SP))){
    data=ForPlot_SP[[j]]
    data$Treatments=sapply(strsplit(data$Treatments,"_", fixed = TRUE), "[", 1) 
    if(length(data$Treatments)>=6){
      p[j]=wilcox.test(Ratio~Treatments,data)$p.value
    }else{
      p[j]=NA
    }
  }
  p = p.adjust(p, method = "fdr")
  names(p)=names(ForPlot_SP)
  p[which(p<=0.1)]="*"
  p[which(p>=0.1)]=""
  return(p)
  }

############

x=adjust_x(length(unique(ForPlot_DOM$Functions))*2)
x=c(x[1:10],6,7,7.78,8.23)##Manually correct....= =
ForPlot_DOM_a=ForPlot_DOM[which(ForPlot_DOM$Replicates==2),]$DEgene_counts
astrick=WilcoxTest(ForPlot_DOM)

astrick_y=numeric()
for(i in 1:8){
  astrick_y[i]=max(split(ForPlot_DOM,ForPlot_DOM$Functions)[[i]]$Ratio)*100
}

TITLE="(B) Falsiroseomonas sp."
Ros_DOM=Plot(ForPlot_DOM,TITLE) + 
  annotate("text", x = x, y = -0.05, label = ForPlot_DOM_a, size = 2.5,col="black") + 
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")
Ros_DOM
c(1,2,0.5,1,0.8,0,0,0.5)
x=adjust_x(length(unique(ForPlot_DOM_Cur$Functions))*2)
ForPlot_DOM_a=ForPlot_DOM_Cur[which(ForPlot_DOM_Cur$Replicates==2),]$DEgene_counts
astrick=WilcoxTest(ForPlot_DOM_Cur)

TITLE="(A) Curvibacter sp."
Cur_DOM=Plot_Cur(ForPlot_DOM_Cur,TITLE) + 
  annotate("text", x = c(1,2,3,4), y = -0.01, label = ForPlot_DOM_a, size = 2,col="black")

ForPlot_DOM_Cur=data.frame(ForPlot_DOM_Cur,"Curvibacter sp.")
ForPlot_DOM=data.frame(ForPlot_DOM,"Falsiroseomonas sp.")
colnames(ForPlot_DOM)[6]="Bact"
colnames(ForPlot_DOM_Cur)[6]="Bact"

StackedBarPLot=rbind(ForPlot_DOM_Cur,ForPlot_DOM)

StackedBarPLot=StackedBarPLot %>%
  group_by(Treatments,Bact) %>%
  mutate(standardized_value = Ratio / sum(Ratio))

Proportion=ggplot(StackedBarPLot, aes(x = Treatments, y = standardized_value, fill = Functions, group)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~Bact , scales = "free", nrow = 1) +  xlab(NULL)+  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.text.x = element_blank()    # Remove x-axis ticks
  )+labs(x = "", y = "Relative abundance proportion", fill = "")+ 
  scale_x_discrete(expand = c(0, 0)) +
      annotate("text", x = 2, y = -0.04, label = "Treatment F", size = 4)+
  annotate("text", x = 5, y = -0.04, label = "Treatment I", size = 4)  

library(ggpubr)
Figure4=ggarrange(Cur_DOM,Ros_DOM,nrow=1)
FigureS5=Proportion
##################################################ß


##### Figure S #################################################
Sec_a=Sec[which(Sec$Replicates==2),]$DEgene_counts
astrick=WilcoxTest(Sec)

astrick_y=numeric()
for(i in 1:length(unique(Sec$Functions))){
  astrick_y[i]=max(split(Sec,Sec$Functions)[[i]]$Ratio)*100
}
TITLE="(A) Biosynthesis of other secondary metabolites"
C=Plot(Sec,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

astrick=WilcoxTest(AA)
astrick_y=numeric()
for(i in 1:length(unique(AA$Functions))){
  astrick_y[i]=max(split(AA,AA$Functions)[[i]]$Ratio)*100
}
TITLE="(B) Amino acid metabolism"
D=Plot(AA,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

astrick=WilcoxTest(Cofactors)
astrick_y=numeric()
for(i in 1:length(unique(Cofactors$Functions))){
  astrick_y[i]=max(split(Cofactors,Cofactors$Functions)[[i]]$Ratio)*100
}
TITLE="(C) Metabolism of cofactors and vitamins"
E=Plot(Cofactors,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

ggarrange(C,D,E,nrow=1) #Figure S5

###
astrick=WilcoxTest(Chl_Sec)

astrick_y=numeric()
for(i in 1:length(unique(Chl_Sec$Functions))){
  astrick_y[i]=max(split(Chl_Sec,Chl_Sec$Functions)[[i]]$Ratio)*100
}

TITLE="(A) Biosynthesis of other secondary metabolites"
Chl_A=Plot_Chl(Chl_Sec,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

astrick=WilcoxTest(Chl_AA)

astrick_y=numeric()
for(i in 1:length(unique(Chl_AA$Functions))){
  astrick_y[i]=max(split(Chl_AA,Chl_AA$Functions)[[i]]$Ratio)*100
}

TITLE="(B) Amino acid metabolism"
Chl_B=Plot_Chl(Chl_AA,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

astrick=WilcoxTest(Chl_Cofactors)

astrick_y=numeric()
for(i in 1:length(unique(Chl_Cofactors$Functions))){
  astrick_y[i]=max(split(Chl_Cofactors,Chl_Cofactors$Functions)[[i]]$Ratio)*100
}
TITLE="(C) Metabolism of cofactors and vitamins"
Chl_C=Plot_Chl(Chl_Cofactors,TITLE)+
  annotate("text", x = 1:length(astrick), y = astrick_y, label = astrick, size = 4,col="blue")

ggarrange(Chl_A,Chl_B,Chl_C,nrow=1) #Figure 5

######## Curvibacter only ####
AA=read.csv("Bact_KEGG_RatioSum_FroPlot_AA_Cur.csv")
Sec=read.csv("Bact_KEGG_RatioSum_FroPlot_Sec_Cur.csv")
Cofactors=read.csv("Bact_KEGG_RatioSum_FroPlot_Cofactors_Cur.csv")
Sec_a=Sec[which(Sec$Replicates==2),]$DEgene_counts
astrick=WilcoxTest(Sec)
TITLE="(A) Biosynthesis of other secondary metabolites"
B=Plot(Sec,TITLE)+
  annotate("text", x = 1:length(astrick), y = 0, label = astrick, size = 6,col="blue")

astrick=WilcoxTest(AA)
TITLE="(B) Amino acid metabolism"
C=Plot(AA,TITLE)+
  annotate("text", x = 1:length(astrick), y = 0, label = astrick, size = 6,col="blue")

astrick=WilcoxTest(Cofactors)
TITLE="(C) Metabolism of cofactors and vitamins"
D=Plot(Cofactors,TITLE)+
  annotate("text", x = 1:length(astrick), y = 0, label = astrick, size = 6,col="blue")

ggarrange(B,C,D,ncol=3,nrow=1) #Figure S6
