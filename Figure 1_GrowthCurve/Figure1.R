#### Figure 1: Bacteria and algal growth curve #####
library(readxl)
library(ggplot2)
library(cowplot)
library(ggpubr)

Count=read_excel("IsolateCount.xlsx")
Count$Density=as.numeric(Count$Density)/1000

TH=Count[Count$Treatments=="D",]
TM=Count[Count$Treatments=="H",]
Mix_TH=Count[Count$Treatments=="D" & Count$Combination=="Mix",]
Mix_TM=Count[Count$Treatments=="H" & Count$Combination=="Mix",]

Single_H=ggplot(data = TH, aes(x = Day, y = Density,col=Bacteria, shape=Bacteria)) +
  geom_point() +  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+
  geom_vline(xintercept=5,linetype=3,color="blue")+
  annotate(geom = "text", x = 4, y = 11.7, 
           label = "RNA-seq", color = "blue",
           angle = 90, cex=2)+ 
  scale_color_manual(values=c("black", "pink"),labels = c("Curvibacter sp.", "Falsiroseomonas sp.")) + 
  scale_shape_manual(values=c(16,17),labels = c("Curvibacter sp.", "Falsiroseomonas sp.")) + 
  theme_bw()+
  ggtitle("(A) Treatment I: Single") +
  theme(legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.7, 0.8),
        plot.title = element_text(size=8),
        text = element_text(size=8),
        legend.text = element_text(face = "italic"))+ylim(0,13)

Single_M=ggplot(data = TM, aes(x = Day, y = Density,col=Bacteria, shape=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+
  geom_vline(xintercept=12,linetype=3,color="blue")+
  annotate(geom = "text", x = 11, y = 11.7, 
           label = "RNA-seq", color = "blue",
           angle = 90, cex=2)+ 
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(C) Treatment F: Single") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)

Mix_H=ggplot(data = Mix_TH, aes(x = Day, y = Density,col=Bacteria, shape=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(B) Treatment I: Mixed") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)

Mix_M=ggplot(data = Mix_TM, aes(x = Day, y = Density,col=Bacteria, shape=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),  fun.y=mean, geom="line")+ 
  ylab("Million of CFUs/ml")+xlab("Day")+
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(D) Treatment F: Mixed") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)


#####Plot algal growth####
library(gdata)
library(reshape2)
Chla = list.files(pattern="^D\\d+\\.xls$") 
Bacteria=read_excel("PlateArrangment.xlsx")

CoDataChla=list()
for (i in 1:length(Chla)){
  rrrawD = lapply(Chla[i], read_excel) ###Just make a list of the plate with different time
  rrawD =as.data.frame(rrrawD)[71,3:50]
  r=as.matrix(rrawD)
  rawD=as.numeric(r) ##Converting data... from dataframe to matrix to numeric...
  Day = sapply(strsplit(Chla[i], ".x"), "[", 1)
  DayRep=rep(Day,length(rawD))
  DayNumeric=as.numeric(rep(sapply(strsplit(Day, "D", fixed = TRUE), "[", 2),length(rawD)))
  x=which(colnames(Bacteria)==Day)
  D=data.frame(Bacteria[,x],rawD,DayNumeric,DayRep)
  colnames(D)=c("Bacteria","Chla","Day_N","Day")
  D=D[D$Bacteria!="Empty",]
  CoDataChla[[i]]=D
  
}
CoMassiveDataChla=melt(CoDataChla,id=1:4)
write.csv(CoMassiveDataChla,"CoMassiveData_Chla.csv")

Algae=ggplot(data = CoMassiveDataChla, aes(x = Day_N, y = Chla, col=Bacteria, shape=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+ 
  ylab("RFU of chlorophyll-a")+xlab("Day")+
  scale_shape_manual(values=c(5,16,4,17),labels = c("Axenic","+Curvibacter sp.","+Mixed symbionts","+Falsiroseomonas sp."))+ 
  scale_color_manual(values=c("darkgreen", "black","brown","pink"),labels = c("Axenic","+Curvibacter sp.","+Mixed symbionts","+Falsiroseomonas sp."))+
  theme_bw()+
  ggtitle("(E) C. sorokiniana growth") +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.25, 0.8),
        plot.title = element_text(size=8),
        text = element_text(size=8),
        legend.text = element_text(face = "italic"))+
  geom_vline(xintercept=12,linetype=3,color="blue")+
  annotate(geom = "text", x = 11, y = 1400, 
           label = "RNA-seq", color = "blue",
           angle = 90, cex=2)

empty_plot = ggdraw() + theme_void()
ggarrange(Single_H, Mix_H,empty_plot,Single_M,Mix_M,Algae,ncol=3, nrow=2)

