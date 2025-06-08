#### Figure 1: Bacteria and algal growth curve #####
library(readxl)
library(ggplot2)
library(cowplot)
library(ggpubr)

Count=read_excel("IsolateCountResult.xlsx")
Count$Day=as.factor(Count$Day)
Count$Density=as.numeric(Count$Density)

ggplot(data = Count, aes(x = Day, y = Density/1000,col=Treatments)) +
  geom_point() +  stat_summary(aes(group=Treatments),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+ theme_bw()+
  ggtitle("") +
  theme(legend.title = element_blank(),
        legend.key.size = unit(3.5, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.1, 0.95),
        plot.title = element_text(size=7),
        text = element_text(size=7),
        legend.text = element_text(face = "italic"))+
facet_wrap(~ BacteriaID, scales = "free_y",ncol = 2, strip.position = "top")+
  scale_color_manual(values=c("darkgreen", "orange")) 

Single_M=ggplot(data = TM, aes(x = Day, y = Density,col=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+
  geom_vline(xintercept=12,linetype=3,color="blue")+
  annotate(geom = "text", x = 11, y = 11.7, 
           label = "RNA-seq", color = "blue",
           angle = 90, cex=2)+ 
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(c) Treatment F: Single") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)

Mix_H=ggplot(data = Mix_TH, aes(x = Day, y = Density,col=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),fun.y=mean, geom="line")+
  ylab("Million of CFUs/ml")+xlab("Day")+
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(b) Treatment I: Mixed") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)

Mix_M=ggplot(data = Mix_TM, aes(x = Day, y = Density,col=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),  fun.y=mean, geom="line")+ 
  ylab("Million of CFUs/ml")+xlab("Day")+
  scale_color_manual(values=c("black", "pink")) + theme_bw()+
  ggtitle("(d) Treatment F: Mixed") +
  theme(plot.title = element_text(size=8))+
  theme(text = element_text(size=8),
        legend.position = "none") +ylim(0,13)


#####Plot algal growth####
library(gdata)
library(reshape2)
Chla = list.files(pattern="^D\\d+\\.xls$") 
Bacteria=read_excel("PlateArrangment_Run2.xlsx")

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

Algae=ggplot(data = CoMassiveDataChla, aes(x = Day_N, y = Chla, col=Bacteria)) +
  geom_point() +  
  stat_summary(aes(group=Bacteria),  fun.y=mean, geom="line")+ 
  ylab("RFU of chlorophyll-a")+xlab("Day")+
  scale_color_manual(values=c("black", "green","brown","pink")) + theme_bw()+
  ggtitle("(e) C. sorokiniana growth") +
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

