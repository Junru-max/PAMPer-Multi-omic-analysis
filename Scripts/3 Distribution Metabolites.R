library(dplyr)
library(tidyr)
library(gdata)
library(forcats)
library(foreign)
library(Seurat)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(FSA)
library(ggpubr)
library(rstatix)
set.seed(1)
##Distribution for metabolites in PAMPer and MM
##For all metabolites
PAMPer.Meta.stat<-Meta.pdata %>% group_by(SUPER.PATHWAY) %>% count(SUPER.PATHWAY)
PAMPer.Meta.stat<-as.data.frame(PAMPer.Meta.stat)
PAMPer.Meta.stat<-PAMPer.Meta.stat[order(PAMPer.Meta.stat$n),]
PAMPer.Meta.stat$Per<-(100*PAMPer.Meta.stat$n)/sum(PAMPer.Meta.stat$n)
lab<-paste0(PAMPer.Meta.stat$SUPER.PATHWAY, "(",c(1,2,3,3,4,4,24,29,29),"%)")
PAMPer.Meta.stat$SUPER.PATHWAY<-factor(PAMPer.Meta.stat$SUPER.PATHWAY,levels = PAMPer.Meta.stat$SUPER.PATHWAY)
ggpie(PAMPer.Meta.stat, "n", label = lab,fill = "SUPER.PATHWAY", lab.pos = "in")+theme(aspect.ratio = 1)

PAMPer.Meta.stat<-Meta.pdata.MM %>% group_by(SUPER.PATHWAY) %>% count(SUPER.PATHWAY)
PAMPer.Meta.stat<-as.data.frame(PAMPer.Meta.stat)
PAMPer.Meta.stat<-PAMPer.Meta.stat[order(PAMPer.Meta.stat$n),]
PAMPer.Meta.stat$Per<-(100*PAMPer.Meta.stat$n)/sum(PAMPer.Meta.stat$n)
lab<-paste0(PAMPer.Meta.stat$SUPER.PATHWAY, "(",c(1,2,3,4,4,5,25,27,31),"%)")
PAMPer.Meta.stat$SUPER.PATHWAY<-factor(PAMPer.Meta.stat$SUPER.PATHWAY,levels = PAMPer.Meta.stat$SUPER.PATHWAY)
ggpie(PAMPer.Meta.stat, "n", label = lab,fill = "SUPER.PATHWAY", lab.pos = "in")+theme(aspect.ratio = 1)
####Distribution for metabolites in PAMPer For 7 modules
plist<-list()
palette1<-hue_pal()(9)
names(palette1)<-Meta.pdata$SUPER.PATHWAY %>% unique()
for (i in 1:7) {
  PAMPer.Meta.stat<-Meta.pdata %>% filter(module==i) %>% group_by(SUPER.PATHWAY) %>% count(SUPER.PATHWAY)
  PAMPer.Meta.stat<-as.data.frame(PAMPer.Meta.stat)
  PAMPer.Meta.stat<-PAMPer.Meta.stat[order(PAMPer.Meta.stat$n),]
  PAMPer.Meta.stat$Per<-(100*PAMPer.Meta.stat$n)/sum(PAMPer.Meta.stat$n) %>% round(3)
  PAMPer.Meta.stat$Per<- PAMPer.Meta.stat$Per %>% round(1)
  lab<-paste0(PAMPer.Meta.stat$SUPER.PATHWAY, "(",PAMPer.Meta.stat$Per,"%)")
  PAMPer.Meta.stat$SUPER.PATHWAY<-factor(PAMPer.Meta.stat$SUPER.PATHWAY,levels = PAMPer.Meta.stat$SUPER.PATHWAY)
  plist[[i]]<-ggpie(PAMPer.Meta.stat, label =lab ,"Per",fill = "SUPER.PATHWAY", lab.pos = "in",palette=palette1)+theme(aspect.ratio = 1)
  
}
ggarrange(plotlist = plist, 
          ncol = 4, nrow = 2,common.legend = T,align = "v")

