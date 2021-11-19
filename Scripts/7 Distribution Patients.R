library(dplyr)
library(tidyr)
library(circlize)
library(dendsort)
library(FSA)
library(ggpubr)
###
Pamper.pdata$Clinical<-1
Pamper.pdata$EC<-0
Pamper.pdata$EC[which(Pamper.pdata$PAMPID %in% biomarker.pdata$PAID)]<- 1
Pamper.pdata$Cytokines<-0
Pamper.pdata$Cytokines[which(Pamper.pdata$PAMPID %in% Luminex.pdata$PAID)]<- 1
Pamper.pdata$Metabolome<-0
Pamper.pdata$Metabolome[which(Pamper.pdata$PAMPID %in% PAMPer.Metabolome$PAID)]<- 1
Pamper.pdata$Lipidome<- Pamper.pdata$Metabolome

Pamper.pdata$Proteome<-0
Pamper.pdata$Proteome[which(Pamper.pdata$PAMPID %in% PAMPer.Proteome$PAID)]<- 1

Pamper.pdata$All5_layer<- 0
Pamper.pdata$All5_layer[which(Pamper.pdata$PAMPID %in% Layer5_merge$PAID) ]<- 1

Pamper.pdata$All6_layer<- 0
Pamper.pdata$All6_layer[which(Pamper.pdata$All5_layer==1 & Pamper.pdata$Proteome==1)]<- 1
###
Pamper.pdata$iss[is.na(Pamper.pdata$iss)]<-median(Pamper.pdata$iss,na.rm = T)
Pamper.pdata$Severity<-Pamper.pdata$iss %>% cut(c(-1,14,24,100),labels = c("Mild","Mod","Severe"))
table(is.na(Pamper.pdata$Severity))

###
Var<-c("Severity","Clinical","EC","Cytokines","Lipidome","Metabolome","Proteome","All5_layer","All6_layer")
Matrix<-Pamper.pdata[,Var]

M2<-Matrix %>% group_by(Severity) %>% summarise_each(sum) %>% pivot_longer(!Severity,names_to="Layers",values_to="Nubmers")
M2$Layers<-factor(M2$Layers,levels = c("Clinical","EC","Cytokines","Lipidome","Metabolome","Proteome","All5_layer","All6_layer"))
# Stacked barplot with multiple groups
ggplot(data=M2, aes(x=Layers, y=Nubmers, fill=Severity)) + geom_bar(stat="identity")+ggplot2::theme(legend.position = "bottom")+theme(aspect.ratio=1)+
  geom_text(aes(y=Nubmers, label=Nubmers), vjust=1, 
            color="black", size=3.5)
###

dim(MM.pdata2)
table(Pamper.pdata$Severity)
head(M2)
