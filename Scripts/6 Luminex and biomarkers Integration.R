library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(stringr)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(dendsort)

####Luminex_pdata
PAMPer.luminex<-read.csv(file = "rawdata/PAMPer luminex .csv",row.names = 1)
dim(PAMPer.luminex)
PAMPer.luminex[is.na(PAMPer.luminex)]<-0
PAMPer.luminex<-PAMPer.luminex[-which(apply(PAMPer.luminex, 1, sum)==0),]
timepoint<-rownames(PAMPer.luminex) %>% str_extract("(?<=\\-).+?(?=\\-)")
timepoint<- timepoint %>% gsub(pattern = "\\b0H\\b",replacement = "0HR") %>% gsub(pattern = "hr",replacement = "HR") %>% gsub(pattern = "HR",replacement = "") 
PAID<-rownames(PAMPer.luminex) %>% str_extract(".+?\\-") %>% str_extract(".+?\\-") %>% str_sub(5,-2)
Luminex.pdata<-data.frame(rownames =rownames(PAMPer.luminex),PAID=PAID, timepoint=timepoint,PAID_TIME=paste(PAID,timepoint,sep = "_"))
Luminex.pdata$outcome<-Pamper.pdata$outcome[match(Luminex.pdata$PAID,Pamper.pdata$PAMPID)]
Luminex.pdata$outcome_time<-paste(Luminex.pdata$outcome,Luminex.pdata$timepoint,sep = "_")
##biomarker_pdata
PAMper.EC.marker<-read.csv(file = "rawdata/Data syndecan + thrombomodulin_PAMper.csv")
PAMper.Adiponectin.marker<-read.csv(file = "rawdata/data_Adiponectin_suPAR_Cell death_s100A10_PAMper.csv")
PAMper.VEGF<-read.csv(file = "rawdata/data_VEGF_R1_PAMper.csv")
PAMPer.biomarker<-inner_join(PAMper.EC.marker,PAMper.Adiponectin.marker,by="ï..Sample.ID") %>% inner_join(PAMper.VEGF,by="ï..Sample.ID")
PAMPer.biomarker<-PAMPer.biomarker[!duplicated(PAMPer.biomarker$ï..Sample.ID),]
timepoint<-PAMPer.biomarker$ï..Sample.ID %>% str_extract("(?<=\\-).+?(?=\\-)")
timepoint<- timepoint %>% gsub(pattern = "\\b0H\\b",replacement = "0HR") %>% gsub(pattern = "hr",replacement = "HR") %>% gsub(pattern = "Hr",replacement = "HR") %>% gsub(pattern = "27",replacement = "24") %>% gsub(pattern = "HR",replacement = "")
PAID<-PAMPer.biomarker$ï..Sample.ID %>% str_extract(".+?\\-") %>% str_extract(".+?\\-") %>% str_sub(5,-2)
PAMPer.biomarker<-PAMPer.biomarker[,-c(3,5,7,9,11,13,15,16,17)]

PAMPer.biomarker$Thrombomodulin_ng.ml<-as.numeric.factor(PAMPer.biomarker$Thrombomodulin_ng.ml)
PAMPer.biomarker$S100A10_ng.mL<-as.numeric.factor(PAMPer.biomarker$S100A10_ng.mL)
PAMPer.biomarker$suPAR_ng.mL<-as.numeric.factor(PAMPer.biomarker$suPAR_ng.mL)
PAMPer.biomarker$Cell.death_<-as.numeric.factor(PAMPer.biomarker$Cell.death_)
PAMPer.biomarker$Adiponectin_ng.mL<-as.numeric.factor(PAMPer.biomarker$Adiponectin_ng.mL)
PAMPer.biomarker$sFLT.1.sVEGFR1_pg.mL<-as.numeric.factor(PAMPer.biomarker$sFLT.1.sVEGFR1_pg.mL)
PAMPer.biomarker[is.na(PAMPer.biomarker)]<-0
rownames(PAMPer.biomarker)<-PAMPer.biomarker$ï..Sample.ID
biomarker.pdata<-data.frame(rownames=PAMPer.biomarker$ï..Sample.ID,PAID=PAID, timepoint=timepoint)
biomarker.pdata$PAID_TIME<-paste(biomarker.pdata$PAID,biomarker.pdata$timepoint,sep = "_")
biomarker.pdata<-filter(biomarker.pdata,!is.na(timepoint))
biomarker.pdata$outcome<-Pamper.pdata$outcome[match(biomarker.pdata$PAID,Pamper.pdata$PAMPID)]
biomarker.pdata$outcome_time<-paste(biomarker.pdata$outcome,biomarker.pdata$timepoint,sep = "_")
dim(biomarker.pdata)

####Layer4_merge
####
Idents(Layer3_merge)<-Layer3_merge$outcome
Layer4_merge<-subset(Layer3_merge,idents="Baseline",invert=T)

###Common_pdata
Common_PAID_TIME<-Layer4_merge$PAID_TIME %>% intersect(Luminex.pdata$PAID_TIME)
###Megrge_data for common samples
Idents(Layer4_merge)<-Layer4_merge$PAID_TIME
Layer4_merge<-subset(Layer4_merge,idents=Common_PAID_TIME)

###Luminex exploration
sample_ID<-Luminex.pdata$rownames[match(Layer4_merge$PAID_TIME,Luminex.pdata$PAID_TIME)] %>% as.character()
table(duplicated(sample_ID))
PAMPer.luminex_4<-t(PAMPer.luminex[sample_ID,])
colnames(PAMPer.luminex_4)<-colnames(Layer4_merge)
Layer4_merge[["Luminex"]]<-CreateAssayObject(counts = PAMPer.luminex_4)
DefaultAssay(Layer4_merge)<-"Luminex"
Layer4_merge <- SetAssayData(Layer4_merge,slot = "data", new.data = log2(PAMPer.luminex_4+1))
Layer4_merge <- FindVariableFeatures(Layer4_merge, selection.method = "vst", nfeatures = nrow(Layer4_merge))
Layer4_merge <- ScaleData(Layer4_merge,do.scale = T,do.center = T)

###
###heatmap
DefaultAssay(Layer4_merge)<-"Luminex"
Luminex.heatmap<-FetchData(Layer4_merge,slot = "scale.data",vars = c(rownames(Layer4_merge),"outcome_time"))
Luminex.heatmap2<-Luminex.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-Luminex.heatmap2$outcome_time
Luminex.heatmap2<-Luminex.heatmap2[,-1]
Luminex.heatmap2<-t(Luminex.heatmap2)
colnames(Luminex.heatmap2)<-timepoint
col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(Luminex.heatmap2,name = "Exp",col =col_fun,
            row_split = 2, show_row_names = T,cluster_columns = F,
            width = unit(10, "cm"),height = unit(10, "cm"))
ht
######

#####Layer5_merge
####
Idents(Layer3_merge)<-Layer3_merge$outcome
Layer5_merge<-subset(Layer3_merge,idents="Baseline",invert=T)

###Common_pdata
Common_PAID_TIME<-Layer5_merge$PAID_TIME %>% intersect(Luminex.pdata$PAID_TIME) %>% intersect(biomarker.pdata$PAID_TIME)
###Megrge_data for common samples
Idents(Layer5_merge)<-Layer5_merge$PAID_TIME
Layer5_merge<-subset(Layer5_merge,idents=Common_PAID_TIME)

###Luminex layer
sample_ID<-Luminex.pdata$rownames[match(Layer5_merge$PAID_TIME,Luminex.pdata$PAID_TIME)] %>% as.character()
table(duplicated(sample_ID))
PAMPer.luminex_5<-t(PAMPer.luminex[sample_ID,])
colnames(PAMPer.luminex_5)<-colnames(Layer5_merge)
Layer5_merge[["Luminex"]]<-CreateAssayObject(counts = PAMPer.luminex_5)
DefaultAssay(Layer5_merge)<-"Luminex"
Layer5_merge <- SetAssayData(Layer5_merge,slot = "data", new.data = log2(PAMPer.luminex_5+1))
Layer5_merge <- FindVariableFeatures(Layer5_merge, selection.method = "vst", nfeatures = nrow(Layer5_merge))
Layer5_merge <- ScaleData(Layer5_merge,do.scale = T,do.center = T)
###EC layer
sample_ID<-biomarker.pdata$rownames[match(Layer5_merge$PAID_TIME,biomarker.pdata$PAID_TIME)] %>% as.character()
table(duplicated(sample_ID))
colnames(PAMPer.biomarker)
PAMPer.biomarker_5<-t(PAMPer.biomarker[sample_ID,2:8])
colnames(PAMPer.biomarker_5)<-colnames(Layer5_merge)
Layer5_merge[["EC"]]<-CreateAssayObject(counts = PAMPer.biomarker_5)
DefaultAssay(Layer5_merge)<-"EC"
Layer5_merge <- SetAssayData(Layer5_merge,slot = "data", new.data = log2(PAMPer.biomarker_5+1))
Layer5_merge <- FindVariableFeatures(Layer5_merge, selection.method = "vst", nfeatures = nrow(Layer5_merge))
Layer5_merge <- ScaleData(Layer5_merge,do.scale = T,do.center = T)

###heatmap
DefaultAssay(Layer5_merge)<-"EC"
EC.heatmap<-FetchData(Layer5_merge,slot = "scale.data",vars = c(rownames(Layer5_merge),"outcome_time"))
EC.heatmap2<-EC.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-EC.heatmap2$outcome_time
EC.heatmap2<-EC.heatmap2[,-1]
EC.heatmap2<-t(EC.heatmap2)
colnames(EC.heatmap2)<-timepoint
col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(EC.heatmap2,name = "Exp",col =col_fun,
            row_split = 2, show_row_names = T,cluster_columns = F,
            width = unit(10, "cm"),height = unit(10, "cm"))
ht
######
