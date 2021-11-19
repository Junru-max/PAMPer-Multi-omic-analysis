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
####Global funtion
##factor to numberic
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
###detection for coloumn containing miss value
na.detect<-function(object){ length.tmp<-ncol(object)
       na.number<-list()
  for (i in 1:length.tmp) {
    na.number[i]<-nrow(object)-table(is.na(object[,i]))[1]
    
  }
       na.number<-unlist(na.number)
       na.table<-data.frame(name=colnames(object),number=na.number)
       na.table<-na.table %>% filter(number!=0)
       return(na.table)
}
##meidan imputaion for missing value
median.impute<-function(object,imp.variable){
  length.tmp<-length(imp.variable)
  for (i in 1:length.tmp) {  object[,imp.variable[i]][is.na(object[,imp.variable[i]])]<-
    median(object[,imp.variable[i]],na.rm = T)
    
  }
  return(object)
}
###Plot UMAP\PCA\Heatmap
plot.combo<-function(object,Select.var){
  plot.list<-list()
  Select.var.sym<-rlang::sym(Select.var)
  myPalette<-hue_pal()(length(unique(object@meta.data[,Select.var])))
  plot.list[["umap.group"]]<-DimPlot(object, reduction = "umap", label = T,pt.size = 1,ncol = 2,group.by = Select.var)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette)
  plot.list[["umap.split"]]<-DimPlot(object, reduction = "umap", label = F,pt.size = 1,ncol = 4,group.by = Select.var,split.by = Select.var)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette)
  plot.list[["pca.group"]]<-DimPlot(object, reduction = "pca", label = T,pt.size = 1,ncol = 2,group.by = Select.var)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette)
  plot.list[["pca.split"]]<-DimPlot(object, reduction = "pca", label = F,pt.size = 1,ncol = 4,group.by = Select.var,split.by = Select.var)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette)
  PAMPer.Meta.heatmap<-FetchData(object,slot = "scale.data",vars = c(rownames(object),Select.var))
  PAMPer.Meta.heatmap2<-PAMPer.Meta.heatmap %>% group_by(!!Select.var.sym) %>% summarise_each(funs(mean)) 
  timepoint<-pull(PAMPer.Meta.heatmap2,!!Select.var.sym)
  PAMPer.Meta.heatmap2<-PAMPer.Meta.heatmap2[,-1]
  PAMPer.Meta.heatmap2<-t(PAMPer.Meta.heatmap2)
  colnames(PAMPer.Meta.heatmap2)<-timepoint
  #rownames(PAMPer.Meta.heatmap2)<-Meta.pdata$BIOCHEMICAL[match(rownames(PAMPer.Meta.heatmap2),Meta.pdata$COMP.ID)]
  col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))
  plot.list[["heatmap"]]<-Heatmap(PAMPer.Meta.heatmap2,name = "Z-score",col =col_fun,
                                  row_gap = unit(2, "mm"),
                                  show_row_names = F,cluster_columns = F,width = unit(10, "cm"),height = unit(10, "cm"),
                                  cluster_rows = T)
  return(plot.list)
}


####object claim
##Pamper.pdata.raw##Patients raw meta.data
##Pamper.pdata.imp##Patients imputed meta.data(median)

##Save object
saveRDS(Pamper.pdata.raw,file = "Rdata_review/Pamper.pdata.raw.rdata")
saveRDS(Pamper.pdata.imp,file = "Rdata_review/Pamper.pdata.imp.rdata")

##read object
##NA Detection and Imputation
Pamper.pdata.imp<-Pamper.pdata.raw
na.table<-na.detect(Pamper.pdata.imp)
var.to.imp<-c("INR","vitals_sbp","vitals_hr","icu_los","trans_icu_los","PH_time","ais_head","ais_face","ais_chest","ais_extremity","ais_external","ais_abdomen","ed_coagulopathy")
Pamper.pdata.imp<-median.impute(Pamper.pdata.imp,imp.variable = var.to.imp)
Pamper.pdata.imp$ed_coagulopathy_grade[is.na(Pamper.pdata.imp$ed_coagulopathy_grade)]<-0
Pamper.pdata.imp$Arms<-factor(Pamper.pdata.imp$plasma_on_helicopter,levels = c(2,1),labels = c("Standard","Plasma"))
##    
#Fig1CDE Fig2ABC sFig2 AB

####pdata only in metabolomics 
pdata<-read.csv(file = "rawdata/PAMPer_Meta_pdata.csv",row.names = 1) %>% as.data.frame() 
pdata<-pdata %>% mutate(ID=rownames(pdata)) %>% filter(TREATMENT !="") %>% drop.levels()
level_key <- c("N/A" = "Baseline")
for (i in c("TREATMENT","TRAUMATIC.BRAIN.INJURY","ALIVE.AT.30.DAYS","TIME.POINT")) {
  pdata[,i]<-pdata[,i] %>% recode_factor(!!!level_key)
}

pdata$PAMPER.ID.NUMBER<-as.character(pdata$PAMPER.ID.NUMBER)
pdata$PAMPER.ID.NUMBER[which(pdata$TREATMENT=="Baseline")]<-as.character(pdata$CLIENT.SAMPLE.ID[which(pdata$TREATMENT=="Baseline")])
pdata$outcome<-Pamper.pdata$outcome[match(pdata$PAMPER.ID.NUMBER,Pamper.pdata$PAMPID)] %>% as.character()
pdata$outcome[which(pdata$TREATMENT=="Baseline")]<-"Baseline"
pdata$outcome<-factor(pdata$outcome,levels = c("Baseline","Resolving","Non-resolving","Early-Nonsurvivors"))
pdata$coagulopathy<-Pamper.pdata$ed_coagulopathy[match(pdata$PAMPER.ID.NUMBER,Pamper.pdata$PAMPID)] %>% as.character()
pdata$TIME_POINT<-pdata$TIME.POINT
pdata$PAID<-pdata$PAMPER.ID.NUMBER
pdata$AGEGROUP<-pdata$AGE  %>% cut(breaks=c(0,40,60,200),labels = c("young","middle","old")) 
pdata$Severity<-pdata$INJURY.SEVERITY.SCORE.1  %>% cut(breaks=c(0,14,24,75),labels = c("Mild","Mod","Severe")) %>% as.character()
pdata$Severity[is.na(pdata$Severity)]<-"Baseline"
pdata$Severity<-factor(pdata$Severity,levels = c("Baseline","Mild","Mod","Severe"))
pdata$TREATMENT<-pdata$TREATMENT %>% fct_relevel("Baseline","standard care","prehospital plasma")
pdata$ALIVE.AT.30.DAYS<-pdata$ALIVE.AT.30.DAYS %>% fct_relevel("Baseline","yes","no")
for (i in c("TREATMENT","TRAUMATIC.BRAIN.INJURY","ALIVE.AT.30.DAYS","TIME.POINT","Severity","AGEGROUP","outcome","PAMPER.ID.NUMBER","GENDER","coagulopathy")) {
  pdata[,paste(i,"time",sep = "_")]=interaction(pdata[,i],pdata[,"TIME.POINT"],sep = "_",lex.order = F,drop = T) 
  
}
Filter_ID<-names(which(table(pdata[which(pdata$outcome %in% c("Non-resolving")),"PAMPER ID NUMBER"])<3))
pdata<-pdata %>% filter(!outcome_time %in% c("Early-Nonsurvivors_24","Early-Nonsurvivors_72")) %>% filter(!duplicated(PAMPER.ID.NUMBER_time)) 
pdata$outcome_time<-pdata$outcome_time %>% droplevels()
rownames(pdata)<-pdata$ID
###Annotaion for metabolites
Meta.pdata<-read.csv(file = "Rawdata/meta_pdata.csv",header = T)
Meta.pdata$COMP.ID<-as.character(Meta.pdata$COMP.ID)
######Object create and dimension reduction 
Impdata<-read.csv(file = "rawdata/PAMPer.Metabolimics1.csv",header = T,row.names = 1)
pdata<-pdata[intersect(rownames(Impdata),rownames(pdata)),]
Impdata<-t(Impdata[rownames(pdata),])
Impdata[1:5,1:5]
rownames(Impdata)<-Meta.pdata$COMP.ID

PAMPer.Metabolome<-CreateSeuratObject(counts = Impdata,meta.data = pdata,assay="Metabolome")
PAMPer.Metabolome <- FindVariableFeatures(PAMPer.Metabolome, selection.method = "vst", nfeatures = nrow(PAMPer.Metabolome))
PAMPer.Metabolome <- ScaleData(PAMPer.Metabolome,do.scale = T,do.center = T)
PAMPer.Metabolome <- RunPCA(PAMPer.Metabolome, npcs = 30, verbose = FALSE)
PAMPer.Metabolome <- RunUMAP(PAMPer.Metabolome, reduction = "pca", dims = 1:20)

PAPMPer.Meta.plot<-plot.combo(PAMPer.Metabolome,"outcome_time")
####
saveRDS(PAMPer.Metabolome,file = "Rdata/PAMPer.Metabolome.rds")
