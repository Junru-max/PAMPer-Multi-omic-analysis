library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(stringr)
library(dplyr)
library(missRanger)
library(ggplot2)
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

###Lipidomics
PAMPer.Lipidome<-readRDS(file = "Rdata/PAMPer.lipidome.rds")
PAMPer.Lipidome$PAID_TIME<-PAMPer.Lipidome$`PAMPER ID NUMBER_time`
###
PAMPer.Metabolome$PAID_TIME<-PAMPer.Metabolome$PAMPER.ID.NUMBER_time
Common_PAID_TIME<-intersect(PAMPer.Metabolome$PAID_TIME,PAMPer.Lipidome$PAID_TIME)
###
Idents(PAMPer.Metabolome)<-PAMPer.Metabolome$PAID_TIME
Idents(PAMPer.Lipidome)<-PAMPer.Lipidome$PAID_TIME

Layer3_merge<-subset(PAMPer.Metabolome,idents=Common_PAID_TIME)

PAMPer.Lipidome.data<-GetAssayData(subset(PAMPer.Lipidome,idents=Common_PAID_TIME),slot = "counts")
colnames(PAMPer.Lipidome.data)<-colnames(Layer3_merge)

Layer3_merge[["Lipidome"]]<-CreateAssayObject(counts =PAMPer.Lipidome.data )
Layer3_merge$LRS<-PAMPer.Lipidome$Lipid.sig.score[match(Layer3_merge$PAID_TIME,PAMPer.Lipidome$PAID_TIME)]
Layer3_merge$Lipid.conc<-PAMPer.Lipidome$Total.conc[match(Layer3_merge$PAID_TIME,PAMPer.Lipidome$PAID_TIME)]

###

DefaultAssay(Layer3_merge) <- 'Metabolome'

VariableFeatures(Layer3_merge) <- rownames(Layer3_merge[["Metabolome"]])
Layer3_merge <-  Layer3_merge %>% ScaleData() %>% RunPCA(reduction.name = 'pca',npcs = 30)


DefaultAssay(Layer3_merge) <- 'Lipidome'

VariableFeatures(Layer3_merge) <- rownames(Layer3_merge[["Lipidome"]])
Layer3_merge <-  Layer3_merge %>% ScaleData() %>% RunPCA(reduction.name = 'apca')

Layer3_merge <- FindMultiModalNeighbors(
  Layer3_merge, reduction.list = list("pca", "apca"),k.nn = 10, 
  dims.list = list(1:20, 1:20), modality.weight.name = "Metabolome.weight"
)

Layer3_merge <- RunUMAP(Layer3_merge, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Layer3_merge <- FindClusters(Layer3_merge, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
p1 <- DimPlot(Layer3_merge, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
p2 <- DimPlot(Layer3_merge, reduction = 'wnn.umap', group.by = 'outcome_time', label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
p1 + p2
DimPlot(Layer3_merge, reduction = "wnn.umap", label = F,pt.size = 1,ncol = 4,group.by = "outcome_time",split.by = "outcome_time")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)

VlnPlot(Layer3_merge, features = "Metabolome.weight", group.by = 'outcome_time', sort = TRUE, pt.size = 0.1) +
  NoLegend()
VlnPlot(Layer3_merge, features = "Metabolome.weight", group.by = 'outcome_time', sort = F, pt.size = 0.1) +
  NoLegend()

hist(Layer3_merge$Metabolome.weight,breaks = 100)
cor(Layer3_merge$Metabolome.weight,Layer3_merge$INJURY.SEVERITY.SCORE.1)
table(Layer3_merge$outcome_time[which(Layer3_merge$Metabolome.weight<0.3)])

