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
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

###K-means clustering for 7 modules
PAMPer.Meta.heatmap<-FetchData(PAMPer.Metabolome,slot = "scale.data",vars = c(rownames(PAMPer.Metabolome),"outcome_time"))
PAMPer.Meta.heatmap2<-PAMPer.Meta.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.Meta.heatmap2$outcome_time
PAMPer.Meta.heatmap2<-PAMPer.Meta.heatmap2[,-1]
PAMPer.Meta.heatmap2<-t(PAMPer.Meta.heatmap2)
colnames(PAMPer.Meta.heatmap2)<-timepoint


col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))
ht<-Heatmap(PAMPer.Meta.heatmap2,name = "Z-score",col =col_fun,
            show_row_names = F,cluster_columns = F,width = unit(10, "cm"),height = unit(10, "cm"),
            row_km  = 7, row_km_repeats  = 1000)
ht
row_order<-row_order(ht)
row_order <- row_order[order(names(row_order))]

###For 7moudle score in PAMPer
module<-list()
Meta.pdata$module<-0

module.data<-FetchData(PAMPer.Metabolome,vars = c("outcome_time","TIME.POINT","outcome","PAID","SAMPLE.NAME"))
multimodule<- function(df, n) {
  
  varname <- paste("module", n , sep=".")
  
  mutate(df, !!varname := apply(PAMPer.Metabolome@assays$Metabolome@scale.data[module[[i]],],2,mean))
  
}
for (i in 1:7) {module[[i]]<-rownames(PAMPer.Meta.heatmap2)[row_order[[i]]]
module.data<-module.data %>% multimodule(n = i)
PAMPer.Metabolome<-AddMetaData(PAMPer.Metabolome,metadata = module.data[,paste0("module.",i)],col.name = paste0("module.",i))
Meta.pdata$module[which(Meta.pdata$COMP.ID %in% module[[i]])]<- i
}

###Heatmap for 7 module score
select.var<-paste0("module.",1:7)
PAMPer.Module.heatmap<-FetchData(PAMPer.Metabolome,slot = "scale.data",vars = c(select.var,"outcome_time"))
PAMPer.Module.heatmap2<-PAMPer.Module.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.Module.heatmap2$outcome_time
PAMPer.Module.heatmap2<-PAMPer.Module.heatmap2[,-1]
PAMPer.Module.heatmap2<-t(PAMPer.Module.heatmap2)
colnames(PAMPer.Module.heatmap2)<-timepoint
col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))
ht<-Heatmap(PAMPer.Module.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F,width = unit(10, "cm"),height = unit(10, "cm"),
            cluster_rows = F)
ht
####Lineplot for 7 module score
####single plot
PAMPer.Module.lineplot<-FetchData(PAMPer.Metabolome,slot = "scale.data",vars = c(select.var,"TIME.POINT","outcome"))
ggline(PAMPer.Module.lineplot, x = "TIME.POINT", y = "module.1", 
       add = c("mean_sd"),size = 1,point.size = 0.2,
       color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")
####combined  plots
module.lineplot<-list()
for(i in c(1:7)){  
  module.name<-colnames(PAMPer.Module.lineplot)[i]
  p.tmp<-ggline(module.data, x = "TIME.POINT", y = module.name, 
                add = c("mean_sd"),size = 1,point.size = 0.01,
                color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")
  module.lineplot[[module.name]]<-p.tmp
}
ggarrange(plotlist = module.lineplot, 
          ncol = 7, nrow = 1,common.legend = T,align = "v")

#######statistical analysis For module scores across timepoint and outcome groups
test.var<-"module.1"
test.var.sym<-rlang::sym(test.var)
#######Comprasion in 0,24,72h for 2 groups of outcome
###table prepration and assumption test
PAMPer.Module.lineplot$outcome_TIME.POINT<-paste(PAMPer.Module.lineplot$outcome,PAMPer.Module.lineplot$TIME.POINT,sep = "_")
PAMPer.Module.lineplot.ex.nons<-PAMPer.Module.lineplot[which(PAMPer.Module.lineplot$outcome %in% c("Resolving","Non-resolving")),] %>% droplevels()
PAMPer.Module.lineplot.ex.nons %>%
  group_by(outcome,TIME.POINT) %>%
  shapiro_test(test.var.sym)
ggqqplot(PAMPer.Module.lineplot.ex.nons, test.var, ggtheme = theme_bw()) +
  facet_grid(TIME.POINT ~ outcome, labeller = "label_both")
#####2-way anova
PAMPer.Module.lineplot.ex.nons<-tibble(PAMPer.Module.lineplot.ex.nons)
res.aov <-PAMPer.Module.lineplot.ex.nons %>% anova_test(
  get(test.var) ~ outcome*TIME.POINT )
get_anova_table(res.aov)
model <- lm(get(test.var) ~ TIME.POINT * outcome, data = PAMPer.Module.lineplot.ex.nons)
PAMPer.Module.lineplot.ex.nons %>%
  group_by(TIME.POINT) %>%
  anova_test(get(test.var) ~ outcome, error = model)
# pairwise comparisons
pwc <- PAMPer.Module.lineplot.ex.nons %>% 
  group_by(TIME.POINT) %>%
  emmeans_test(get(test.var) ~ outcome, p.adjust.method = "BH") 
pwc

#######Comprasion in 0h for 3 groups of outcome
####kruskal.test
PAMPer.Module.lineplot.ex.nons<-PAMPer.Module.lineplot[which(PAMPer.Module.lineplot$outcome_TIME.POINT %in% c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0")),] %>% droplevels()
kruskal.test(get(test.var) ~ outcome_TIME.POINT,data = PAMPer.Module.lineplot.ex.nons)
### Dunn test
dunnTest(get(test.var) ~ outcome_TIME.POINT,
         data=PAMPer.Module.lineplot.ex.nons,
         method="bh")    # Can adjust p-values;
###export data
saveRDS(Meta.pdata,file = "Rdata/Meta.pdata.rds")
write.csv(Meta.pdata,file = "Rdata/Meta.pdata.csv")

