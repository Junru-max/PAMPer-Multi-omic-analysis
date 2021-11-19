library(Seurat)
library(dplyr)
library(tidyr)
library(FELLA)
library(org.Mm.eg.db)
library(KEGGREST)
library(igraph)
library(magrittr)
set.seed(1)
##
table(duplicated(Meta.pdata$COMP.ID))
###For pathway for contrast metabolites
Idents(PAMPer.Metabolome)<-PAMPer.Metabolome$outcome_time
table(PAMPer.Metabolome$outcome_time)
DE.lipdi.0h<-FindMarkers(PAMPer.Metabolome,ident.1 = "Early-Nonsurvivors_0",ident.2 = c("Resolving_0","Non-resolving_0"),test.use = "LR",logfc.threshold = 0,latent.vars = c("GENDER","AGE","TREATMENT"))
DE.lipdi<-DE.lipdi.0h[which(DE.lipdi.0h$p_val<0.05 & DE.lipdi.0h$avg_logFC>0),]
Meta.pdata.DE<-Meta.pdata[which(Meta.pdata$COMP.ID %in% rownames(DE.lipdi)),]
DE.Meta<-as.character(Meta.pdata$KEGG[which(Meta.pdata$module==1)])

table(rownames(DE.lipdi.0h) %in% Meta.pdata$COMP.ID)
####
 # Filter overview pathways
graph <- buildGraphFromKEGGREST( organism = "hsa",
     filter.path = c("01100", "01200", "01210", "01212", "01230"))

tmpdir <- paste0(tempdir(), "/my_database")
# Make sure the database does not exist from a former vignette build
# Otherwise the vignette will rise an error
# because FELLA will not overwrite an existing database
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
  keggdata.graph = graph,
   databaseDir = tmpdir,internalDir = FALSE,
   matrices = "diffusion",
   normality = "diffusion",
  niter = 50)
fella.data <- loadKEGGdata(
   databaseDir = tmpdir,
   internalDir = FALSE,
   loadMatrix = "diffusion")
fella.data  
cat(getInfo(fella.data))
##################For pathway enriched in 7 modules
DE.Meta<-list()
for (i in 1:7) {DE.Meta[[i]]<-as.character(Meta.pdata$KEGG[which(Meta.pdata$module==i)])
analysis.pamper <- defineCompounds(
  compounds = DE.Meta[[i]],
  data = fella.data)
getInput(analysis.pamper)
getExcluded(analysis.pamper)

analysis.pamper <- runDiffusion(
  object = analysis.pamper,
  data = fella.data,
  approx = "normality")
analysis.pamper

nlimit <- 200
vertex.label.cex <- 0.5
plot(
  analysis.pamper,
  method = "diffusion",
  data = fella.data,
  nlimit = nlimit,
  vertex.label.cex = vertex.label.cex)

tmpfile <- tempfile()
exportResults(format = "csv",
              file = paste(i,"M_Enriched_Pathway.csv",sep = ""),
              method = "diffusion",
              object = analysis.pamper,
              data = fella.data,nlimit=1000,
              LabelLengthAtPlot=100)

  
}
###