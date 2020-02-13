library("DoubletFinder",lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library("modes",lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library("Seurat")

args <- commandArgs(trailingOnly = TRUE)

rds <- as.character(args[1])
sample = basename(rds)
rdata = as.character(args[2])
load(rdata)


doublets <-function(dfso){
  dfso <- NormalizeData(dfso, verbose = F)
  dfso <- ScaleData(dfso,verbose = F)
  dfso <- FindVariableFeatures(dfso, do.plot=FALSE)
  dfso <- RunPCA(dfso, pc.genes = dfso@var.genes, pcs.print = 0,verbose = F,npcs =30)
  npcs = 10
  dfso <- RunUMAP(dfso, verbose=TRUE,dims = 1:npcs)
  
	  
  sweep.res.list_kidney <- paramSweep_v3(dfso,PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
#  sweep.res.list_kidney <- paramSweep_v3(dfso, PCs = 1:10, sct = FALSE)
#  gt.calls <- dfso@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
#  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
#  bcmvn_kidney <- find.pK(sweep.stats_kidney)


  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(dfso$annot)
  perc = 0.008 * (length(colnames(dfso))/1000)
  nExp_poi <- round(perc*length(colnames(dfso)))#dfso@cell.names
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE,PCs = 1:10)
  pAAN=tail(names(dfso@meta.data),2)[1]
  dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pAAN,PCs = 1:10)
 
  ## Plot results --------------------------------------------------------------------------------------------------------------
  DF1=tail(names(dfso@meta.data),2)[1]
  DF2=tail(names(dfso@meta.data),2)[2]
  
  dfso@meta.data[,"DF_hi.lo"] <- dfso[[DF1]]
  dfso@meta.data$DF_hi.lo[which(dfso@meta.data$DF_hi.lo == "Doublet" & dfso@meta.data[[DF2]] == "Singlet")] <- "Doublet_lo"
  dfso@meta.data$DF_hi.lo[which(dfso@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
  return(dfso$DF_hi.lo)
}

#so = readRDS(file)
so$annot = so$HPCA 
doublets = doublets(so)
so$DF_hi.lo = doublets
#saveRDS(so,"x2.rds")
so=SubsetData(so,cells=names(so$DF_hi.lo)[so$DF_hi.lo =="Singlet"])
saveRDS(so,rds)
save(list=c("sce_BC","sce_filt","so"),file = rdata)


