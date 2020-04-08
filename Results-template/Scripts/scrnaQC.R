#.libPaths(c("/data//CCBR_Pipeliner/db/PipeDB/scrna_lib/R-3.6.1/library",.libPaths()[1],.libPaths()[2],.libPaths()[3]))
#.libPaths(c("/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA",.libPaths()[1],.libPaths()[2],.libPaths()[3]))
.libPaths("/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.1_scRNA")
print(.libPaths())

args <- commandArgs(trailingOnly = TRUE)

file <- as.character(args[1])
sample = basename(file)
outRDS = as.character(args[2])
outRdata = as.character(args[3])
species = as.character(args[4])
resolution = as.character(args[5])
clusterAlg =  as.numeric(args[6])
annotDB = as.character(args[7])
citeseq = as.character(args[8])
#matrix = paste0(file,"/outs/filtered_feature_bc_matrix.h5")

matrix=file

resolution = as.numeric(strsplit(gsub(",+",",",resolution),split=",")[[1]]) #remove excess commas, split into numeric vector

#library(future)
#plan("multiprocess", workers = 8) 


#library(DropletUtils,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
#library(scran,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
#library(scater,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library(BiocGenerics)
library(Biobase)#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library(farver)#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library(S4Vectors)#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library("Seurat")#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library("DoubletFinder")#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library(AnnotationDbi)
library("modes")#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")

library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library("URD")#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library(Routliers)#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")
library(dplyr)
library(Matrix) 
library(reshape2)
library(tools)
library("DoubletFinder")#,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/")



###Run Seurat Clustering 
seuratClustering = function(so){
	
	tail(strsplit(file,"/")[[1]],1)
	so$Sample=gsub(".h5","",tail(strsplit(file,"/")[[1]],1))
	so = SCTransform(so)
	
	so <- RunPCA(object = so, features = VariableFeatures(object = so), do.print = TRUE, pcs.print = 1:5,genes.print = 0,verbose=F,npcs = 30)
	
	if(ncol(so)<50000) {
		mat1 <- t(as.matrix(FetchData(object = so,slot = "counts",vars = rownames(so))))
		urdObj <- createURD(count.data = mat1, min.cells=3, min.counts=3)
		varGenes = so@assays$RNA@var.features[so@assays$RNA@var.features %in% rownames(urdObj@logupx.data)]
		pcs <- calcPCA(urdObj, mp.factor = 1,pcs.store = 30,genes.use=varGenes)
		npcs =  sum(pcs@pca.sig)
	}
	else {
		npcs = 30
	}
	
	so <- FindNeighbors(so,dims = 1:npcs)
	for (i in 1:length(resolution)){
		so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = resolution[i],algorithm = clusterAlg)
	}
	so <- RunUMAP(so,dims = 1:npcs)
	
	return(so)
}

doublets <-function(dfso){
	dfso <- SCTransform(dfso)
	dfso <- RunPCA(dfso, pc.genes = dfso@var.genes, pcs.print = 0,verbose = F,npcs =30)
	npcs = 10
	dfso <- RunUMAP(dfso, verbose=TRUE,dims = 1:npcs)
	
			
	sweep.res.list_kidney <- paramSweep_v3(dfso,PCs = 1:10, sct = T)
	sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
	#print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
	bcmvn_kidney <- find.pK(sweep.stats_kidney)
	## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
	#sweep.res.list_kidney <- paramSweep_v3(dfso, PCs = 1:10, sct = FALSE)
	#gt.calls <- dfso@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
	#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
	#bcmvn_kidney <- find.pK(sweep.stats_kidney)

	## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
	homotypic.prop <- modelHomotypic(dfso$annot)
	perc = 0.008 * (length(colnames(dfso))/1000)
	nExp_poi <- round(perc*length(colnames(dfso)))#dfso@cell.names
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
	## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
	dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE,PCs = 1:10,sct = T)
	pAAN=tail(names(dfso@meta.data),2)[1]
	dfso <- doubletFinder_v3(dfso, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pAAN,PCs = 1:10,sct = T)
 
  ## Plot results --------------------------------------------------------------------------------------------------------------
	#DF1=tail(names(dfso@meta.data),2)[1]
	#DF2=tail(names(dfso@meta.data),2)[2]
	#dfso@meta.data[,"DF_hi.lo"] <- dfso[[DF2]]
	#dfso@meta.data$DF_hi.lo[which(dfso@meta.data$DF_hi.lo == "Doublet" & dfso@meta.data[[DF2]] == "Singlet")] <- "Doublet_lo"
	#dfso@meta.data$DF_hi.lo[which(dfso@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
	return(dfso)
}



convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}

if (citeseq=="Yes"){
	fileInput = Read10X_h5(matrix)
	so_BC = CreateSeuratObject(fileInput[[1]])
	so_BC[['CITESeq']] = CreateAssayObject(counts=fileInput[[2]])
	so_BC = NormalizeData(so_BC,assay="CITESeq",normalization.method="CLR")
	so_BC = ScaleData(so_BC,assay="CITESeq")
}else{
	so_BC = Read10X_h5(matrix) 
	so_BC = CreateSeuratObject(so_BC)
}

if(species=="human"){so_BC[["percent.mt"]] <- PercentageFeatureSet(so_BC, pattern = "^MT-")}
if(species=="mouse"){so_BC[["percent.mt"]] <- PercentageFeatureSet(so_BC, pattern = "^mt-")}

nCount_out = outliers_mad(so_BC$nCount_RNA,threshold = 3)$LL_CI_MAD
nFeature_out = outliers_mad(so_BC$nFeature_RNA,threshold = 3)$LL_CI_MAD
mt_out = outliers_mad(so_BC$percent.mt,threshold = 3)$UL_CI_MAD

so_filt <- subset(so_BC, subset = nFeature_RNA > nFeature_out & nCount_RNA > nFeature_out & percent.mt < mt_out)

load("/data/CCBR_Pipeliner/db/PipeDB/scrna_lib/immgenMod.RData")

#rownames(sce_norm) = gsub("_","-",rownames(sce_norm))
#rownames(sce_norm) = make.names(sce_norm@rowRanges@elementMetadata$Symbol,unique = T)
so = seuratClustering(so_filt)

if(species=="human"){
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
}

if(species=="mouse"){
	s.genes <- convertHumanGeneList(cc.genes$s.genes)
	g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
}


so = CellCycleScoring(so,s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

runSingleR = function(obj,refFile,fineORmain){
	sce = as.SingleCellExperiment(obj,assay = "SCT")
	ref = refFile
	s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]])
	#return(s$pruned.labels)
	return(s)
	print(head(s$pruned.labels))
}

if(species == "human"){
	so@misc$SingleR$HPCA_main <- runSingleR(so,HumanPrimaryCellAtlasData(),"label.main")
	so$HPCA_main = so@misc$SingleR$HPCA_main$pruned.labels
	so@misc$SingleR$HPCA <-  runSingleR(so,HumanPrimaryCellAtlasData(),"label.fine")
	so$HPCA = so@misc$SingleR$HPCA$pruned.labels
	so@misc$SingleR$BP_encode_main <-  runSingleR(so,BlueprintEncodeData(),"label.main")
	so$BP_encode_main=so@misc$SingleR$BP_encode_main$pruned.labels
	so@misc$SingleR$BP_encode <-  runSingleR(so,BlueprintEncodeData(),"label.fine")
	so$BP_encode=so@misc$SingleR$BP_encode$pruned.labels
	so@misc$SingleR$monaco_main <-  runSingleR(so,MonacoImmuneData(),"label.main")
	so$monaco_main=so@misc$SingleR$monaco_main$pruned.labels
	so@misc$SingleR$monaco <-	 runSingleR(so,MonacoImmuneData(),"label.fine")
	so$monaco = so@misc$SingleR$monaco$pruned.labels
	so@misc$SingleR$dice_main <-  runSingleR(so,DatabaseImmuneCellExpressionData(),"label.main")
	so$dice_main = so@misc$SingleR$dice_main$pruned.labels
	so@misc$SingleR$dice <- runSingleR(so,DatabaseImmuneCellExpressionData(),"label.fine")
	so$dice = so@misc$SingleR$dice$pruned.labels
	so@misc$SingleR$Novershtern_main <-  runSingleR(so,NovershternHematopoieticData(),"label.main")
	so$Novershtern_main = so@misc$SingleR$Novershtern_main$pruned.labels
	so@misc$SingleR$Novershtern <- runSingleR(so,NovershternHematopoieticData(),"label.fine")
	so$Novershtern = so@misc$SingleR$Novershtern$pruned.labels
}

if(species == "mouse"){
	so@misc$SingleR$immgen_main <-  runSingleR(so,ImmGenData(),"label.main")
	so$immgen_main = so@misc$SingleR$immgen_main$pruned.labels
	so@misc$SingleR$immgen <- runSingleR(so,ImmGenData(),"label.fine")
	so$immgen = so@misc$SingleR$immgen$pruned.labels
	so@misc$SingleR$mouseRNAseq_main <-  runSingleR(so,MouseRNAseqData(),"label.main")
	so$mouseRNAseq_main = so@misc$SingleR$mouseRNAseq_main$pruned.labels
	so@misc$SingleR$mouseRNAseq <- runSingleR(so,MouseRNAseqData(),"label.fine")
	so$mouseRNAseq = so@misc$SingleR$mouseRNAseq$pruned.labels
}

so$annot = so[[paste0(annotDB,"_main")]] 

so = doublets(so)
so$DF_hi.lo = so[[tail(names(so@meta.data),1)]]
so_noDoublet=subset(so,cells=names(so$DF_hi.lo)[so$DF_hi.lo =="Singlet"])

saveRDS(so_noDoublet,outRDS)
save.image(outRdata)
