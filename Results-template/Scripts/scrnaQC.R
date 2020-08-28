#.libPaths(c("/data//CCBR_Pipeliner/db/PipeDB/scrna_lib/R-3.6.1/library",.libPaths()[1],.libPaths()[2],.libPaths()[3]))
#.libPaths(c("/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA",.libPaths()[1],.libPaths()[2],.libPaths()[3]))
.libPaths("/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.1_scRNA")
print(.libPaths())

args <- commandArgs(trailingOnly = TRUE)

file <- as.character(args[1])
sample = basename(file)
outRDS = as.character(args[2])
#outRdata = as.character(args[3])
outImageDir=as.character(args[3])
species = as.character(args[4])
resolution = as.character(args[5])
clusterAlg =  as.numeric(args[6])
annotDB = as.character(args[7])
citeseq = as.character(args[8])
#matrix = paste0(file,"/outs/filtered_feature_bc_matrix.h5")
outRDS.doublets=gsub("\\.rds","_doublets\\.rds",outRDS)
matrix=file
sample = gsub("\\.h5$","",sample)
resolutionString = as.character(strsplit(gsub(",+",",",resolution),split=",")[[1]])
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
library(VennDiagram)
library(ggplot2)
library(scales)
library(cluster)

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
	so@misc$npcs=npcs
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
fileInput = Read10X_h5(matrix)

if (citeseq=="Yes"){
	so_BC = CreateSeuratObject(fileInput[[1]])
	so_BC[['CITESeq']] = CreateAssayObject(counts=fileInput[[2]])
	so_BC = NormalizeData(so_BC,assay="CITESeq",normalization.method="CLR")
	so_BC = ScaleData(so_BC,assay="CITESeq")
}else{
	if (is.list(fileInput)==FALSE){
		so_BC = CreateSeuratObject(fileInput)
	}	else{
		so_BC = CreateSeuratObject(fileInput[[1]])
	}
}

if(species=="human"){so_BC[["percent.mt"]] <- PercentageFeatureSet(so_BC, pattern = "^MT-")}
if(species=="mouse"){so_BC[["percent.mt"]] <- PercentageFeatureSet(so_BC, pattern = "^mt-")}

#NW EDIT 2020 MAY: REIMPLEMENTED LOG TRANSFORM ON THRESHOLDS
nCount_out = outliers_mad(log2(so_BC$nCount_RNA),threshold = 3)$LL_CI_MAD
nFeature_out = outliers_mad(log2(so_BC$nFeature_RNA),threshold = 3)$LL_CI_MAD
mt_out = outliers_mad(log2(so_BC$percent.mt),threshold = 3)$UL_CI_MAD

cellsToRemove.Feature= colnames(so_BC)[which(log2(so_BC$nFeature_RNA)<nFeature_out)]
cellsToRemove.Count = colnames(so_BC)[which(log2(so_BC$nCount_RNA)<nCount_out)]
cellsToRemove.Mito = colnames(so_BC)[which(log2(so_BC$percent.mt)>mt_out)]

so_filt <- subset(so_BC,cells=unique(c(cellsToRemove.Feature,cellsToRemove.Count,cellsToRemove.Mito)),invert=T)

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

#### IMAGE OUTPUTS #####
#Preliminary statistics
png(paste0(outImageDir,"/filterStats_",sample,".png"),height=3000,width=3000,res=500,units="px")
par(mar=c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.05, y=0.95, paste0("QC Metrics for ",sample,":"),
     cex = 1.5, col = "black", family="serif", font=2, adj=0)

text(x=0.075, y=0.9, paste0("Cells before filtering: ", ncol(so_BC)),
     cex = 1.25, col = "black", family="serif", font=3, adj=0)
text(x=0.075, y=0.865, paste0("Cells after filtering: ", ncol(so_noDoublet)),
     cex = 1.25, col = "black", family="serif", font=3, adj=0)

text(x=0.1, y=0.815,paste0("Cells removed for nFeature (number of genes): ",length(cellsToRemove.Feature)),
		 cex = 1, col="dimgray",family="serif", font = 1, adj =0)
text(x=0.1, y=0.785,paste0("Cells removed for nCount (number of genes): ",length(cellsToRemove.Count)),
		 cex = 1, col="dimgray",family="serif", font = 1, adj =0)
text(x=0.1, y=0.755,paste0("Cells removed for mitochondrial percentage: ",length(cellsToRemove.Mito)),
		 cex = 1, col="dimgray",family="serif", font = 1, adj =0)
text(x=0.1, y=0.725,paste0("Cells removed as doublets: ",length(which(so$DF_hi.lo =="Doublet"))),
		 cex = 1, col="dimgray",family="serif", font = 1, adj =0)

text(x=0.1, y=0.68, paste0("Number of principal components used: ",so@misc$npcs),
		 cex = 1, col="dimgray",family="serif", font = 1, adj =0)

dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)
#Venn Diagram of initial filtered cells
colorList=hue_pal()(3)
venn.diagram(x=list(cellsToRemove.Feature,cellsToRemove.Count,cellsToRemove.Mito),
			 filename=paste0(outImageDir,"/cellsRemovedVenn_",sample,".png"),
			 category.names=c("Below feature threshold","Below count threshold","Above mitochondrial threshold"),
			 main = paste(sample,"filtered cells:",
						  length(unique(c(cellsToRemove.Feature,cellsToRemove.Count,cellsToRemove.Mito),sep=" "))),
			 cat.dist=c(-0.07,-0.07,-0.07),
			 fill = colorList,
			 alpha = 0.5
			)

#Violin Plots
pdf(paste0(outImageDir,"/nFeature_preFilter_",sample,".pdf"))#,height=3000,width=3000,res=500,units="px")
gplot=ggplot(so_BC@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(nFeature_RNA))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[1],labels=c("nFeature"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+geom_hline(yintercept=nFeature_out,linetype="dashed")+labs(x=paste(sample,"Pre-filter"), y="Log2(nFeature_RNA)")+
	labs(fill = paste(sample,"Pre-filter"),title="Pre-filter nFeature_RNA")
dev.off()
#Post-filter nFeature
pdf(paste0(outImageDir,"/nFeature_postFilter_",sample,".pdf"))
gplot=ggplot(so_filt@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(nFeature_RNA))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[1],labels=c("nFeature"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+labs(x=sample, y="Log2(nFeature_RNA)")+
	labs(fill = paste(sample,"Post-filter"),title="Post-filter nFeature_RNA")
dev.off()
#Pre-filter nCount
pdf(paste0(outImageDir,"/nCount_preFilter_",sample,".pdf"))
gplot=ggplot(so_BC@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(nCount_RNA))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[2],labels=c("nCount"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+geom_hline(yintercept=nCount_out,linetype="dashed")+labs(x=paste(sample,"Pre-filter"), y="Log2(nCount_RNA)")+
	labs(fill = paste(sample,"Pre-filter"),title="Pre-filter nCount_RNA")
dev.off()
#Post-filter nCount
pdf(paste0(outImageDir,"/nCount_postFilter_",sample,".pdf"))
gplot=ggplot(so_filt@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(nCount_RNA))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[2],labels=c("nCount"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+labs(x= paste(sample,"Post-filter"), y="Log2(nCount_RNA)")+
	labs(fill =  paste(sample,"Post-filter"),title="Post-filter nCount_RNA")
dev.off()
#Pre-filter mitochondrial
pdf(paste0(outImageDir,"/mitoPct_preFilter_",sample,".pdf"))
gplot=ggplot(so_BC@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(percent.mt))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[3],labels=c("Percent Mito"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+geom_hline(yintercept=mt_out,linetype="dashed")+labs(x=paste(sample,"Pre-filter"), y="Log2(Percent Mito)")+
	labs(fill = paste(sample,"Pre-filter"),title="Pre-filter Percent Mito")
dev.off()
#post-filter mitochondrial
pdf(paste0(outImageDir,"/mitoPct_postFilter_",sample,".pdf"))
gplot=ggplot(so_filt@meta.data,aes(x=orig.ident,y=log2(as.numeric(as.character(percent.mt))),fill=orig.ident))
gplot+geom_violin(trim=F)+scale_fill_manual(values=hue_pal()(3)[3],labels=c("Percent Mito"))+
	geom_boxplot(width=0.1,fill="white")+
	scale_x_discrete(labels=NULL)+
	theme_classic()+labs(x= paste(sample,"Post-filter"), y="Log2(Percent Mito)")+
	labs(fill =  paste(sample,"Post-filter"),title="Post-filter Percent Mito")
dev.off()

#Cluster analysis: UMAPs and silhouettes
for (res in resolutionString){
	pdf(paste0(outImageDir,"/clusterResolution_",res,"_",sample,".pdf"))
	resMod=as.numeric(gsub("\\.0$","",res))
	clusterPlot=DimPlot(so_noDoublet,group.by=paste0("SCT_snn_res.",resMod),label=T,repel=T)
	print(clusterPlot + labs(title = paste0("Sample: ",sample," at resolution ",res)))
	dev.off()
	
	pdf(paste0(outImageDir,"/silhouetteResolution_",res,"_",sample,".pdf"))

	Idents(so_noDoublet)=paste0("SCT_snn_res.",resMod)
	coord=Embeddings(so_noDoublet,reduction='pca')[,1:so@misc$npcs]
	clusters=Idents(so_noDoublet)
	d = dist(coord,method="euclidean")
	sil=silhouette(as.numeric(as.character(clusters)),dist=d)
	palette=alpha(colour=hue_pal()(length(unique(Idents(so_noDoublet)))),alpha=0.7)
	print(plot(sil, col=palette[as.factor(clusters[order(clusters,decreasing=F)])],
	main=paste0("Silhouette plot of clustering resolution ", res), lty=2,
	sub=paste("Average silhouette width:",format(round(mean(sil[,3]), 4), nsmall = 4))))
  
	abline(v=mean(sil[,3]), col="red4", lty=2)
	dev.off()

}

#Primary cell type annotation
pdf(paste0(outImageDir,"/primaryAnnotation_",sample,".pdf"))
singleRPlot=DimPlot(so_noDoublet,group.by="annot",label=T,repel=T)
print(singleRPlot + labs(title = paste0(sample," annotations by ",annotDB)))
dev.off()

pdf(paste0(outImageDir,"/primaryAnnotationCounts_",sample,".pdf"))
library(gridExtra)
df=data.frame(sort(table(so_noDoublet$annot),decreasing=T))
nTypes=length(table(so_noDoublet$annot))
names(df)=c("CellAnnot","Count")
if (nTypes <= 15){
	grid.table(df,rows=NULL)
}else{
	if ((nTypes%%2)==1){
		addition = as.data.frame(matrix(ncol=2,nrow=1))
		names(addition)=names(df)
		addition[]=""
		df = rbind(df,addition)
	}
	cutoff=nrow(df)/2
	grid.table(cbind(df[1:cutoff,],df[(cutoff+1):nrow(df),]),rows=NULL)
}
dev.off()

#All cell type annotations ####NEED TO FIGURE BEST WAY TO COUNT CELL TYPES-grid.table with modulus to max out at 4 columns between labels and counts
if(species == "human"){
	annotList=c("HPCA","BP_encode","monaco","dice","Novershtern")
	for (annot in annotList){
		pdf(paste0(outImageDir,"/cellAnnotUMAP_",annot,"_",sample,".pdf"))
		singleRPlot=DimPlot(so_noDoublet,group.by=paste0(annot,"_main"),label=T,repel=T)
		print(singleRPlot + labs(title = paste0(sample," annotations by ",toupper(annot))))
		dev.off()
		pdf(paste0(outImageDir,"/cellAnnotVln_",annot,"_",sample,".pdf"))
		print(plotScoreDistribution(so@misc$SingleR[[paste0(annot,"_main")]],show="delta.med",dots.on.top=F,size = 0.01,show.nmads=3))
		dev.off()
		
	}
}
if (species=="mouse"){
	annotList=c("immgen","mouseRNAseq")
	for (annot in annotList){
		pdf(paste0(outImageDir,"/cellAnnotUMAP_",annot,"_",sample,".pdf"))
		singleRPlot=DimPlot(so_noDoublet,group.by=paste0(annot,"_main"),label=T,repel=T)
		print(singleRPlot + labs(title = paste0(sample," annotations by ",toupper(annot))))
		dev.off()
		pdf(paste0(outImageDir,"/cellAnnotVln_",annot,"_",sample,".pdf"))
		plotScoreDistribution(so@misc$SingleR[[paste0(annot,"_main")]],show="delta.med",dots.on.top=F,size = 0.01,show.nmads=3)
		dev.off()
	}
}

#Doublet locations
pdf(paste0(outImageDir,"/doublets_",sample,".pdf"))
doubletPlot=DimPlot(so,group.by="DF_hi.lo")
print(doubletPlot+labs(title = paste0(sample," Doublets")))
dev.off()

#Cell cycle annotations
pdf(paste0(outImageDir,"/cellCycle_",sample,".pdf"))
cellCyclePlot=DimPlot(so_noDoublet,group.by="Phase")
print(cellCyclePlot+labs(title = paste0(sample," Cell Cycle")))
dev.off()


#### FINAL OUTPUT FILES ####
saveRDS(so,outRDS.doublets)
saveRDS(so_noDoublet,outRDS)
#save.image(outRdata)
