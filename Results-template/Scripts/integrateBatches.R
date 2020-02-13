library(Biobase,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library(farver,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library(S4Vectors,lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library("SingleR",lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library(scRNAseq)
library(SingleCellExperiment)
library("destiny",lib.loc="/data/CCBR_Pipeliner/db/PipeDB/Rlibrary_3.6.0_scRNA")
library("URD",lib.loc="/data/CCBR_Pipeliner/db/PipeDB/scrna_lib")
library("Seurat")
library(dplyr)
library(Matrix) 
library(tools)
library(stringr)



args <- commandArgs(trailingOnly = TRUE)

matrix <- as.character(args[1])
#output = as.character(args[2])
outDirSeurat = as.character(args[2])
outDirMerge = as.character(args[3])
specie = as.character(args[4])
resolution = as.numeric(args[5])
clusterAlg =  as.numeric(args[6])
annotDB = as.character(args[7])
nAnchors = as.numeric(args[8])
groups = as.character(args[9])
contrasts = as.character(args[10])



file.names <- dir(path = matrix,pattern ="rds")

if (groups == "YES") {
    
groupFile = read.delim("groups.tab",header=F,stringsAsFactors = F)
groupFile=groupFile[groupFile$V2 %in% stringr::str_split_fixed(contrasts,pattern = "-",n = Inf)[1,],]
#groupFile = groupFile[groupFile$V2 == strsplit(contrasts,"-")[[1]][1] | groupFile$V2 == strsplit(contrasts,"-")[[1]][2] ,] 

splitFiles = gsub(".rds","",file.names)#str_split_fixed(file.names,pattern = "[.rd]",n = 2) 
file.names=file.names[match(groupFile$V1,splitFiles,nomatch = F)]
print(groupFile$V1)
print(splitFiles)
print(file.names)
}

readObj = list()
for (obj in file.names) {
  Name=strsplit(obj,".rds")[[1]][1]
  assign(paste0("S_",Name),readRDS(paste0(matrix,"/",obj)))
  readObj = append(readObj,paste0("S_",Name))  
 }

#for (obj in readObj) {
 # print(obj)  	
  #sample = CreateSeuratObject(GetAssayData(object = eval(parse(text =obj)),slot = "counts")[,colnames(x = eval(parse(text =obj)))])
  #sample@assays$RNA@data = GetAssayData(object = eval(parse(text =obj)),slot = "data")[,colnames(x = eval(parse(text =obj)))]
  #sample@meta.data =eval(parse(text =obj))@meta.data
  #assign(obj,sample)

#}


combinedObj.list=list()
i=1
for (p in readObj){
  combinedObj.list[[p]] <- eval(parse(text = readObj[[i]]))
  i <- i + 1
 }


reference.list <- combinedObj.list[unlist(readObj)]
print(reference.list)
for (i in 1:length(x = reference.list)) {
#    reference.list[[i]] <- NormalizeData(object = reference.list[[i]], verbose = FALSE)
    reference.list[[i]] <- FindVariableFeatures(object = reference.list[[i]], 
        selection.method = "vst", nfeatures = nAnchors, verbose = FALSE)
}

print(length(reference.list))
print(reference.list)
combinedObj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30,anchor.features = nAnchors)
combinedObj.integrated <- IntegrateData(anchorset = combinedObj.anchors, dims = 1:30)


DefaultAssay(object = combinedObj.integrated) <- "integrated"
combinedObj.integrated <- ScaleData(object = combinedObj.integrated, verbose = FALSE)



combinedObj.integratedRNA = combinedObj.integrated
DefaultAssay(object = combinedObj.integratedRNA) <- "SCT"

combinedObj.integratedRNA = FindVariableFeatures(combinedObj.integratedRNA,mean.cutoff = c(0.0125, 3),dispersion.cutoff = c(0.5, Inf),selection.method = "vst")
combinedObj.integratedRNA <- ScaleData(object = combinedObj.integratedRNA, verbose = FALSE)


mat1 <- t(as.matrix(FetchData(object = combinedObj.integratedRNA,slot = "counts",vars = rownames(combinedObj.integratedRNA))))
urdObj <- createURD(count.data = mat1, min.cells=3, min.counts=3)
varGenes_batch = VariableFeatures(combinedObj.integrated)[VariableFeatures(combinedObj.integrated) %in% rownames(urdObj@logupx.data)]
varGenes_merge = VariableFeatures(combinedObj.integratedRNA)[VariableFeatures(combinedObj.integratedRNA) %in% rownames(urdObj@logupx.data)]

pcs_batch <- calcPCA(urdObj, pcs.store = 50,genes.use=varGenes_batch,mp.factor = 1.2)
npcs_batch =  sum(pcs_batch@pca.sig)
print(npcs_batch)

pcs_merge <- calcPCA(urdObj, pcs.store = 50,genes.use=varGenes_merge,mp.factor = 1.2)
npcs_merge =  sum(pcs_merge@pca.sig)
print(npcs_merge)


runInt = function(obj,res,npcs){
obj <- RunPCA(object = obj, npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj,dims = 1:npcs)
obj <- FindClusters(obj, reduction.type = "pca", dims.use = 1:npcs, save.SNN = T,resolution = res,algorithm = clusterAlg)

obj <- RunUMAP(object = obj, reduction = "pca", 
                                  dims = 1:npcs,n.components = 3)

obj$groups = groupFile$V2[match(obj$Sample,  groupFile$V1,nomatch = F)]



runSingleR = function(obj,refFile,fineORmain){
avg = AverageExpression(obj,assays = "SCT")
avg = as.data.frame(avg)
ref = refFile
s = SingleR(test = as.matrix(avg),ref = ref,labels = ref[[fineORmain]])

clustAnnot = s$labels
names(clustAnnot) = colnames(avg)
names(clustAnnot) = gsub("SCT.","",names(clustAnnot))

obj$clustAnnot = clustAnnot[match(obj$seurat_clusters,names(clustAnnot))]
return(obj$clustAnnot)
}

if(annotDB == "HPCA"){
obj$clustAnnot <- runSingleR(obj,HumanPrimaryCellAtlasData(),"label.main")
obj$clustAnnotDetail <-  runSingleR(obj,HumanPrimaryCellAtlasData(),"label.fine")
}

if(annotDB == "BP_encode"){
obj$clustAnnot <-  runSingleR(obj,BlueprintEncodeData(),"label.main")
obj$clustAnnotDetail <-  runSingleR(obj,BlueprintEncodeData(),"label.fine")
}
if(annotDB == "monaco"){
obj$clustAnnot <-  runSingleR(obj,MonacoImmuneData(),"label.main")
obj$clustAnnotDetail <-     runSingleR(obj,MonacoImmuneData(),"label.fine")
}
if(annotDB == "immu_cell_exp"){
obj$clustAnnot <-  runSingleR(obj,DatabaseImmuneCellExpressionData(),"label.main")
obj$clustAnnotDetail <- runSingleR(obj,DatabaseImmuneCellExpressionData(),"label.fine")
}

if(annotDB == "immgen"){
obj$clustAnnot <-  runSingleR(obj,ImmGenData(),"label.main")
obj$clustAnnotDetail <- runSingleR(obj,ImmGenData(),"label.fine")
}
if(annotDB == "mouseRNAseq"){
obj$clustAnnot <-  runSingleR(obj,MouseRNAseqData(),"label.main")
obj$clustAnnotDetail <- runSingleR(obj,MouseRNAseqData(),"label.fine")
}
return(obj)
}

combinedObj.integrated = runInt(combinedObj.integrated,0.3,npcs_batch)
saveRDS(combinedObj.integrated, paste0(outDirSeurat,"_0.3.rds"))
combinedObj.integrated = runInt(combinedObj.integrated,0.6,npcs_batch)
saveRDS(combinedObj.integrated, paste0(outDirSeurat,"_0.6.rds"))
combinedObj.integrated = runInt(combinedObj.integrated,0.8,npcs_batch)
saveRDS(combinedObj.integrated, paste0(outDirSeurat,"_0.8.rds"))
combinedObj.integrated = runInt(combinedObj.integrated,1.0,npcs_batch)
saveRDS(combinedObj.integrated, paste0(outDirSeurat,"_1.0.rds"))
combinedObj.integrated = runInt(combinedObj.integrated,1.2,npcs_batch)
saveRDS(combinedObj.integrated, paste0(outDirSeurat,"_1.2.rds"))



combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,0.3,npcs_merge)
saveRDS(combinedObj.integratedRNA, paste0(outDirMerge,"_0.3.rds"))

combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,0.6,npcs_merge)
saveRDS(combinedObj.integratedRNA, paste0(outDirMerge,"_0.6.rds"))

combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,0.8,npcs_merge)
saveRDS(combinedObj.integratedRNA, paste0(outDirMerge,"_0.8.rds"))

combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,1.0,npcs_merge)
saveRDS(combinedObj.integratedRNA, paste0(outDirMerge,"_1.0.rds"))

combinedObj.integratedRNA = runInt(combinedObj.integratedRNA,1.2,npcs_merge)
saveRDS(combinedObj.integratedRNA, paste0(outDirMerge,"_1.2.rds"))


