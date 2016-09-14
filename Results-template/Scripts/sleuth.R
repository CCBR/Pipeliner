args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
PIPERLIB <- args[2]
SAMTAB <- args[3]
CONTRASTS <- args[4]
SPECIES <- args[5]

setwd(paste(DIR,"/salmonrun",sep=""))

library(sleuth)
#library("wasabi", lib.loc=PIPERLIB)
library("wasabi")
library(biomaRt)

s2c<-read.table(file.path(paste(DIR,"/",SAMTAB,sep="")), header=TRUE, stringsAsFactors = FALSE)
s2c<-dplyr::select(s2c, sample=sampleName, condition)
s2c$sample<-gsub("\\..*","",s2c$sample)

salmon_dir<-getwd()
salmon_dir
# Get all sample dirs
sample_id <- s2c$sample
sample_id
sample_dirs <- sapply(sample_id, function(id) file.path(salmon_dir,id))
sample_dirs
s2c <- dplyr::mutate(s2c, path = sample_dirs)

for(sdir in s2c$path){
     prepare_fish_for_sleuth(sdir)
}

if (SPECIES == "hg19") {
   mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
} else {
   mart <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
}

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
  "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

contras=unlist(strsplit(CONTRASTS, split=" "))
for(i in seq(1, length(contras), by = 2)){
	s2c_pair = s2c[which(s2c$condition==as.character(contras[i])|s2c$condition==as.character(contras[i+1])),]
	pairname = paste(contras[i],"_vs_",contras[i+1],sep="")

	so <- sleuth_prep(s2c_pair, ~ condition, target_mapping = t2g)
	 
	so <- sleuth_fit(so)

	beta <- colnames(design_matrix(so))[2]

	so <- sleuth_wt(so, which_beta = beta)

	#save(so, file=paste(pairname,"_sleuth_object.R",sep=""))

	results_table <- sleuth_results(so, beta)
	write.csv(results_table, file=paste(pairname,"_sleuthResults_Tx.txt",sep=""))

	genetable<-sleuth_gene_table(so,beta,test_type ='wt',which_group="ens_gene")
	abundance<-kallisto_table(so, use_filtered = FALSE, normalized = FALSE,include_covariates = FALSE)
	write.csv(abundance, file=paste(pairname,"_sleuth_abundance.txt",sep=""))

	pdf(file=paste(pairname,"_isoform_DEG.pdf",sep=""), title=pairname)
	print(plot_pca(so))
	print(plot_sample_heatmap(so))
	print(plot_group_density(so))
	print(plot_ma(so,beta))
	print(plot_qq(so,beta))
	#FOR SOME REASON THIS COMMAND DOES NOT WORK, IT SEEMS TO BE BECAUSE OF ILLEGAL CHARACTERS IN DIR NAMES?
	#plot_scatter(so,sample_x = so$sample_to_covariates$sample[1],sample_y = so$sample_to_covariates$sample[2])
	print(plot_volcano(so,beta))
	dev.off()	

}

write.csv(s2c, file="sleuth_completed.txt")

