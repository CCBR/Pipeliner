args <- commandArgs(trailingOnly = TRUE)

dir <- as.character(args[1])
setwd(dir)
#sample=strsplit(dir,".Rdat")[[1]][1]
#outDir = as.character(args[2])
rmarkdown::render("../Scripts/scrnaQC.Rmd", output_file="samples_QC.html", params = list(
  dir=dir
))
