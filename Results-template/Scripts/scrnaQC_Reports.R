args <- commandArgs(trailingOnly = TRUE)

inputDir <- as.character(args[1])
outputDir <- as.character(args[2])
resolution <- as.character(args[3])
sample = as.character(args[4])
#sample=strsplit(dir,".Rdat")[[1]][1]
#outDir = as.character(args[2])
rmarkdown::render("Scripts/scrnaQC.Rmd", output_file=paste0("QC_Report_",sample,".html"),output_dir = outputDir,
        params = list(dir=inputDir, resolution = resolution, sample= sample)
)
