# Example usage: Rscript runDiffBind.R 'directory' 'outfilename.html' 'input.csv' 'WT_vs_KO' 'macs_narrow' 'projectID' 'Enter CCBR Project Description and Notes here.'

args <- commandArgs(trailingOnly = TRUE)

DIR <- args[1]
setwd(DIR)
outHtml <- args[2]

Sys.setenv(RSTUDIO_PANDOC="/usr/local/apps/rstudio/rstudio-1.1.447/bin/pandoc/")

rmarkdown::render("DiffBind_pipeliner.Rmd",output_file=outHtml, params = list(
	csvfile = args[3],
	contrasts = args[4],
	peakcaller = args[5],
	projectID = args[6],
	projectDesc = args[7]
))