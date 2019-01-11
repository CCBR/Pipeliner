
args = commandArgs(trailingOnly=TRUE)

library(sequenza)
library(parallel)

seqz_file = args[1]
out_dir = args[2]
n_cores = as.numeric(args[3])

if (! file.exists(seqz_file)) {
    stop(paste0("Can't find this SEQZ output file: ", seqz_file))
}
if (is.na(n_cores)) {
    n_cores = 1
}
if (is.na(out_dir)) {
    out_dir = dirname(seqz_file)
}


print(paste0("Using ",n_cores," cores..."))
#date()
#print("Reading seqz file...")
#seqz.data <- read.seqz(seqz_file)
#str(seqz.data, vec.len = 2)

date()
print("Extracting seqz data...")
seqzdata <- sequenza.extract(seqz_file, min.reads = 30, min.reads.normal= 10)

date()
print("Fitting model...")
CP.example <- sequenza.fit(seqzdata, mc.cores = n_cores))

print("Printing results...")
sequenza.results(sequenza.extract = seqzdata,cp.table = CP.example, sample.id = out_dir, out.dir=out_dir)

date()
print("Done")
