args <- commandArgs(trailingOnly = TRUE)

sample = args[1]

T = read.table("out")

chrs = c("M", 1:22, "X", "Y")
seg = list()

for (i in chrs) {
    theC = paste0("chr", i)
    if (any(T[, 1] == theC)) {
        seg[[i]] = which(T[, 1] == theC)
    } else {
        seg[[i]] = 0
    }
}

window = 10000

for (j in 1:length(nchar(chrs))) {
    
    if (seg[[j]][1] == 0) {
        next
    }
    t = T[seg[[j]], ]
    C = matrix(ncol = 2)
    for (i in seq(1, dim(t)[1], window)) {
        k = sum(t[i:(i + window), 3])/window
        
        C = rbind(C, c(t[i, 2], k))
        
    }
    
    png(paste0(sample, ".chr", chrs[j], ".depth.png"))
    plot(C, type = "b", main = paste(sample, "chr", chrs[j]), xlab = paste("Base Position"), 
        ylab = "Sequencing Depth", col = "brown")
    dev.off()
    
    png(paste0(sample, ".chr", chrs[j], ".coverage.png"))
    hist(t[, 3], main = paste(sample, "chr", chrs[j]), xlab = "Sequencing Depth", 
        ylab = "Frequency", col = "#fdfadd")
    dev.off()
    
}

 
