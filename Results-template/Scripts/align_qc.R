args <- commandArgs(trailingOnly = TRUE)


stats = matrix()
stats2 = matrix()

for (F in args) {
    f = read.table(paste0(F, ".flags"), fill = TRUE)
    stats = cbind(stats, f["V1"])
}

stats = stats[, -1]
stats = stats[-12, ]

for (F in args) {
    f = read.table(paste0(F, ".dedup.flags"), fill = TRUE)
    stats2 = cbind(stats2, f["V1"])
}

stats2 = stats2[, -1]
stats2 = stats2[-12, ]

names(stats) = args
row.names(stats) = c("Total Reads", "Duplicates", "Mapped", "Paired in Sequencing", 
    "Read 1", "Read 2", "Properly Paired", "Mapped with Mate Mapped", "Singletons", 
    "Mate Mapped to Different Chr", "Mate Mapped to Different chr (mapQ>5)")

names(stats2) = names(stats)
row.names(stats2) = row.names(stats)

write.table(stats, "align_qc.tab", sep = "\t")
write.table(stats2, "align_qc2.tab", sep = "\t")


png("align_qc.png")

barplot(as.matrix(stats[c(1, 3, 7), ]), xlab = "Read ID", ylab = "Total---------------------------------Mapped-------------------------Properly Paired", 
    main = "Total, Mapped, and Properly Paired Reads")
 
