#!/usr/bin/env Rscript

##1. COLLECT ARGUMENTS
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      makeHistogramCounts.R Plots the accumulative distribution from the counting of kmer mappings
 
      Arguments:
      --countFile=kmerCount.txt        - kmerCount, File where are stored the kmer mapping countings FORMAT [chrom start end sequence counts]
      --pngPlot=distributionPlot.png   - Png output file, for the accumalative distribution plot
      --histoPlot=histogramPlot.png    - Png output file, for the histogram plot
      --kmerOverThreshold=kmerCountOver- output Kmer bed file of regions over a given number of mapping hits
      --threshold=20                   - Number of mappings to consider a overrepresented kmer        
      --help              - print this text
      
      Example:
      ./makeHistogramCounts.R --countFile=kmerCounts.bed --pngPlot=distributionPlot.png --histoPlot=histogramPlot.png --kmerOverThreshold=kmerOver20.bed --threshold=20 \n\n")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
countFile <- argsL$countFile
pngPlot <- argsL$pngPlot
histoPlot <- argsL$histoPlot
outputFile <- argsL$kmerOverThreshold
threshold <- argsL$threshold

##2 IMPORT AND EDIT 
## 2.1 change to the new directory
counts <- read.delim(countFile, h = F)
names(counts) <- c("chrom", "start", "end", "seq", "counts")

## How many kmers?
length(counts[,1])

## Cumulative distribution
counts.table <- data.frame(table(counts$counts))
names(counts.table) <- c("counts", "num.kmers")
counts.table <- transform(counts.table, freq = counts.table$num.kmers / sum(counts.table$num.kmers))
counts.table$cum.freq <- rep(0, length(counts.table[, 1]))
for (i in 1:length(counts.table[, 1])) {
        counts.table[i, ]$cum.freq <- sum(counts.table[1:i, ]$freq)
}


## Plot
png(pngPlot, h = 700, w = 700, res = 100)
plot(counts.table$counts, counts.table$cum.freq, xlim = c(0, 100),
        ylab = "Cumulative percentage",
        xlab = "# placements for a given K-mer",
        main = "(K-mer size = 36, step = 5)")
abline(v = 19, col = "blue", lwd = 2)
legend("topright", "20 placements", lty = 1, lwd = 2, col = "blue", bty = "n")
dev.off()

## Select those kmers with at least 20 placements
write.table(counts[counts$counts >= threshold, ], outputFile, sep = "\t", quote = F, row.names = F, col.names = F)

hist(counts$counts, breaks = 1e5, xlim = c(0, 100))
png(histoPlot, h = 960, w = 960, res = 100)
#x11(w = 12, h = 10)
par(mfrow = c(2, 1))
plot(counts.table$counts, counts.table$num.kmers,
        ylab = "# K-mers",
        main = "All K-mers",
        xlab = "# placement for a given K-mer")
abline(h = mean(counts.table$num.kmers), col = "orange", lwd = 2)
legend("topright", "Threshold: mean", lty = 1, lwd = 2, col = "orange", bty = "n")
plot(counts.table$counts, counts.table$num.kmers,
        ylab = "# K-mers",
        main = "K-mers with <50 counts",
        xlab = "# placements for a given K-mer",
        xlim = c(0, 50))
abline(h = mean(counts.table$num.kmers), col = "orange", lwd = 2)
legend("topright", "Threshold: mean", lty = 1, lwd = 2, col = "orange", bty = "n")
dev.off()


