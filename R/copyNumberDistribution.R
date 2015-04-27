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
      copyNumberDistributions.R Calculate genome-wide copy number distribution, copy number distribuion in control regions, non-control regions, autosomic non control regions.
      Gain/Loss CutOffs: Using just control regions, normalize copy number -using an StDev excluding windows with the 1% most extreme copy number- and define cut offs as 
      +-3 StDev of the normalized distribution
 
      Arguments:
      -sample=sampleName      - Sample Name
      -cw=copyWindowsBed      - Path to bed copy windows file
      -cn=copyNumberBed       - Path to the copy number bed file
      -rData=RDataFile        - Output path to RData results file
      -cutOffs=cutOffFile     - Output path to cutt offs results file [sample, mean, stDev_99, wins.excl, p.wins.excl, gain_cutoff, loss_cutoff]
      -plotFile=plot.png      - Output path R plot file
      -distribution=distribution.txt - Output path to control regions distribution file
      -help                      - print this text
 
       Example:
      ./copyNumberDistribution.R -cw=sample.calls.cw_norm.bed  -cn=sample.calls.copynumber.bed -rData=/CNdistribution/copy_number_distribution.RData 
                                 -cutOffs=/CNdistribution/sample.specific.cutoffs.txt -plotFile=/CNdistribution/controlRegions_density_plot.png 
                                 -distribution=/CNdistribution/controlRegions_distrib.txt\n\n")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^-", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
##############################################################################
## Data preparation
##############################################################################

## Samples
samples <-  argsL$sample

## Colors
c.cols <- colors()[grep("orange", colors())][1:length(c)]
my.cols <- c.cols

## Define paths
path.copyWindows <- argsL$cw
path.copyNumber <- argsL$cn
path.rdata <- argsL$rData
path.cutOffs <- argsL$cutOffs
path.plotFile <- argsL$plotFile
path.distribution <- argsL$distribution

#dir.create(file.path(path.out), recursive = TRUE)

## Import read depth and copy number per window
cn <- list()
cn.dist <- list()
cn.control.dist <- list()
cn.non.control.dist <- list()
cn.non.control.autosomes.dist <- list()
for (i in 1:length(samples)) {
	print(samples[i])
	## Import	
	path.import1 <- paste(path.copyWindows, sep = "")
	path.import2 <- paste(path.copyNumber, sep = "")
	x <- read.delim(path.import1)
	y <- read.delim(path.import2)
	names(x)[c(1, 4)] <- c("CHROM", "GC")
	names(y)[c(1, 4)] <- c("CHROM", "GC")
	# Merge
	z <- merge(x, y, merge = name(x)[1:4], sort = F)
	## Calculate copy number distribution in non-control regions but only on autosomes
	cn.non.control.autosomes.dist[[i]] <- data.frame(table(cut(z[z$CHROM != "chrUn" & z$CHROM != "chrX", ]$COPYNUMBER,
		breaks = c(seq(0, 100, 0.1), 1e3, 1e4, 1e5), right = T, include.lowest = T)))
	row.names(cn.non.control.autosomes.dist[[i]]) <- c(seq(0, 100, 0.1), 1e3, 1e4)

	cn[[i]] <- z
}

# Save data
save.image(path.rdata)

# estimate the maximum copy number for control/non-control in each dataset
k1 <- c()
k2 <- c()
for (i in 1:length(cn)) {
	k1 <- c(k1, max(cn[[i]][cn[[i]]$IS_CONTROL ==  "Y" , ]$COPYNUMBER))
	k2 <- c(k2, max(cn[[i]][cn[[i]]$IS_CONTROL ==  "N" , ]$COPYNUMBER))
}
max(k1)
max(k2)

##############################################################################
# Gain/loss cut-offs
##############################################################################

# Using only control regions, normalize copy number (using an Stdev excluding windows with the 1% most extreme copy number) and define cut-offs as +-3 Stdev of the normalized distribution
cutoffs <- data.frame(matrix(ncol = 7, nrow = 0))
names(cutoffs) <- c("sample", "mean", "stDev_99", "wins.excl", "p.wins.excl", "gain_cutoff", "loss_cutoff")
for (i in 1:length(samples)) {
	print(samples[i])
	# Subset control regions
	x <- cn[[i]]
	control <- x[x$IS_CONTROL == "Y", ]
	# Normalization with parameters
	# mean
	mu <- mean(control$COPYNUMBER)
	# Standard deviation (excluding windows with the 1% most extreme copy number)
	cn.q99 <- quantile(control$COPYNUMBER, probs = c(0.99))[[1]]
	stDev <- sd(control[control$COPYNUMBER < cn.q99, ]$COPYNUMBER)
	wins.excl <- length(control[control$COPYNUMBER > cn.q99, ]$COPYNUMBER)
	p.wins.excl <- (wins.excl * 100) / length(control$COPYNUMBER)
	control$COPYNUMBER_NORM <- (control$COPYNUMBER - mu) / stDev
	#print(c(mu, stDev))
 	# Cut-offs
	gain.cutoff <- (3 * stDev) + mu
	loss.cutoff <- (-3 * stDev) + mu
	y <- data.frame(samples[i], mu, stDev, wins.excl, p.wins.excl, gain.cutoff, loss.cutoff)
	names(y) <- names(cutoffs)
	cutoffs <- rbind(cutoffs, y)
}

cutoffs$species <- c(1:length(samples))

data4plot <- cutoffs
data4plot <- data4plot[length(samples):1, ]

# Export
write.table(data4plot, path.cutOffs, sep = "\t", quote = F, row.names = F)


##############################################################################
## Belen Lorente Plots
##############################################################################

## Samples
samplesLorente <- argsL$sample

## Define paths
path.densityPlot <- paste(path.plotFile, sep = "")
path.controlRegions <- paste(path.distribution, sep = "")


distrCR<-data.frame(Sample=c(),Min=c(),Mean=c(),Median=c(),StDev=c(),Max=c())
df0<-data.frame(CopyNumber=c(),Sample=c())

for (i in 1:length(samplesLorente)) 
{
	print(samplesLorente[i])
        read.table(paste(c(path.copyWindows), collapse=""),header=F)->calls
        CR<-calls[calls$V6=="Y" & calls$V1 != "chrX",]
        distrCR<-rbind(distrCR,data.frame(Sample=samplesLorente[i],Min=min(CR$V5),Mean=mean(CR$V5),Median=median(CR$V5),StDev=sd(CR$V5),Max=max(CR$V5)))
        df0<-rbind(df0,data.frame(CopyNumber=CR$V5/median(CR$V5)*2,Sample=rep(samples[i],times=length(CR$V5))))
}

png(path.densityPlot, h = 700, w = 700, res = 100)
library(ggplot2)
m=ggplot(data=df0,aes(x=CopyNumber,color=Sample)) + geom_density() + theme_bw()+ coord_cartesian(xlim=c(0.7,3.3))
m

dev.off()

write.table(distrCR,path.controlRegions,quote=F,sep="\t",row.names=F)

