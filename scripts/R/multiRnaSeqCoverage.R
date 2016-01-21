# Script to generate a normalized coverage vs. position along transcript plot and overlay for multiple samples
# Based off of rnaseqCoverage.R from Picard Tool's CollectRNASeqMetrics

# @author Ryan Golhar

# Parse the arguments:
# Arg0 - ArgN-2 Metrics files
# ArgN-2 output metric txt file
# ArgN output png file
args <- commandArgs(trailing = TRUE)
metricsFiles <- vector(mode="character", length=length(args)-2)
for (i in 1:length(args)-2) {
  metricsFiles[i]  <- args[i]
}
outputTxt    <- args[length(args)-1]
outputFile   <- args[length(args)] # Get the last argument as the output file

# Figure out where the metrics and the histogram are in the first file and parse them out
print(paste("Reading", metricsFiles[1], sep=" "))
startFinder <- scan(metricsFiles[1], what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder)) {
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1
                } else {
                        secondBlankLine=i+1
                        break
                }
        }
}

data <- read.table(metricsFiles[1], header=T, sep="\t", skip=secondBlankLine, check.names=FALSE)

# The histogram has a normalized_position and normalized_coverage column for each metric "level"
# This code parses out the distinct levels so we can output one graph per level
headers <- sapply(sub(".normalized_coverage","",names(data),fixed=TRUE), "[[" ,1)
colnames(data) <- c("normalized_position", head(unlist(strsplit(tail(unlist(strsplit(metricsFiles[1], "/")), n=1), "[.]")), n=1))

## Duplicated header names cause this to barf. KT & Yossi report that this is going to be extremely difficult to
## resolve and it's unlikely that anyone cares anyways. Trap this situation and avoid the PDF so it won't cause
## the workflow to fail
#if (any(duplicated(headers))) {
#  print(paste("Not creating insert size PDF as there are duplicated header names:", headers[which(duplicated(headers))]))
  # TBD: Exit here
#}

# Read in the metrics for the rest of the samples
sampleMetrics <- data

for (i in 2:length(metricsFiles)) {
    print(paste("Reading", metricsFiles[i], sep=" "))
    data <- read.table(metricsFiles[i], header=T, sep="\t", skip=secondBlankLine, check.names=FALSE)
    colnames(data) <- c("normalized_position", head(unlist(strsplit(tail(unlist(strsplit(metricsFiles[i], "/")), n=1), "[.]")), n=1))
    sampleMetrics[i+1] <- data[2]
}

# Write out the metrics
write.table(sampleMetrics, file=outputTxt, quote = FALSE, sep="\t", row.names = FALSE)

# Determine the overall y-limit
ylim <- range(0, sampleMetrics[,2:length(sampleMetrics)])


levels <- c()
for (i in 2:length(headers)) {
    if (!(headers[i] %in% levels)) {
        levels[length(levels)+1] <- headers[i]
    }
}

# Some constants that are used below
COLORS = c("royalblue", "#FFAAAA", "palegreen3");

#pdf(outputFile)
png(outputFile)

# For each level, plot of the normalized coverage by GC
#for (i in 1:length(levels)) {
i <- 1
    # Reconstitutes the histogram column header for this level
    nc <- paste(levels[i], "normalized_coverage", sep=".")

    # plot the first sample, then we will iteratively add the others
    plot(x=sampleMetrics$normalized_position, y=as.matrix(sampleMetrics[2]),
         type="o",
         xlab="Normalized Distance Along Transcript",
         ylab="Normalized Coverage",
         xlim=c(0, 100),
         ylim=ylim,
         col="royalblue",
         main=paste("RNA-Seq Coverage vs. Transcript Position\n", levels[i], sep=""))

    # Add a horizontal line at coverage=1
    abline(h=1, col="lightgrey")

    for (j in 3:length(sampleMetrics)) {
      lines(x=sampleMetrics$normalized_position, y=as.matrix(sampleMetrics[j]), type='o', col="royalblue")
    }
#}

# redraw the first sample in a different color (because its a good sample)
lines(x=sampleMetrics$normalized_position, y=as.matrix(sampleMetrics[2]), type='o', col="red")

dev.off()
