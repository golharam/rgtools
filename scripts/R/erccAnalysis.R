# get the list of names
rnaSamples <- read.table("rnaSamples.txt")
colnames(rnaSamples) <- c("sample", 'fq1', 'fq2', 'readCount1', 'readCount2', 'readlen1', 'readlen2')
sampleNames <- levels(rnaSamples[,1])

# Read in first sample to get chromosome names
sampleData <- read.table(paste("analysis/",sampleNames[1],"/",sampleNames[1],".idxStats.txt", sep=""))
chrNames <- levels(sampleData$V1)
erccNames <- chrNames[95:186]
data <- sampleData[95:186,3]
erccData <- data.frame(sample=data)
rownames(erccData) <- erccNames
colnames(erccData) <- gsub("-", ".", sampleNames[1])

for (i in 2:length(sampleNames)) {
  sampleData <- read.table(paste("analysis/",sampleNames[i],"/",sampleNames[i],".idxStats.txt", sep=""))
  data <- sampleData[95:186,3]
  erccData[i] <- data
  colnames(erccData)[i] <- gsub("-", ".", sampleNames[i
}

# Write out ERCC matrix
write.table(erccData, file="analysis/erccData.txt", quote=FALSE, sep="\t")

# read in rsem gene counts
rsem.genes.counts <- read.table("analysis/genes.counts.matrix")

# Combine the data frames by sample name
exprData <- rbind(rsem.genes.counts, erccData)

# Save matrix
write.table(exprData, file="analysis/exprERCC.txt", quote=FALSE, sep="\t")