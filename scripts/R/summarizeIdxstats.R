# get list of sample names
# TBD: Make NA_samples.txt a parameter
samples <- read.table("NA_samples.txt")
colnames(samples) <- c('sample', 'fq1', 'fq2')
sampleNames <- as.vector(samples$sample)

# Read in first sample to get chromosome names
allSampleData <- read.table(paste("analysis/",sampleNames[1],"/",sampleNames[1],".idxStats.txt", sep=""))
colnames(allSampleData) <- c('chr', 'length', paste(sampleNames[1],'mappedReads',sep='.'), paste(sampleNames[1],'unmappedReads',sep='.'))

for (i in 2:length(sampleNames)) {
  # Read in sample
  sampleData <- read.table(paste("analysis/",sampleNames[i],"/",sampleNames[i],".idxStats.txt", sep=""))
  mappedReadsColName = paste(sampleNames[i],'mappedReads',sep='.')
  unmappedReadsColName = paste(sampleNames[i],'unmappedReads',sep='.')
  colnames(sampleData) <- c('chr', 'length', mappedReadsColName, unmappedReadsColName)
  # Append sample to master table
  allSampleData[, mappedReadsColName] <- sampleData[,mappedReadsColName]
  allSampleData[, unmappedReadsColName] <- sampleData[,unmappedReadsColName]
}

# Write out matrix
write.table(allSampleData, file="analysis/idxStats.txt", quote=FALSE, sep="\t", row.names=FALSE)


