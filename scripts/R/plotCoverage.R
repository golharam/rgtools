# Get a list of samples and count of targets
samples <- read.table('samples.txt', header=FALSE)
x <- read.table(paste(as.matrix(samples$V1)[1], "/", as.matrix(samples$V1)[1], ".per_target_metrics.txt", sep=""), header=TRUE)


nc <- as.data.frame(matrix(nrow=nrow(x), ncol=nrow(samples)))
colnames(nc) <- samples$V1
rownames(nc) <- x$name

for (s in samples$V1) {
	x <- read.table(paste(s, "/", s, ".per_target_metrics.txt", sep=""), header=TRUE)
	nc[s] <- x$normalized_coverage
}

write.table(nc, file="norm_cov.txt", quote=FALSE, sep="\t")

# Plot
if (!require("gplots")) {
	install.packages("gplots", dependencies = TRUE)
	library(gplots)
}
if (!require("RColorBrewer")) {
	install.packages("RColorBrewer", dependencies = TRUE)
	library(RColorBrewer)
}

mat_data <- data.matrix(nc)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

pdf("coverage_per_target.pdf")
#heatmap.2(mat_data, 		# data
#	  sepcolor="white",
#	  sepwidth=c(0.05,0.05),
#	  col=my_palette)	# color palette to use
heatmap.2(mat_data, 		# data
	  col=my_palette,	# color palette to use
	  cexRow=0.1,		# row label font size
	  cexCol=0.5,		# col label font size
	  srtCol=45		# angle col labels
	 )
dev.off()

