# Get a list of samples and count of targets
# This script will generate 2 plots:
# Mean Coverage - Which can be used to see how samples compare with each other
# Normalized coverage - Which can be used to see how target within samples compare (only within a sample)

samples <- read.table('../samples.txt', header=FALSE)
x <- read.table(paste(as.matrix(samples$V1)[1], ".dedup_n6dupsRemoved_sorted.per_target_metrics.txt", sep=""), header=TRUE)

nc <- as.data.frame(matrix(nrow=nrow(x), ncol=nrow(samples)))
colnames(nc) <- samples$V1
rownames(nc) <- x$name

mc <- as.data.frame(matrix(nrow=nrow(x), ncol=nrow(samples)))
colnames(mc) <- samples$V1
rownames(mc) <- x$name

for (s in samples$V1) {
	x <- read.table(paste(s, ".dedup_n6dupsRemoved_sorted.per_target_metrics.txt", sep=""), header=TRUE)
	mc[s] <- x$mean_coverage
	nc[s] <- x$normalized_coverage
}

write.table(mc, file="mean_cov.txt", quote=FALSE, sep="\t")
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

nc_mat_data <- data.matrix(nc)
mc_mat_data <- data.matrix(mc)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

pdf("norm_coverage_per_target.pdf")
heatmap.2(nc_mat_data, 		# data
	  col=my_palette,	# color palette to use
	  cexRow=0.1,		# row label font size
	  cexCol=0.5,		# col label font size
	  srtCol=45		# angle col labels
	 )
dev.off()

pdf("mean_coverage_per_target.pdf")
heatmap.2(mc_mat_data,          # data
          col=my_palette,       # color palette to use
          cexRow=0.1,           # row label font size
          cexCol=0.5,           # col label font size
          srtCol=45             # angle col labels
         )
dev.off()

