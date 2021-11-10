## Read the data (output of filtering_tests.sh)
filt <- read.table("summary_filtering.txt", h=T)
col.names(filt) <- c("DP", "max_miss", "kept")

## Create an empty plot
plot(0, type="n", xlim=range(filt$max_miss), ylim=c(0, max(filt$kept)), las=1, xlab="Proportion of individuals whith non-missing data", ylab="Number of loci retained")

## Add a set of points for each value of depth
for (i in unique(filt$DP)) {
    points(filt$kept[which(filt$DP==i)] ~ filt$max_miss[which(filt$DP==i)], col=seasun(17)[i-min(filt$DP)+1], pch=19, type="b")
}

## Add a legend
legend("topright", pch=19, col = seasun(17), legend=unique(filt$DP), bty="n", ncol=2)
