pi.data <- data.frame(genus = character(0), sp = character(0), threshold = numeric(0), pi = numeric(0))

# Format the dataset:
# For each file, record genus, species and threshold from the file name
# Then estimate genome-wide average pi
# as the sum of differences divided by the sum of comparisons.
# Then append to the dataset
# Remove "FsiM_PB_v5_scf3", which is the social chromosome of Formica.
# It doesn't have any effect outside Formica since they don't have it
for (file in list.files("./pi/", pattern = "_pi.txt")) {
  data <- read.table(paste0("./pi/", file), h=T)
  genus <- unlist(strsplit(file, "_"))[1]
  sp <- paste(genus, unlist(strsplit(file, "_"))[2], sep="_")
  threshold <- 1 - as.numeric(unlist(strsplit(file, "_"))[3])
  pi <- sum(data$count_diffs[which(data$chromosome != "FsiM_PB_v5_scf3")], na.rm = T) / sum(data$count_comparisons[which(data$chromosome != "FsiM_PB_v5_scf3")], na.rm = T)
  df <- data.frame(genus = genus, sp = sp, threshold = threshold, pi = pi)
  pi.data <- rbind(pi.data, df)
}

pi.data$col <- NA
pi.data$col[which(pi.data$genus == "Camp")] <- "#29235C"
pi.data$col[which(pi.data$genus == "Form")] <- "#A3195B"
pi.data$col[which(pi.data$genus == "Lasi")] <- "#00A19A"
pi.data$col[which(pi.data$genus == "Myrm")] <- "#3BBDDC"
pi.data$col[which(pi.data$genus == "Tapi")] <- "#F39200"
pi.data$col[which(pi.data$genus == "Temn")] <- "#95C11F"
pi.data$col[which(pi.data$genus == "Tetr")] <- "#2D2E83"


pi.data$col.transparent = paste0(substr(pi.data$col, 1, 7), "30")


# Define which species hybridized
pi.data$hybrid <- F
for (sp in levels(as.factor(pi.data$sp))) {
  lines <- which(pi.data$sp == sp)
  if (length(lines)>1) {
    pi.data$hybrid[lines] <- T
  }
}

##################################################################################################
### This is the overly complicated version that was not retained in the paper, but here you go ###
##################################################################################################

# Define relative_pi as pi with the given threshold
# divided by pi with 0.9999 or 0.9998 (in Formica) as threshold for the same species
# (i.e. with pure individuals only)
pi.data$relative_pi <- NA
for (i in 1:length(pi.data$pi)) {
  sp <- pi.data$sp[i]
  pi.data$relative_pi[i] <- pi.data$pi[i] / pi.data$pi[which(pi.data$sp == sp & pi.data$threshold < 0.0002)]
}

par(mfrow=c(1,2))

# Plot pi as a function of the threshold for each species
plot(pi.data$pi~pi.data$threshold, pch=NA, las=1, ylab="pi", xlab="Introgression threshold")
for (sp in levels(as.factor(pi.data$sp))) {
  lines <- which(pi.data$sp == sp)
  if (sum(pi.data$relative_pi[lines] == 1)!=4) {
    with(pi.data[which(pi.data$sp == sp),], points(pi ~ threshold, type="l", pch=20, col=col, lwd=1.5))
  }
}

# Plot relative pi as a function of the threshold for each species
plot(pi.data$relative_pi~pi.data$threshold, pch=NA, las=1, ylab="Relative pi", xlab="Introgression threshold")
for (sp in levels(as.factor(pi.data$sp))) {
  lines <- which(pi.data$sp == sp)
  if (sum(pi.data$relative_pi[lines] == 1)!=4) {
    with(pi.data[which(pi.data$sp == sp),], points(relative_pi ~ threshold, type="l", pch=20, col=col, lwd=1.5))
  }
}



######################
### Easier version ###
######################
pure.pi <- numeric(0)
max.pi <- numeric(0)
max.relative.pi <- numeric(0)
col.pi <- character(0)

for (sp in levels(as.factor(pi.data$sp))) {
  if (pi.data$hybrid[which(pi.data$sp == sp)][1] == T) {
    pi.sp <- pi.data[which(pi.data$sp == sp),]
    pure.pi <- c(pure.pi, pi.sp$pi[which(pi.sp$threshold == min(pi.sp$threshold))])
    max.pi <- c(max.pi, pi.sp$pi[which(pi.sp$threshold == max(pi.sp$threshold))])
    max.relative.pi <- c(max.relative.pi, pi.sp$relative_pi[which(pi.sp$threshold == max(pi.sp$threshold))])
    col.pi <- c(col.pi, pi.sp$col[which(pi.sp$threshold == min(pi.sp$threshold))])
  }
}

col_trans.pi <- paste0(substr(col.pi, 1, 7), "80")


jitter = jitter(rep(1, length=length(pure.pi)), factor = 5)
jitter2 = 1 + jitter

par(mfrow=c(1,1))
boxplot(pure.pi, max.pi, col=F, border=F, las=1)
points(pure.pi ~ jitter, pch=19, col = col.pi)
points(max.pi ~ jitter2, pch=19, col = col.pi)
segments(x0=jitter, x1 = jitter2, y0 = pure.pi, y1 = max.pi, col = col.pi)

plot(density(max.pi/pure.pi, bw = 0.005), col=F, main = "", yaxt="n", ylab="", xlab="")
polygon(density(max.pi/pure.pi, bw = 0.005), col="lightgrey", border=F)
mean(max.pi/pure.pi)
range(max.pi/pure.pi)
abline(v=mean(max.pi/pure.pi))
