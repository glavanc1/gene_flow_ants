library(vcfR)
library(adegenet)
library(phytools)
library(MASS)

prop.het<-function(vector) {
  return(sum(vector=="0/1" | vector=="1/0", na.rm=T)/sum(is.na(vector)==F))
}

count.het<-function(vector) {
  return(sum(vector=="0/1" | vector=="1/0", na.rm=T))
}

count.na<-function(x) sum(is.na(x))

'%!in%' <- Negate('%in%')

vcf <- read.vcfR("Lasius_1snp_DP8_meanDP200_mac2_miss75_ADR_filtered_discard0.05_correct0.25.vcf")
gt <- extract.gt(vcf, "GT")
dp <- extract.gt(vcf, "DP", as.numeric=T)
gl <- vcfR2genlight(vcf)

individuals <- dimnames(gt)[[2]]

meta.all <- read.csv("Metadata_combined_public_scientific.csv", h=T, sep=";")

meta <- meta.all[match(dimnames(gt)[[2]], meta.all$CATALOGUENUMBER),]#[which(is.na(meta$CATALOGUENUMBER)==F),]


meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_alie")] <- "#27AAE1"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_brun")] <- "#006838"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_emar")] <- "#8DC63F"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_flav")] <- "#F9ED32"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_fuli")] <- "#1B75BC"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_mixt")] <- "#2B3990"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_myop")] <- "#92278F"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_nige")] <- "#EC008C"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_para")] <- "#F15A29"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_plat")] <- "#A97C50"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_psam")] <- "#00A79D"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Lasi_umbr")] <- "#C2B59B"


legend_Lasius <- c("L. niger", "L. platythorax", "L. emarginatus", "L. alienus", "L. paralienus", "L. fuliginosus", "L. brunneus", "L. psammophilus", "L. flavus", "L. myops", "L. mixtus", "L. umbratus", "not identified")
col_Lasius <- c("#EC008C", "#A97C50", "#8DC63F", "#27AAE1", "#F15A29", "#1B75BC", "#006838", "#00A79D", "#F9ED32", "#92278F", "#2B3990", "#C2B59B", "grey")

tree.not.pruned <- read.tree(".Lasius.tree")
tree <- drop.tip(phy = tree.not.pruned, tip = tree.not.pruned$tip.label[tree.not.pruned$tip.label %in% meta$CATALOGUENUMBER == F], trim.internal = T)

######################
### Heterozygosity ###
######################

meta$percent_missing <- apply(gt, 2, count.na) / dim(gt)[1]

meta$het <- apply(gt, 2, prop.het)

par(mfrow=c(1,1))
plot(meta$het ~ meta$percent_missing, col=meta$col.morpho, pch=19, las=1, xlab="Proportion of missing data", ylab="Relative heterozygosity")
legend("topleft", legend = legend_Lasius, pch=15, col=col_Lasius, bty="n", ncol=2)

###################
### mtDNA phylo ###
###################

tree.meta <- data.frame(ID = tree.not.pruned$tip.label, sp = NA, col = "grey", pch = NA, legend = NA)
for (i in tree.meta$ID) {
  n <- which(tree.meta$ID == i)
  if (i %in% meta$CATALOGUENUMBER) {
    j <- which(meta$CATALOGUENUMBER == i)
    tree.meta$sp[n] <- meta$SPECIESSUMMARY[j]
    tree.meta$col[n] <- meta$col.morpho[j]
    tree.meta$pch[n] <- 16
    tree.meta$legend[n] <- meta$CATALOGUENUMBER[j]
  } 
}

plot.phylo(tree.not.pruned, cex=0.2, show.node.label = F, show.tip.label = F, align.tip.label = T, no.margin = T )
tiplabels(pch=tree.meta$pch[match(tree.not.pruned$tip.label, tree.meta$ID)], col=tree.meta$col[match(tree.not.pruned$tip.label, tree.meta$ID)], cex=0.05, frame="none", offset=0.002)
tiplabels(text = tree.meta$legend[match(tree.not.pruned$tip.label, tree.meta$ID)], cex=0.05, frame="none", offset=0.005, adj = 0)
add.scale.bar()
legend("bottomleft", legend_Lasius, col=col_Lasius,  pch=15, bty="n")

###########
### MDS ###
###########

par(mfrow=c(1,1))
dist_all <- dist(gl)
MDS_all <- isoMDS(dist_all, k = 2, maxit = 100, p = 1)
plot(MDS_all$points[,2] ~ MDS_all$points[,1], col=meta$col.morpho, pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="Lasius - MDS", cex=10 * sqrt(meta$het))
legend("bottomright", legend=legend_Lasius, pch=19, col=col_Lasius, bty="n", ncol = 1)

nig_grp <- dimnames(gt)[[2]][which(MDS_all$points[,1] < -10)]
nig.vcf <- vcf[, c(0, which(MDS_all$points[,1] < -10)) + 1]
nig.gl <- vcfR2genlight(nig.vcf)
dist_nig <- dist(nig.gl)
MDS_nig <- isoMDS(dist_nig, k = 2, maxit = 100, p = 1)
plot(MDS_nig$points[,2] ~ MDS_nig$points[,1], col=meta$col.morpho[which(MDS_all$points[,1] < -10)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="L. niger group - MDS", cex=5 * sqrt(meta$het))
text(x=MDS_nig$points[,1], y=MDS_nig$points[,2], labels = nig_grp, cex=0.3)


#################
### admixture ###
#################

### Choose best K

CV <- read.table("admixture/CV_summary_across_rep.txt", h=F)
CV <- CV[order(CV$V1),]
dimnames(CV)[[1]] <- 1:12
K <- dim(CV)[1] - 1

plot(NULL, xlim=range(CV$V1), ylim=range(CV[,2:dim(CV)[2]]), las=1, ylab="Cross-Validation error", xlab = "K")
for (i in 2:dim(CV)[2]) {
  points(CV[,i] ~ CV$V1, type="b", pch=20, col=rgb(0, 0, 0, 0.2))
}

CV.mean <- apply(CV[,-1], 1, mean)
CV.sd <- apply(CV[,-1], 1, sd)
plot(CV.mean, pch=20, cex=0.75, col=NULL, las=1)
points(CV.mean, type="l", col="lightgrey")
for (i in 1:K) segments(x0=i, x1=i, y0=CV.mean[i] - CV.sd[i], y1=CV.mean[i] + CV.sd[i])
# we chose K = 9 based on these results
# the replicate with the lowest CV error was replicate 6

k9 <- read.table("admixture/rep6/cleaned.9.Q")
row.names(k9) <- dimnames(gt)[[2]]

par(mfrow=c(3,1), mai=c(0.1, 0.6, 0.1, 0.1))
barplot(meta$het[match(tree$tip.label, individuals)], las=1, ylab="Relative heterozygosity", border=NA)
barplot(t(as.matrix(k9))[,match(tree$tip.label, individuals)], las=2, border=NA, space=0, cex.names=0.6, col=c("#BAB023", "#F9ED32", "#A97C50", "#F15A29", "#EC008C", "#2B3990", "#8DC63F", "#006838", "#1B75BC"))
plot.phylo(tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
legend("bottomleft", legend=legend_Lasius, pch=15, col=col_Lasius, bty="n")
tiplabels(pch = 15, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])

order.structure <- function(Q.matrix){
  k <- dim(Q.matrix)[2]
  whichismax <- function (x) which(x == max(x, na.rm = T))
  assignment <- apply(Q.matrix, 1, whichismax)
  order <- numeric(0)
  for (i in 1:k){
    this_k <- which(assignment == i)
    order <- c(order, this_k[order(Q.matrix[this_k,i])])
  }
  return(Q.matrix[order,])
}

k9_ordered <- as.matrix(order.structure(k9))
par(mfrow=c(3,1))

barplot(meta$percent_missing[match(row.names(k9_ordered), individuals)], las=1, ylab="Missing data", border=NA)
barplot(meta$het[match(row.names(k9_ordered), individuals)], las=1, ylab="Relative heterozygosity", border=NA)
barplot(t(k9_ordered), col=c("#1B75BC", "#F15A29", "#A97C50", "#C2B59B", "#8DC63F", "#006838", "#F9ED32", "#92278F", "#EC008C"), las = 2, cex.names = 0.1, space = 0, border = NA)

#########################
### haploid admixture ###
#########################

choose.best.K <- function(CV_summary) {
  # Order CV summary by 1st column (K)
  # so they are ordered by value, not alphabetically
  CV <- CV_summary[order(CV_summary[,1]),]
  
  # plot all values from all replicates
  par(mfrow=c(1,2))
  plot(NULL, xlim=range(CV$V1), ylim=range(CV[,2:dim(CV)[2]]), las=1, ylab="Cross-Validation error", xlab = "K")
  for (i in 2:dim(CV)[2]) {
    points(CV[,i] ~ CV$V1, type="b", pch=20, col=rgb(0, 0, 0, 0.2))
  }
  
  # plot average +/- SD across replicates
  CV.mean <- apply(CV[,-1], 1, mean)
  CV.sd <- apply(CV[,-1], 1, sd)
  plot(CV.mean, pch=20, cex=0.75, col=NULL, las=1, ylab="Cross-Validation error (mean +/- SD)")
  points(CV.mean, type="l", col="lightgrey")
  for (i in 1:K) segments(x0=i, x1=i, y0=CV.mean[i] - CV.sd[i], y1=CV.mean[i] + CV.sd[i])
  
}


hap.CV <- read.table("haploid_admixture//CV_summary_across_rep.txt")
choose.best.K(hap.CV)

###############################
### Understanding 3-hybrids ###
###############################

# Is it caused by missing data?

plot(meta$percent_missing ~ jitter(apply(k9, 1, max), factor = 1000), pch=20, col="#00000030", las=1, xlab="Assignment to main cluster (+ noise)", ylab="Missing data")
boxplot(meta$percent_missing[which(apply(k9, 1, max) >= 0.9999)], meta$percent_missing[which(apply(k9, 1, max) < 0.9999)])

get2nd <- function(vector) return(sort(vector, decreasing = T)[2])
plot(meta$percent_missing ~ apply(k9, 1, get2nd))

get3rd <- function(vector) return(sort(vector, decreasing = T)[3])
plot(meta$percent_missing ~ apply(k9, 1, get3rd))
not_hybrids <- which(apply(k9, 1, max) >= 0.9999)
single_hybrids <- which(apply(k9, 1, get3rd) < 0.0001 & apply(k9, 1, get2nd) > 0.0001)
triple_hybrids <- which(apply(k9, 1, get3rd) > 0.0001)


boxplot(meta$percent_missing[not_hybrids], meta$percent_missing[single_hybrids], meta$percent_missing[triple_hybrids], las=1, outline = F, col = NULL, ylim=c(0, 0.55), ylab="Missing data", names=c("Pure species", "Hybrids: 2 spp.", "Hybrid: 3 spp."))
points(meta$percent_missing[not_hybrids] ~ jitter(rep(1, length = length(not_hybrids)), amount = 0.3), pch=20, col="#00000030")
points(meta$percent_missing[single_hybrids] ~ jitter(rep(2, length = length(single_hybrids)), amount = 0.3), pch=20, col="#00000030")
points(meta$percent_missing[triple_hybrids] ~ jitter(rep(3, length = length(triple_hybrids)), amount = 0.3), pch=20, col="#00000030")
text(x = c(1,2,3), y = 0.54, label = c("a", "ab", "b"))


meta$hybrid <- NA
meta$hybrid[which(apply(k9, 1, max) >= 0.9999)] <- "no"
meta$hybrid[which(apply(k9, 1, get3rd) < 0.0001 & apply(k9, 1, get2nd) > 0.0001)] <- "2sp"
meta$hybrid[which(apply(k9, 1, get3rd) > 0.0001)] <- "3sp"
aov_miss <- aov(meta$percent_missing ~ meta$hybrid)
summary(aov_miss)
plot(aov_miss)
kruskal.test(meta$percent_missing ~ meta$hybrid)
TukeyHSD(aov_miss)


# Distribution?
meta$X[which(meta$X==0)] <- NA
meta$Y[which(meta$Y==0)] <- NA
plot(meta$Y ~ meta$X, pch=20, col="#00000030", las=1, main="Imagine there's a map background")
points(meta$Y[which(meta$hybrid=="2sp")] ~ meta$X[which(meta$hybrid=="2sp")], col="blue")
points(meta$Y[which(meta$hybrid=="3sp")] ~ meta$X[which(meta$hybrid=="3sp")], col="red")
legend("topleft", col=c("blue", "red"), legend=c("2 species", "3 species"), bty="n", pch=1)

##################
### genetic ID ###
##################
meta$sp.nucl <- NA
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V1 + k9$V2 > 0.6)])] <- "Lasi_flav"
# 4 individuals are in fact myops, they need to be removed
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V2 > 0.6 & k9$V6 > 0.2)])] <- NA
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V3 > 0.6)])] <- "Lasi_plat"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V4 > 0.6)])] <- "Lasi_para"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V5 > 0.6)])] <- "Lasi_nige"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V6 > 0.6)])] <- "Lasi_mixt"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V7 > 0.6)])] <- "Lasi_emar"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V8 > 0.6)])] <- "Lasi_brun"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k9)[which(k9$V9 > 0.6)])] <- "Lasi_fuli"
#meta$sp.nucl[which(meta$CATALOGUENUMBER == "20350")] <- psammophilus

write.table(meta, "metadata_filled.csv", quote = F, row.names = F, sep=";")

################
### 2nd step ###
################

meta$retained_for_step_2 <- F
meta$retained_for_step_2[which(is.na(meta$sp.nucl) == F)] <- T

Lasius.pop <- data.frame(CATALOGUENUMBER=meta$CATALOGUENUMBER[which(meta$retained_for_step_2 == T)], pop = rep("-", sum(meta$retained_for_step_2)))
#Lasius.pop <- read.table("2nd_admixture/2nd_Lasius.pop")
names(Lasius.pop) <- c("CATALOGUENUMBER", "pop")
nref = 5
cutoff = 0.9995
pure <- row.names(k9[which(apply(k9, 1, max) > cutoff),])

# for each species, take the nref individuals with lowest missing data among
for (i in (unique(meta$sp.nucl[which(is.na(meta$sp.nucl)==F)]))) {
  Lasius.pop$pop[which(Lasius.pop$CATALOGUENUMBER %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)][order(meta$percent_missing[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)])][1:nref])] <- i
}

write.table(Lasius.pop, "Lasius.pop", row.names = F, col.names = F, quote = F)


##########################
### Hybridization rate ###
##########################

### Choose best K again

CV2 <- read.table("2nd_admixture/CV_results.txt")
which(CV2$V4 == min(CV2$V4))
# the replicate with the lowest CV error was replicate 6

k8_2 <- read.table("2nd_admixture/rep6/cleaned.8.Q")
dimnames(k8_2)[[1]] <- Lasius.pop$CATALOGUENUMBER


# I'm naming the columns of the admixture results
# based on the species ID of the individual with the maximum assignment to that column
# admixture is a dataframe output by admixture (*.Q) , with individuals as rows, clusters as columns
# meta is the metadata. It must have a "CATALOGUENUMBER" field with individual ID 
# and a "sp.nucl" field with species identification.

rename.admixture.columns <- function(admixture, meta) {
  for (i in 1:length(colnames(admixture))) {
    colnames(admixture)[i] <- meta$sp.nucl[which(meta$CATALOGUENUMBER == row.names(admixture)[which(admixture[,i] == max(admixture[,i]))[1]])]
  }
  return(admixture)
}

k8_renamed <- rename.admixture.columns(k8_2, meta)

barplot(t(as.matrix(order.structure(k8_renamed))), col=c("#1B75BC", "#8DC63F", "#2B3990", "#A97C50", "#F15A29", "#006838", "#F9ED32", "#EC008C"), las = 2, cex.names = 0.1, space = 0, border = NA)

# select only hybrids
hybrids2nd <- which(apply(k8_renamed, 1, max) < 0.999)

barplot(t(as.matrix(order.structure(k8_renamed[hybrids2nd,]))), col=c("#1B75BC", "#8DC63F", "#2B3990", "#A97C50", "#F15A29", "#006838", "#F9ED32", "#EC008C"), las = 2, cex.names = 1, space = 0, border = NA)

######################################################
### Variation in admixture proportions across runs ###
######################################################

# Create a new dataframe with all individuals and all runs
Q_replicates <- data.frame(ind = Lasius.pop$CATALOGUENUMBER)

# iterate across the 20 admixture replicates
for (rep in 1:20){
  # build the path to the Q matrix for one replicate admixture run
  path <- paste0("2nd_admixture/rep", rep, "/cleaned.8.Q")
  # load the corresponding Q matrix
  Q <- as.matrix(read.table(path))
  # get the major cluster assignment for each individual
  this_Q <- apply(Q, 1, max)
  # add it to the main dataframe
  Q_replicates[rep+1] <- this_Q
  # add the replicate number as column name
  names(Q_replicates)[rep+1] <- paste0("rep", rep)
}

# subset only hybrids
# this assumes that any individual with 
Q_hybrids <- subset(Q_replicates, Q_replicates$rep6 < max(Q_replicates$rep6))
Q_hybrids <- Q_hybrids[order(Q_hybrids$rep6),]
Q_hybrids$mean <- apply(Q_hybrids[,2:21], 1, mean)

# Plot the individual admixture proportion for the selected replicate
plot(1-Q_hybrids$rep6,
     log = "y",
     las = 1, 
     col = "#00A19A",
     pch = 20,
     cex=2,
     ylab="Individual admixture proportion")

# add points for the other replicates
for (i in c(2:6, 14:21)) {
  points(1-Q_hybrids[,i],
         pch=20,
         col="#0000001A")
}

abline(h = 2e-4, col="red")


### Based on the above, we "correct" all values 
### between 1e-5 (the minimum) and 2e-4 and assign them to 1e-5

tmp <- as.matrix(k8_renamed)
tmp[which(tmp < 2e-4)] <- 1e-5
# In Lasius, 0.99993 is the maximum value 
# (because all K-1 other clusters have 1e-5)
tmp[which(tmp > 1-2e-4)] <- 0.99993
k8_renamed <- as.data.frame(tmp)

################
### Figure 1 ###
################

## Step 1: Figure out amount of introgression per species

# create vector for introgression amount
intr_amount <- numeric(dim(k8.renamed)[2])
names(intr_amount) <- dimnames(k8.renamed)[[2]]

# Get maximum contribution of each species
max_admixture <- apply(k8.renamed, 1, max)

# For each species i
for (i in 1:dim(k8.renamed)[2]) {
  # Get which individuals belong to this species
  ismax <- k8.renamed[,i] == max_admixture
  
  # estimate the amout of introgression in that species as
  # the sum of the values of other columns for individuals of that species
  # divided by the number of individuals of that species
  # i.e. average amount of introgression in the species
  intr_amount[i] <- sum(k8.renamed[ismax,-i]) / (sum(ismax))
}

sort(intr_amount)

## Step 2: plot "own species" in one colour
## and "introgressed" (all species mixed) in another
# 
# # assign a species name to each sample
# whichismax <- function (x) which(x == max(x, na.rm = T))
# species_column <- apply(k8.renamed, 1, whichismax)
# species_ID = colnames(k8.renamed)[species_column]


# create new df with value for "own species" and "introgressed"
for (i in 1:dim(k8.renamed)[[1]]){
  newadmix$introgressed[i] <- sum(k8.renamed[i,]) - max(k8.renamed[i,])
  newadmix$own[i] <- max(k8.renamed[i,])
}


# Compute maximum assignment probability for all individuals
max_Q <- apply(k8.renamed, 1, max)

# Identify the majority cluster name for each individual
major_cluster <- names(intr_amount)[apply(k8.renamed, 1, which.max)]

## Compute sum of probabilities for all non-majority clusters
sum_other <- rowSums(k8.renamed) - max_Q

## Build new data frame 
Q.for.figure <- data.frame(
  introgressed = sum_other,
  pure = max_Q,
  majority_cluster = major_cluster
)

# Set order of individuals based on species (in order of increasing introgression)
# then by increasing amount of introgression within species
priority_map <- setNames(rank(intr_amount, ties.method = "first"), names(intr_amount))

## ---- 7. Order individuals by:
##      (a) cluster in the order of intr_amount
##      (b) within cluster, decreasing pure
order_idx <- order(match(Q.for.figure$majority_cluster, names(intr_amount)),-Q.for.figure$pure)
order_idx <- order(priority_map[Q.for.figure$majority_cluster], -Q.for.figure$pure)

Q.for.figure <- Q.for.figure[order_idx, ]

Q.for.figure$pure2 <- 0
Q.for.figure$pure2[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10)])] <- Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10,12)])]
Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10)])] <- 0


## Step 3: Plot
par(mfrow=c(1,1), mai=c(0.1, 0.5, 0.1, 0.1))
barplot(t(as.matrix(Q.for.figure[,c(1,2,4)])), las = 2, cex.names = 1, space = 0, border = NA, col=c("#00A19A", "grey80", "grey90"))


####################################
### OLD: Investigating 3-hybrids ###
####################################

### Do 3-hybrids have private alleles from all 3 species ?

# This function returns 1 if two groups of individuals (provided) are fixed for opposite alleles;
# it returns 0 otherwise.
# gt.line is one line of a genotype file (obtained with extract.gt(vcf, "GT"))
# one is a vector containing the positions of the columns of a species of interest
# two is a vector containing the positions of the columns of the species for comparison
# fixed.opposite() returns 1 if alternative alleles are fixed in ones and in two, 0 otherwise.

fixed.opposite <- function(gt.line, one, two) {
  fixed.opposite <- 0
  if (sum(!is.na(gt.line[one]))>0 & sum(!is.na(gt.line[two]))>0) {
    ones <- unique(gt.line[one][which(!is.na(gt.line[one]))])
    if (sum(ones %in% c("0/0", "1/1")) == 1) {
      if (sum(!is.na(gt.line[two]))>0) {
        twos <- unique(gt.line[two][which(!is.na(gt.line[two]))])
        if (length(unique(twos))==1 & length(unique(ones))==1 ) {
          if (sum(twos %in% c("0/0", "1/1") == 1) & ones != twos) fixed.opposite <- 1
        }
      }
    }
  }
  return(fixed.opposite)
}

pure_emar <- which(meta$sp.nucl == "Lasi_emar" & meta$CATALOGUENUMBER %in% row.names(k8.renamed[which(k8.renamed$Lasi_emar > 0.999),]))
pure_nige <- which(meta$sp.nucl == "Lasi_nige" & meta$CATALOGUENUMBER %in% row.names(k8.renamed[which(k8.renamed$Lasi_nige > 0.999),]))
pure_plat <- which(meta$sp.nucl == "Lasi_plat" & meta$CATALOGUENUMBER %in% row.names(k8.renamed[which(k8.renamed$Lasi_plat > 0.999),]))
# Checking how many loci are fixed for opposite alleles in species pairs in each dataset
sum(FO_emar_plat <- apply(gt, 1, fixed.opposite, pure_emar, pure_plat))
sum(FO_emar_nige <- apply(gt, 1, fixed.opposite, pure_emar, pure_nige))
sum(FO_nige_plat <- apply(gt, 1, fixed.opposite, pure_nige, pure_plat))

# Much more clever: private alleles to species A = fixed opposites in A vs (B+C)
sum(private_emar <- apply(gt, 1, fixed.opposite, pure_emar, c(pure_nige, pure_plat)))
sum(private_nige <- apply(gt, 1, fixed.opposite, pure_nige, c(pure_emar, pure_plat)))
sum(private_plat <- apply(gt, 1, fixed.opposite, pure_plat, c(pure_nige, pure_emar)))

# get rid of the flavus/brunneus hybrids
hybrids.emar_nige_plat <- hybrids2nd[- which(names(hybrids2nd) %in% c("10247", "9990815", "9991373"))]

emar_private.gt <- gt[which(private_emar == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))]
nige_private.gt <- gt[which(private_nige == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))]
plat_private.gt <- gt[which(private_plat == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))]

# heterozygous = has a private allele from the corresponding species

N_emar_alleles <- apply(emar_private.gt, 2, count.het)
N_nige_alleles <- apply(nige_private.gt, 2, count.het)
N_plat_alleles <- apply(plat_private.gt, 2, count.het)

par(mfrow=c(4,1), mar = c(4, 5, 0, 0), oma=c(0.5,0.5,0.1,0.1))
barplot(meta$percent_missing[match(row.names(order.structure(k8_renamed[hybrids.emar_nige_plat,])), meta$CATALOGUENUMBER)], las=1, ylab="Missing data", border=NA)
barplot(meta$het[match(row.names(order.structure(k8_renamed[hybrids.emar_nige_plat,])), meta$CATALOGUENUMBER)], las=1, ylab="Relative heterozygosity\n", border=NA)
barplot(rbind(N_emar_alleles, N_nige_alleles, N_plat_alleles)[,match(row.names(order.structure(k8_renamed[hybrids.emar_nige_plat,])), names(N_emar_alleles))], beside = F, col = c("#8DC63F", "#EC008C", "#A97C50"), las=2, ylab="Number of species-specific\nprivate alleles")
legend("topright", pch=15, col=c("#8DC63F", "#EC008C", "#A97C50"), legend=c("L. emarginatus: 75 private alleles", "L. niger: 51 private alleles", "L. platythorax: 22 private alleles"), bty="n")
barplot(t(as.matrix(order.structure(k8_renamed[hybrids.emar_nige_plat,]))), col=c("#8DC63F", "#2B3990", "#F15A29", "#A97C50", "#006838", "#1B75BC", "#EC008C", "#F9ED32"), las = 2, cex.names = 1, space = 0, border = NA)

par(mfrow=c(3,1), mar = c(4, 5, 0, 0), oma=c(0.5,0.5,0.1,0.1))
barplot(meta$percent_missing[match(row.names(order.structure(k8_renamed[hybrids.emar_nige_plat,])), meta$CATALOGUENUMBER)], las=1, ylab="Missing data", border=NA)
barplot(rbind(N_emar_alleles, N_nige_alleles, N_plat_alleles)[,match(row.names(order.structure(k8_renamed[hybrids.emar_nige_plat,])), names(N_emar_alleles))], beside = F, col = c("#8DC63F", "#EC008C", "#A97C50"), las=2, ylab="Number of species-specific\nprivate alleles")
legend("topright", pch=15, col=c("#8DC63F", "#EC008C", "#A97C50"), legend=c("L. emarginatus: 75 private alleles", "L. niger: 51 private alleles", "L. platythorax: 22 private alleles"), bty="n")
barplot(t(as.matrix(order.structure(k8_renamed[hybrids.emar_nige_plat,]))), col=c("#8DC63F", "#2B3990", "#F15A29", "#A97C50", "#006838", "#1B75BC", "#EC008C", "#F9ED32"), las = 2, cex.names = 1, space = 0, border = NA)



## Do these private alleles look like they could have been contaminated?
par(mfrow=c(3,2))
hist(AD.ratio[which(private_emar == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))], main="Triple hybrids", xlab="Distribution of ADR at private alleles", ylab="", col="#8DC63F", las=1)
hist(AD.ratio[which(private_emar == 1), pure_emar], main="Pure individuals", xlab="Distribution of ADR at private alleles", ylab="", col="#8DC63F", las=1)

hist(AD.ratio[which(private_nige == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))], main="", xlab="Distribution of ADR at private alleles", ylab="", col="#EC008C", las=1)
hist(AD.ratio[which(private_nige == 1), pure_nige], main="", xlab="Distribution of ADR at private alleles", ylab="", col="#EC008C", las=1)

hist(AD.ratio[which(private_plat == 1), which(dimnames(gt)[[2]] %in% names(hybrids.emar_nige_plat))], main="", xlab="Distribution of ADR at private alleles", ylab="", col="#A97C50", las=1)
hist(AD.ratio[which(private_plat == 1), pure_plat], main="", xlab="Distribution of ADR at private alleles", ylab="", col="#A97C50", las=1)

###########################
### Hybridization rates ###
###########################

# I compute hybridization rate directionally, from sp2 into sp1.
# as the sum of admixture proportions from sp2 into sp1 divided by the number of individuals in sp1.
# This basically measures the average contribution of sp2 into sp1.
hybridization.rates <- function(admixture_renamed) {
  species <- colnames(admixture_renamed)
  hyb <- matrix(ncol = dim(admixture_renamed)[2], nrow=dim(admixture_renamed)[2], dimnames = list(species, species))
  
  for (sp1 in 1:length(species)) {
    individuals.sp1 <- which(admixture_renamed[,sp1] == apply(admixture_renamed, 1, max))
    for (sp2 in (1:length(species))[-sp1]) {
      # row = sp1, column = sp2
      hyb[sp1, sp2] <- sum(admixture_renamed[individuals.sp1,sp2]) / length(individuals.sp1)
    }
  }
  return(hyb)
}

hybridization_rate.matrix <- hybridization.rates(k8_renamed)
write.table(hybridization_rate.matrix, "Lasius_hybridization_rates.txt", quote = F)

#########################################
### Genetic distances between species ###
#########################################

# This function returns a symmetrical matrix of N x N dimension, where N is the number of species.
# Each cell contains the average distance (same scale as the tree) between individuals of the two corresponding species.
# tree is an object of class = "phylo".
# meta is a dataframe containing at least two columns: 
#   - CATALOGUENUMBER is the name of each individual, with values matching tree$tup.label
#   - sp.nucl is the species of the individual, with N possible values.
get.distances.between.species <- function(tree, meta) {
  spp <- names(table(meta$sp.nucl))
  distances.between.species <- matrix(data = NA, nrow = length(spp), ncol = length(spp), dimnames = list(spp, spp))
  dist.all <- cophenetic(tree)
  for (sp1 in spp) {
    i = which(spp == sp1)
    tips.sp1 <- which(tree$tip.label %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == sp1)])
    for (sp2 in spp) {
      j = which(spp == sp2)
      tips.sp2 <- which(tree$tip.label %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == sp2)])
      if (length(tips.sp1) == 0 | length(tips.sp2) == 0) {} else {
        distances.between.species[i,j] <- mean(dist.all[tips.sp1, tips.sp2])
      }
    }
  }
  diag(distances.between.species) <- NA
  return(distances.between.species)
}

COI_distances <- get.distances.between.species(tree, meta)
COI_distance.reordered <- COI_distances[match(row.names(hybridization_rate.matrix), row.names(COI_distances)), match(colnames(hybridization_rate.matrix), colnames(COI_distances))]
write.table(COI_distance.reordered, "Lasius_COI_distances.txt", quote = F)


########################
### Final formatting ###
########################

matrix.to.df <- function(hybridization_rate.matrix, divergence.matrix) {
  df <- data.frame(sp1 = expand.grid(row.names(hybridization_rate.matrix), colnames(hybridization_rate.matrix))[,2], sp2=expand.grid(row.names(hybridization_rate.matrix), colnames(hybridization_rate.matrix))[,1], hybridization = NA, divergence = NA)
  for (i in 1:dim(df)[1]) {
    df$hybridization[i] <- hybridization_rate.matrix[which(row.names(hybridization_rate.matrix) == df$sp1[i]), which(colnames(hybridization_rate.matrix) == df$sp2[i])]
    df$divergence[i] <- divergence.matrix[which(row.names(hybridization_rate.matrix) == df$sp1[i]), which(colnames(hybridization_rate.matrix) == df$sp2[i])]
  }
  return(df)
}


hybridization.df <- matrix.to.df(hybridization_rate.matrix, COI_distance.reordered)

plot(hybridization.df$hybridization ~ -hybridization.df$divergence, las=1, xlab="COI divergence", ylab="Admixture proportion inherited from other species", pch=19, col=rgb(0,0,0,0.2))
write.table(hybridization.df, "Lasius_hybridization_df.txt", row.names = F, quote = F)




#########################
### Species abundance ###
#########################
table(meta$sp.nucl)
write.table(table(meta$sp.nucl), "spp_abundance.txt", sep = "\t", quote = F)
###



####################
### Co-occurrence ###
####################

co_occurrence <- function(sp, BDM) {
  spp <- names(table(sp))
  co_occurrence <- matrix(data = NA, nrow = length(spp), ncol = length(spp), dimnames = list(spp, spp))
  for (sp1 in spp) {
    i = which(spp == sp1)
    BDM1 <- names(table(meta$BDM[which(meta$sp.nucl == sp1)]))
    for (sp2 in spp) {
      j = which(spp == sp2)
      BDM2 <- names(table(meta$BDM[which(meta$sp.nucl == sp2)]))
      co_occurrence[i,j] <- sum(BDM2 %in% BDM1) / length(BDM1)
    }
  }
  diag(co_occurrence) <- NA
  return(co_occurrence)
}

co_occ <- co_occurrence(meta$sp.nucl, meta$BDM)

write.table(co_occ, "Lasius_co_occurence.txt", row.names = F, quote = F)



###################################
### Pi with and without hybrids ###
###################################

# Defining different thresholds for hybrid exclusion.

introgression.thresholds <- c(0.5, 0.9, 0.99, 0.9999)

# Get amount of "pure" genome per individual
maxQ <- apply(k8_2, 1, max, na.rm=T)

# make directory to store the lists if it doesn't already exist
if (file.exists("maxQ_lists") == F) dir.create("maxQ_lists")

# Making lists of individuals with introgression proportions smaller than the thresholds
# Version with only list of individuals per species
for (threshold in introgression.thresholds) {
  indiv <- names(maxQ)[maxQ > threshold]
  popmap <- data.frame(indiv = indiv, sp = meta$sp.nucl[match(indiv, meta$CATALOGUENUMBER)])
  for (sp in levels(as.factor(popmap$sp))) {
    write(popmap$indiv[which(popmap$sp == sp)], file = paste0("maxQ_lists/", sp, "-maxQlist-", threshold, ".txt"))
  }
}



####################
### Co-occurrence ###
####################
# Large scale
co_occurrence_BDM <- function(sp) {
  spp <- names(table(sp))
  co_occurrence <- matrix(data = NA, nrow = length(spp), ncol = length(spp), dimnames = list(spp, spp))
  for (sp1 in spp) {
    i = which(spp == sp1)
    BDM1 <- names(table(meta$BDM[which(meta$sp.nucl == sp1)]))
    for (sp2 in spp) {
      j = which(spp == sp2)
      BDM2 <- names(table(meta$BDM[which(meta$sp.nucl == sp2)]))
      co_occurrence[i,j] <- sum(BDM2 %in% BDM1) / length(BDM1)
    }
  }
  diag(co_occurrence) <- NA
  return(co_occurrence)
}

co_occ_BDM <- co_occurrence_BDM(meta$sp.nucl)
write.table(co_occ_BDM, "Lasius_co_occurence_BDM.txt", row.names = F, quote = F)


# Fine scale

co_occurrence_plot <- function(sp) {
  spp <- names(table(sp))
  co_occurrence <- matrix(data = NA, nrow = length(spp), ncol = length(spp), dimnames = list(spp, spp))
  for (sp1 in spp) {
    i = which(spp == sp1)
    plot1 <- names(table(meta$PLOT[which(meta$sp.nucl == sp1)]))
    for (sp2 in spp) {
      j = which(spp == sp2)
      plot2 <- names(table(meta$PLOT[which(meta$sp.nucl == sp2)]))
      co_occurrence[i,j] <- sum(plot2 %in% plot1) / length(plot1)
    }
  }
  diag(co_occurrence) <- NA
  return(co_occurrence)
}


co_occ_plot <- co_occurrence_plot(meta$sp.nucl)
write.table(co_occ_plot, "Lasius_co_occurence_fine_scale.txt", row.names = F, quote = F)

#################################################
### Summary table with thresholds -- Table S1 ###
#################################################

intro_table <- data.frame(sp = sort(unique(Q.for.figure$majority_cluster)))

# Set thresholds as average values for all backcross values
# plus some arbitrary thresholds (10%, 5%, 1%, and 2e-4 which is essentially zero)
thresholds = c(0.5^(2:11), 0.1, 0.05, 0.01, 2e-4)
for (t in thresholds) {
  counts <- tapply(Q.for.figure$introgressed > t,
                   Q.for.figure$majority_cluster,
                   sum, na.rm = TRUE)
  intro_table[[paste0("n_above_", t)]] <- as.vector(counts[as.character(intro_table$sp)])
}

write.table(intro_table, "Lasius_introgression_thresholds_summary.txt", quote = F, row.names = F)



###########################
### Mapping per species ###
###########################

mapping <- read.csv("mappingreads_Lasius.csv", h=T)
mapping$sp <- meta$sp.nucl[match(mapping$ID, meta$CATALOGUENUMBER)]
mapping$x <- NA
mapping$x[which(mapping$sp == "Lasi_nig")] <- 0
mapping$x[which(mapping$sp == "Lasi_fuli")] <- 0.013583517
mapping$x[which(mapping$sp == "Lasi_emar")] <- 0.005107297
mapping$x[which(mapping$sp == "Lasi_mixt")] <- 0.017830347
mapping$x[which(mapping$sp == "Lasi_plat")] <- 0.004152741
mapping$x[which(mapping$sp == "Lasi_para")] <- 0.004547326
mapping$x[which(mapping$sp == "Lasi_brun")] <- 0.011446771
mapping$x[which(mapping$sp == "Lasi_flav")] <- 0.017193665


boxplot(mapping$Map.Perc ~ mapping$sp, las=2, xlab="", ylab = "Percentage of reads mapping")
plot(mapping$Map.Perc ~ jitter(mapping$x, factor = 2), pch=20, col="#00000030", las=1, xlab = "dxy", ylab = "Percentage of reads mapping")
plot(mapping$Map.Perc ~ jitter(mapping$x, factor = 1), pch=20, col="#00000030", las=1, xlab = "dxy", ylab = "Percentage of reads mapping", ylim=c(95, 100))



mapping.summary <- data.frame(sp=names(table(mapping$sp)), genus = "Form", dxy = NA, mapping.mean = NA, mapping.SD = NA)
for (i in mapping.summary$sp) {
  mapping.summary$dxy[which(mapping.summary$sp == i)] = mapping$x[which(mapping$sp == i)][1]
  mapping.summary$mapping.mean[which(mapping.summary$sp == i)] = mean(mapping$Map.Perc[which(mapping$sp == i)], na.rm = T)
  mapping.summary$mapping.SD[which(mapping.summary$sp == i)] = sd(mapping$Map.Perc[which(mapping$sp == i)], na.rm = T)
}

write.table(mapping.summary, "Lasi_mapping_summary.txt", quote = F, sep = "\t")

################################
### Missing data per species ###
################################

species = names(table(meta$sp.nucl))
par(mfrow = c(2, length(species)/2), mai=c(0.8, 0.3, 0.2, 0.1))
for (i in species) {
  hist(meta$percent_missing[which(meta$sp.nucl == i)], xlim=c(0, 1), main=i, ylab = "", xlab = "Proportion of missing data", las=1)
}

sum(meta$percent_missing > 0.5) / length(meta$percent_missing)


Q.for.figure$miss <- meta$percent_missing[match(dimnames(Q.for.figure)[[1]], meta$CATALOGUENUMBER)]

Q.for.figure$col <- NA
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_brun")] <- "#006838"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_emar")] <- "#8DC63F"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_flav")] <- "#F9ED32"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_fuli")] <- "#1B75BC"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_mixt")] <- "#2B3990"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_nige")] <- "#EC008C"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_para")] <- "#F15A29"
Q.for.figure$col[which(Q.for.figure$majority_cluster == "Lasi_plat")] <- "#A97C50"

Q.for.figure$introgressed2 <- Q.for.figure$introgressed
Q.for.figure$introgressed2 [which(Q.for.figure$introgressed == 7e-05)] <- 0
par(mfrow = c(1,1), mai=c(0.8, 0.8, 0.2, 0.1))

plot(Q.for.figure$introgressed ~ Q.for.figure$miss, log="y", pch=20, col=Q.for.figure$col, las = 1, ylab = "Introgression proportion", xlab = "proportion of missing data")
(zi.miss <- glmmTMB(introgressed ~ miss + (1|majority_cluster),
                    ziformula = ~ miss + (1|majority_cluster), 
                    family = beta_family(link = "logit"), 
                    data=Q.for.figure))

summary(zi.miss)


########################
### Maps per species ###
########################
# This part was written in large parts by Claude Opus 5.5

library(terra)
library(sf)
library(elevatr)
library(geodata)

crs_ch <- "EPSG:21781"

# ---- 1) Vaud boundary (GADM) ----
ch   <- geodata::gadm("CHE", level = 1, path = tempdir())
vaud <- ch[ch$NAME_1 == "Vaud", ]

# ---- 1) Elevation model download (AWS Terrain Tiles via elevatr) ----
# z = zoom level: 9 ~ 200 m, 10 ~ 100 m, 11 ~ 50 m pixels (higher = slower)
dem <- elevatr::get_elev_raster(sf::st_as_sf(vaud), z = 10, clip = "bbox")
dem <- rast(dem)

# ---- 2) Reproject to Swiss coordinates ----
dem       <- project(dem, crs_ch, method = "bilinear")
vaud_ch   <- project(vaud, crs_ch)
dem       <- mask(crop(dem, vaud_ch), vaud_ch)

# ---- 1) Hillshade ----
slope  <- terrain(dem, "slope",  unit = "radians")
aspect <- terrain(dem, "aspect", unit = "radians")

# angle = sun elevation (degrees above horizon), direction = sun azimuth (0 = N, 90 = E)
hs <- shade(slope, aspect, angle = 40, direction = 315)

# Optional: multidirectional hillshade (average of several sun azimuths)
# hs <- mean(rast(lapply(c(270, 315, 0, 45),
#                        function(d) shade(slope, aspect, angle = 45, direction = d))))

## Define colour palette and min/max introgression values in the dataset
meta$introgressed <- NA
meta$introgressed[match(row.names(Q.for.figure), meta$CATALOGUENUMBER)] <- Q.for.figure$introgressed

palette <- colorRampPalette(c("grey95", "#00A19A"))(100)

lmin <- log10(min(meta$introgressed, na.rm = TRUE))
lmax <- log10(max(meta$introgressed, na.rm = TRUE))
idx  <- (log10(meta$introgressed) - lmin) / (lmax - lmin)   # 0â€“1


# ---- 3) Plot ----
par(mfrow=c(3,3), mar=c(0.05,0,0,0))
for (sp in unique(Q.for.figure$majority_cluster)) {
  plot(hs, col = grey(0:100 / 100), legend = F, axes = F, main=sp)
  
  # add admixture proportion to metadata
  
  meta$col.for.map <- palette[1 + round(idx * 99)]
  meta$pch.for.map <- ifelse(meta$introgressed > 10^lmin, 21, 4)  # 21 = filled circle, 4 = x
  
  points(meta$X[which(meta$sp.nucl == sp)], meta$Y[which(meta$sp.nucl == sp)], pch = meta$pch.for.map[which(meta$sp.nucl == sp)], bg = meta$col.for.map[which(meta$sp.nucl == sp)], col = 1)
}



#######################
### ADR per species ###
#######################

# Function to apply a sum like pmin or pmax
# Written by Sven Hohenstein here:
# https://stackoverflow.com/questions/13123638/there-is-pmin-and-pmax-each-taking-na-rm-why-no-psum
psum <- function(...,na.rm=FALSE) {
  dat <- do.call(cbind,list(...))
  res <- rowSums(dat, na.rm=na.rm)
  idx_na <- !rowSums(!is.na(dat))
  res[idx_na] <- NA
  res
}


# get allelic depth for the reference allele at heterozygous sites
ad.ref <<- extract.gt(vcf, "AD", as.numeric = T)
ad.ref[which(dp < 8)] <- NA
ad.ref.het <<- ad.ref
ad.ref.het[which(gt %in% c("0/0", "1/1"))] <- NA

# and the alternative allele
ad.alt <<- dp - ad.ref
ad.alt[which(dp < 8)] <- NA
ad.alt.het <<- ad.alt
ad.alt.het[which(gt %in% c("0/0", "1/1"))] <- NA

# compute allelic depth ratio (ADR), defined as
# depth at the allele supported by most reads divided by total depth at the locus.
# takes two matrices with allelic depths (integer) as input
get.ad.ratio <<- function(ad.ref, ad.alt) return(pmax(ad.ref, ad.alt) / mapply(psum, ad.ref, ad.alt, na.rm=T))

# use the function on the whole dataset
AD.ratio <- get.ad.ratio(as.matrix(ad.ref.het), as.matrix(ad.alt.het))

# get positions of heterozygote genotypes
het_genotypes <- which(gt %in% c("0/1", "1/0"))

# subset depth values at heterozygote genotypes only
dp.het <- dp[het_genotypes]

# this function samples a number of reads for one of the alleles
# following a binomial distribution and equal probabilities for both alleles
# given the depth at the locus
sim.ad <- function(depth) {
  rbinom(n = 1, size = depth, prob = 0.5)
}

# probably not the most efficient way but it works.
# applies sim.ad.ref on observed depth values at heterozygous loci
sim.ad.ref <- unlist(lapply(X = dp.het, FUN = sim.ad))

# create a dataframe with simulated AD at both alleles
ad.df <- data.frame(A=sim.ad.ref, B=dp.het - sim.ad.ref)

# compute ADR for the simulated AD values
sim.ad.ratio <- get.ad.ratio(apply(ad.df, 1, max), apply(ad.df, 1, min))


### Loop through species

par(mfrow=c(1,4), mai=c(0.5, 0.1, 0.4, 0.1))

for(sp in unique(sort(Q.for.figure$majority_cluster))) {
  # Getthe list of individuals from this species
  individuals <- dimnames(Q.for.figure)[[1]][which(Q.for.figure$majority_cluster == sp)]
  plot(NULL, xlim=c(0.45, 0.7), ylim=c(0,30), xlab="ADR", ylab="", main=sp, yaxt="n", las=1)
  # For each individual of this species
  for (i in which(dimnames(AD.ratio)[[2]] %in% individuals)) {
    this_individual_name <- dimnames(AD.ratio)[[2]][i]
    # Only if there are enough non-NA ADR values
    if (sum(is.na(AD.ratio[,i]) == F) > 10) {
      col.ADR = "#00000030"
      if (Q.for.figure$introgressed[which(dimnames(Q.for.figure)[[1]] == this_individual_name)] > 2e-04) {
        col.ADR <- "#FF000080"
      }
      points(density(AD.ratio[,i][which(AD.ratio[,i] != 1)], na.rm=T) , type="l", col=col.ADR)
    }
  }
  # Add the simulated distribution
  # points(density(sim.ad.ratio, bw = 0.01) , type="l", lwd=3, col="#00A19A")
}
