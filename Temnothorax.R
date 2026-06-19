library(vcfR)
library(adegenet)
library(phytools)
library(MASS)

prop.het<-function(vector) {
  return(sum(vector=="0/1" | vector=="1/0", na.rm=T)/sum(is.na(vector)==F))
}

count.na<-function(x) sum(is.na(x))

'%!in%' <- Negate('%in%')


################
### Datasets ###
################

vcf <- read.vcfR("Temno_DP8_meanDP200_mac2_miss75_ADR_filtered_discard0.05_correct0.25.vcf.gz")
gt <- extract.gt(vcf, "GT")
dp <- extract.gt(vcf, "DP", as.numeric=T)
gl <- vcfR2genlight(vcf)

individuals <- dimnames(gt)[[2]]

meta.all <- read.csv("Metadata_combined_public_scientific.csv", h=T, sep=";")

meta <- meta.all[match(dimnames(gt)[[2]], meta.all$CATALOGUENUMBER),]#[which(is.na(meta$CATALOGUENUMBER)==F),]

meta$col.morpho <- "#A0A0A0"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_affi")] <- "#EC008C"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_nigr")] <- "#3E499E"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_nyla")] <- "#F58122"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_unif")] <- "#74C483"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_cort")] <- "#D2D028"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_parv")] <- "#B92025"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_tube")] <- "#6F3686"
meta$col.morpho[which(meta$SPECIESSUMMARY == "Temn_inte")] <- "#81CEC7"


meta$percent_missing <- apply(gt, 2, count.na) / dim(gt)[1]

legend_Temno <- c("T. affinis", "T. corticalis", "T. interruptus", "T. nigriceps", "T.nylanderi", "T. parvulus", "T. tuberum", "T. unifasciatus")
col_Temno <- c("#EC008C", "#D2D028", "#81CEC7", "#3E499E", "#F58122", "#B92025", "#6F3686", "#74C483")

tree.not.pruned <- read.tree("Temnothorax.tree")
tree <- drop.tip(phy = tree.not.pruned, tip = tree.not.pruned$tip.label[tree.not.pruned$tip.label %in% meta$CATALOGUENUMBER == F], trim.internal = T)

######################
### Heterozygosity ###
######################

meta$het <- apply(gt, 2, prop.het)

par(mfrow=c(1,1))
plot(meta$het ~ meta$percent_missing, col=meta$col.morpho, pch=19)
legend("topright", legend = legend_Temno, pch=15, col=col_Temno, bty="n")


###########
### MDS ###
###########

par(mfrow=c(1,1))
dist_all <- dist(gl)
MDS_all <- isoMDS(d = dist_all, k = 2, maxit = 100)
plot(MDS_all$points[,2] ~ MDS_all$points[,1], col=meta$col.morpho, cex=6*sqrt(meta$het)[match(dimnames(gt)[[2]],  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - all chromosomes")
legend("topright", legend=legend_Temno, pch=15, col=col_Temno, bty="n", ncol = 2)



#################
### admixture ###
#################

### Choose best K

CV <- read.table("admixture/CV_summary_across_rep.txt", h=F)
CV <- CV[order(CV$V1),]
dimnames(CV)[[1]] <- 1:10
K <- dim(CV)[2] - 1

plot(NULL, xlim=range(CV$V1), ylim=range(CV[,2:dim(CV)[2]]), las=1, ylab="Cross-Validation error", xlab = "K")
for (i in 2:dim(CV)[2]) {
  points(CV[,i] ~ CV$V1, type="b", pch=20, col=rgb(0, 0, 0, 0.2))
}

CV.mean <- apply(CV[,-1], 1, mean)
CV.sd <- apply(CV[,-1], 1, sd)
plot(CV.mean, pch=20, cex=0.75, col=NULL, las=1)
points(CV.mean, type="l", col="lightgrey")
for (i in 1:K) segments(x0=i, x1=i, y0=CV.mean[i] - CV.sd[i], y1=CV.mean[i] + CV.sd[i])
# we chose K = 6 based on these results
# the replicate with the lowest CV error was replicate 7

k6 <- read.table("admixture/rep7/cleaned.6.Q")
row.names(k6) <- dimnames(gt)[[2]]

par(mfrow=c(2,1), mar=c(3, 4, 0, 0))
barplot(meta$het[order(meta$col.morpho)], las=2, border=F, cex.names = 0, col=meta$col.morpho[order(meta$col.morpho)], ylab="Relative heterozygosity")
barplot(t(as.matrix(k6[order(meta$col.morpho),])), las=2, col=c("#EC008C", "#74C483", "#3E499E", "#F58122", "black", "#81CEC7"), border=NA, space=0, cex.names=0.6)

par(mfrow=c(2,1), mai=c(0.1, 0.6, 0.1, 0.1))
barplot(t(as.matrix(k6))[,match(tree$tip.label, individuals)], las=2, border=NA, space=0, cex.names=0.6, col=c("#EC008C", "#74C483", "#3E499E", "#F58122", "black", "#81CEC7"))
plot.phylo(tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
legend("bottomleft", legend=legend_Temno, pch=15, col=col_Temno, bty="n")
tiplabels(pch = 15, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])

##################
### genetic ID ###
##################

meta$sp.nucl <- NA
meta$sp.nucl[which(MDS_all$points[,1] > 20)] <- "Temn_nyla"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k6)[which(k6$V2 > 0.6)])] <- "Temn_unif"
meta$sp.nucl[which(MDS_all$points[,2] > 20)] <- "Temn_affi"
meta$sp.nucl[which(MDS_all$points[,2] > 2.5 & MDS_all$points[,2] < 10 & MDS_all$points[,1] < -20)] <- "Temn_nigr"


################
### 2nd step ###
################

meta$retained_for_step_2 <- F
meta$retained_for_step_2[which(is.na(meta$sp.nucl) == F)] <- T

Temnothorax.pop <- data.frame(CATALOGUENUMBER=meta$CATALOGUENUMBER[which(meta$retained_for_step_2 == T)], pop = rep("-", sum(meta$retained_for_step_2)))

nref = 5
cutoff = 0.9995
pure <- row.names(k6[which(apply(k6, 1, max) > cutoff),])
# for each species, take the nref individuals with lowest missing data among
for (i in (unique(meta$sp.nucl[which(is.na(meta$sp.nucl)==F)]))) {
  Temnothorax.pop$pop[which(Temnothorax.pop$CATALOGUENUMBER %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)][order(meta$percent_missing[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)])][1:nref])] <- i
}

write.table(Temnothorax.pop, "Temnothorax.pop", row.names = F, col.names = F, quote = F)



##########################
### Hybridization rate ###
##########################

### Choose best K again

CV2 <- read.table("2nd_admixture/CV_results.txt")
which(CV2$V4 == min(CV2$V4))
# the replicate with the lowest CV error was replicate 9

k4_2 <- read.table("2nd_admixture/rep9/cleaned.4.Q")
dimnames(k4_2)[[1]] <- Temnothorax.pop$CATALOGUENUMBER


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

k4_renamed <- rename.admixture.columns(k4_2, meta)


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

par(mfrow=c(1,1))
barplot(t(as.matrix(order.structure(k4_renamed))), col=c("#74C483", "#F58122", "#3E499E", "#EC008C"), las = 2, cex.names = 0.1, space = 0, border = NA)


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

hybridization_rate.matrix <- hybridization.rates(k4_renamed)
write.table(hybridization_rate.matrix, "Temnothorax_hybridization_rates.txt", quote = F)



################
### Figure 1 ###
################

## Step 1: Figure out amount of introgression per species

# create vector for introgression amount
intr_amount <- numeric(dim(k4_renamed)[2])
names(intr_amount) <- dimnames(k4_renamed)[[2]]

# Get maximum contribution of each species
max_admixture <- apply(k4_renamed, 1, max)

# For each species i
for (i in 1:dim(k4_renamed)[2]) {
  # Get which individuals belong to this species
  ismax <- k4_renamed[,i] == max_admixture
  
  # estimate the amout of introgression in that species as
  # the sum of the values of other columns for individuals of that species
  # divided by the number of individuals of that species
  # i.e. average amount of introgression in the species
  intr_amount[i] <- sum(k4_renamed[ismax,-i]) / (sum(ismax))
}

sort(intr_amount)

## Step 2: plot "own species" in one colour
## and "introgressed" (all species mixed) in another
# 
# # assign a species name to each sample
# whichismax <- function (x) which(x == max(x, na.rm = T))
# species_column <- apply(k4_renamed, 1, whichismax)
# species_ID = colnames(k4_renamed)[species_column]


# create new df with value for "own species" and "introgressed"
for (i in 1:dim(k4_renamed)[[1]]){
  newadmix$introgressed[i] <- sum(k4_renamed[i,]) - max(k4_renamed[i,])
  newadmix$own[i] <- max(k4_renamed[i,])
}


# Compute maximum assignment probability for all individuals
max_Q <- apply(k4_renamed, 1, max)

# Identify the majority cluster name for each individual
major_cluster <- names(intr_amount)[apply(k4_renamed, 1, which.max)]

## Compute sum of probabilities for all non-majority clusters
sum_other <- rowSums(k4_renamed) - max_Q

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
Q.for.figure$pure2[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2)])] <- Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2)])]
Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2)])] <- 0


## Step 3: Plot
par(mfrow=c(1,1), mai=c(0.1, 0.5, 0.1, 0.1))
barplot(t(as.matrix(Q.for.figure[,c(1,2,4)])), las = 2, cex.names = 1, space = 0, border = NA, col=c("#95C11F", "grey80", "grey90"))



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
write.table(COI_distance.reordered, "Temnothorax_COI_distances.txt", quote = F)


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

plot(hybridization.df$hybridization ~ -hybridization.df$divergence, las=1, xlab="Divergence time (My)", ylab="Admixture proportion inherited from other species", pch=19, col=rgb(0,0,0,0.2))
write.table(hybridization.df, "Temnothorax_hybridization_df.txt", row.names = F, quote = F)



###################################
### Pi with and without hybrids ###
###################################

# Defining different thresholds for hybrid exclusion.

introgression.thresholds <- c(0.5, 0.9, 0.99, 0.9999)

# Get amount of "pure" genome per individual
maxQ <- apply(k4_2, 1, max, na.rm=T)

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
write.table(co_occ_BDM, "Temnothorax_co_occurence_BDM.txt", row.names = F, quote = F)


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
write.table(co_occ_plot, "Temnothorax_co_occurence_fine_scale.txt", row.names = F, quote = F)
