library(vcfR)
library(adegenet)
library(phytools)
library(MASS)


prop.het<-function(vector) {
  return(sum(vector=="0/1" | vector=="1/0", na.rm=T)/sum(is.na(vector)==F))
}

count.na<-function(x) sum(is.na(x))

'%!in%' <- Negate('%in%')

###############
### Dataset ###
###############

vcf <- read.vcfR("Campo_1snp_DP8_meanDP200_mac2_miss75_ADR_filtered_discard0.05_correct0.25.vcf.gz")
gt <- extract.gt(vcf, "GT")
dp <- extract.gt(vcf, "DP", as.numeric=T)
gl <- vcfR2genlight(vcf)

individuals <- dimnames(gt)[[2]]

meta.all <- read.csv("Metadata_combined_public_scientific.csv", h=T, sep=";")

meta <- meta.all[match(dimnames(gt)[[2]], meta.all$CATALOGUENUMBER),]

meta$col.morpho <- "grey"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_aeth")] <- "#BE1E2D"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_fall")] <- "#EF3E23"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_herc")] <- "#27AAE1"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_lign")] <- "#F9ED32"
meta$col.morpho[which(meta$SPECIESSUMMARY=="Camp_pice")] <- "#006838"

tree.not.pruned <- read.tree("Camponotus.tree")
tree <- drop.tip(phy = tree.not.pruned, tip = tree.not.pruned$tip.label[tree.not.pruned$tip.label %in% meta$CATALOGUENUMBER == F], trim.internal = T)


legend_Campo <- c("C. ligniperda", "C. herculeanus", "C. aethiops", "C. fallax", "C. piceus")
col_Campo <- c("#F9ED32", "#27AAE1", "#BE1E2D", "#EF3E23", "#006838")

meta$percent_missing <- apply(gt, 2, count.na) / dim(gt)[1]


######################
### Heterozygosity ###
######################

meta$het <- apply(gt, 2, prop.het)

par(mfrow=c(1,1))
plot(meta$het ~ meta$percent_missing, col=meta$col.morpho, pch=19)

###########
### MDS ###
###########
par(mfrow=c(1,1))
dist_all <- dist(gl)
MDS_all <- isoMDS(dist_all, k = 2, maxit = 100)
plot(MDS_all$points[,2] ~ MDS_all$points[,1], col=meta$col.morpho, pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - all markers", cex=6*sqrt(meta$het)[match(dimnames(gt)[[2]],  meta$CATALOGUENUMBER)])
legend("bottomright", legend=legend_Campo, pch=19, col=col_Campo, bty="n", ncol = 1)


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
# we chose K = 4 based on these results
# the replicate with the lowest CV error was replicate 1

k4 <- read.table("admixture/rep1/cleaned.4.Q")
row.names(k4) <- dimnames(gt)[[2]]

par(mfrow=c(3,1), mai=c(0.1, 0.6, 0.1, 0.1))
barplot(meta$het[match(tree$tip.label, individuals)], las=1, ylab="Relative heterozygosity", border=NA)
barplot(t(as.matrix(k4))[,match(tree$tip.label, individuals)], las=2, border=NA, space=0, cex.names=0.6, col=c("#006838", "#F9ED32", "#27AAE1", "#EF3E23"))
plot.phylo(tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
legend("bottomleft", legend=legend_Campo, pch=15, col=col_Campo, bty="n")
tiplabels(pch = 15, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])


##################
### genetic ID ###
##################

meta$sp.nucl <- NA
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k4)[which(k4$V2 > 0.9)])] <- "Camp_lign"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k4)[which(k4$V3 > 0.9)])] <- "Camp_herc"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% rownames(k4)[which(k4$V4 > 0.9)])] <- "Camp_fall"


meta$retained_for_step_2 <- F
meta$retained_for_step_2[which(is.na(meta$sp.nucl) == F)] <- T

Camponotus.pop <- data.frame(CATALOGUENUMBER=meta$CATALOGUENUMBER[which(meta$retained_for_step_2 == T)], pop = rep("-", sum(meta$retained_for_step_2)))

nref = 5
cutoff = 0.9995
pure <- row.names(k5[which(apply(k5, 1, max) > cutoff),])
# for each species, take the nref individuals with lowest missing data among
for (i in (unique(meta$sp.nucl[which(is.na(meta$sp.nucl)==F)]))) {
  Camponotus.pop$pop[which(Camponotus.pop$CATALOGUENUMBER %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)][order(meta$percent_missing[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)])][1:nref])] <- i
}

write.table(Camponotus.pop, "Camponotus.pop", row.names = F, col.names = F, quote = F)
#Camponotus.pop <- read.table("2nd_admixture/2nd_Camponotus.pop")


##########################
### Introgression rate ###
##########################

### Choose best K again

CV2 <- read.table("2nd_admixture/CV_results.txt")
which(CV2$V4 == min(CV2$V4))
# the replicate with the lowest CV error was replicate 8

k3_2 <- read.table("2nd_admixture/rep8/cleaned.3.Q")
rownames(k3_2) <- Camponotus.pop[,1]

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

k3_renamed <- rename.admixture.columns(k3_2, meta)

# This function orders a Q matrix to group samples by species
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
barplot(t(as.matrix(order.structure(k3_renamed))), col=c("#F9ED32", "#EF3E23", "#27AAE1"), las = 2, cex.names = 0.1, space = 0, border = NA)


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

hybridization_rate.matrix <- hybridization.rates(k3_renamed)

write.table(hybridization_rate.matrix, "Camponotus_hybridization_rates.txt", quote = F)



################
### Figure 1 ###
################

## Step 1: Figure out amount of introgression per species

# create vector for introgression amount
intr_amount <- numeric(dim(k3_renamed)[2])
names(intr_amount) <- dimnames(k3_renamed)[[2]]

# Get maximum contribution of each species
max_admixture <- apply(k3_renamed, 1, max)

# For each species i
for (i in 1:dim(k3_renamed)[2]) {
  # Get which individuals belong to this species
  ismax <- k3_renamed[,i] == max_admixture
  
  # estimate the amout of introgression in that species as
  # the sum of the values of other columns for individuals of that species
  # divided by the number of individuals of that species
  # i.e. average amount of introgression in the species
  intr_amount[i] <- sum(k3_renamed[ismax,-i]) / (sum(ismax))
}

sort(intr_amount)

## Step 2: plot "own species" in one colour
## and "introgressed" (all species mixed) in another
# 
# # assign a species name to each sample
# whichismax <- function (x) which(x == max(x, na.rm = T))
# species_column <- apply(k3_renamed, 1, whichismax)
# species_ID = colnames(k3_renamed)[species_column]


# create new df with value for "own species" and "introgressed"
for (i in 1:dim(k3_renamed)[[1]]){
  newadmix$introgressed[i] <- sum(k3_renamed[i,]) - max(k3_renamed[i,])
  newadmix$own[i] <- max(k3_renamed[i,])
}


# Compute maximum assignment probability for all individuals
max_Q <- apply(k3_renamed, 1, max)

# Identify the majority cluster name for each individual
major_cluster <- names(intr_amount)[apply(k3_renamed, 1, which.max)]

## Compute sum of probabilities for all non-majority clusters
sum_other <- rowSums(k3_renamed) - max_Q

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
barplot(t(as.matrix(Q.for.figure[,c(1,2,4)])), las = 2, cex.names = 1, space = 0, border = NA, col=c("#A3195B", "grey80", "grey90"))


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
write.table(COI_distance.reordered, "Camponotus_COI_distances.txt", quote = F)

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

plot(hybridization.df$hybridization ~ -hybridization.df$divergence, las=1, xlab="Divergence time (My)", ylab="Admixture proportion inherited from other species", pch=19, col=rgb(0,0,0,0.2), log="y")
write.table(hybridization.df, "Camponotus_hybridization_df.txt", row.names = F, quote = F)



###################################
### Pi with and without hybrids ###
###################################

# Defining different thresholds for hybrid exclusion.

introgression.thresholds <- c(0.5, 0.9, 0.99, 0.9999)

# Get amount of "pure" genome per individual
maxQ <- apply(k3_2, 1, max, na.rm=T)

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
write.table(co_occ_BDM, "Camponotus_co_occurence_BDM.txt", row.names = F, quote = F)


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
write.table(co_occ_plot, "Camponotus_co_occurence_fine_scale.txt", row.names = F, quote = F)


###########################
### Mapping per species ###
###########################

mapping <- read.csv("mappingreads_Camponotus.csv", h=T)
mapping$sp <- meta$sp.nucl[match(mapping$ID, meta$CATALOGUENUMBER)]
boxplot(mapping$Map.Perc ~ mapping$sp, las=2, xlab="", ylab = "Percentage of reads mapping")
boxplot(mapping$Map.Perc ~ mapping$sp, las=2, xlab="", ylab = "Percentage of reads mapping", ylim=c(80, 100))

mapping$x <- NA
mapping$x[which(mapping$sp == "Camp_fall")] <- 0
mapping$x[which(mapping$sp == "Camp_lign")] <- 0.011404068
mapping$x[which(mapping$sp == "Camp_herc")] <- 0.011510916

plot(mapping$Map.Perc ~ jitter(mapping$x, factor = 2), pch=20, col="#00000030", las=1, xlab = "dxy", ylab = "Percentage of reads mapping")
plot(mapping$Map.Perc ~ jitter(mapping$x, factor = 2), pch=20, col="#00000030", las=1, xlab = "dxy", ylab = "Percentage of reads mapping", ylim=c(80, 100))

mapping.summary <- data.frame(sp=names(table(mapping$sp)), genus = "Camp", dxy = NA, mapping.mean = NA, mapping.SD = NA)
for (i in mapping.summary$sp) {
  mapping.summary$dxy[which(mapping.summary$sp == i)] = mapping$x[which(mapping$sp == i)][1]
  mapping.summary$mapping.mean[which(mapping.summary$sp == i)] = mean(mapping$Map.Perc[which(mapping$sp == i)])
  mapping.summary$mapping.SD[which(mapping.summary$sp == i)] = sd(mapping$Map.Perc[which(mapping$sp == i)])
}

write.table(mapping.summary, "Campo_mapping_summary.txt", quote = F, sep = "\t")


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

palette <- colorRampPalette(c("grey95", "#49196B"))(100)

lmin <- log10(min(meta$introgressed, na.rm = TRUE))
lmax <- log10(max(meta$introgressed, na.rm = TRUE))
idx  <- (log10(meta$introgressed) - lmin) / (lmax - lmin)


# ---- 3) Plot ----
par(mfrow=c(1,3), mar=c(0.05,0,0,0))
for (sp in unique(Q.for.figure$majority_cluster)) {
  plot(hs, col = grey(0:100 / 100), legend = F, axes = F, main=sp)
  
  # add admixture proportion to metadata
  
  meta$col.for.map <- palette[1 + round(idx * 99)]
  meta$pch.for.map <- ifelse(meta$introgressed > 2e-4, 21, 4)  # 21 = filled circle, 4 = x
  
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
  plot(NULL, xlim=c(0.45, 0.7), ylim=c(0,25), xlab="ADR", ylab="", main=sp, yaxt="n", las=1)
  # For each individual of this species
  for (i in which(dimnames(AD.ratio)[[2]] %in% individuals)) {
    this_individual_name <- dimnames(AD.ratio)[[2]][i]
    # Only if there are enough non-NA ADR values
    if (sum(is.na(AD.ratio[,i]) == F) > 10) {
      col.ADR = "#00000030"
      if (Q.for.figure$introgressed[which(dimnames(Q.for.figure)[[1]] == this_individual_name)] > 2e-04) {
        col.ADR <- "red"
      }
      points(density(AD.ratio[,i][which(AD.ratio[,i] != 1)], na.rm=T) , type="l", col=col.ADR)
    }
  }
}
