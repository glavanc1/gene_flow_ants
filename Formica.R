library(vcfR)
library(adegenet)
library(phytools)

prop.het<-function(vector) {
  return(sum(vector=="0/1" | vector=="1/0", na.rm=T)/sum(is.na(vector)==F))
}

count.na<-function(x) sum(is.na(x))

'%!in%' <- Negate('%in%')

meta.all <- read.csv("Metadata_combined_public_scientific.csv", h=T, sep=";")
meta.all$col.morpho <- "grey"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_cine")] <- "#6D6E71"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_sely")] <- "#2B3990"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_fusco")] <- "#27AAE1"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_fusca")] <- "#603913"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_lema")] <- "#1B75BC"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_cuni")] <- "#A97C50"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_rufi")] <- "#00A19A"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_pice")] <- "#009640"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_pres")] <- "#8DC63F"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_exse")] <- "#006633"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_trun")] <- "#E6007E"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY %in% c("Form_lugu/para", "Form_lugubris/paralugubris"))] <- "#29235C"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_rufa")] <- "#D45E4D"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_poly")] <- "#CBBBA0"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_prat")] <- "#D7DF23"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_sang")] <- "#BE1E2D"
meta.all$col.morpho[which(meta.all$SPECIESSUMMARY=="Form_aqui")] <- "#FBB040"


vcf <- read.vcfR("Formica_1snp_DP8_meanDP200_mac2_miss75_ADR_filtered_discard0.05_correct0.25.vcf")

gt <- extract.gt(vcf, "GT")
dp <- extract.gt(vcf, "DP", as.numeric=T)
gl <- vcfR2genlight(vcf)

individuals <- dimnames(gt)[[2]]

meta <- meta.all[match(dimnames(gt)[[2]], meta.all$CATALOGUENUMBER),]
meta$SPECIESIDFINAL[which(meta$SPECIESIDFINAL %in% c("", "Myrm_lobi"))] <- NA
meta$SPECIESSUMMARY[which(meta$SPECIESSUMMARY == "Myrm_lobi")] <- "Form"

tree.not.pruned <- read.tree("Formica.tree")

legend_Formica <- c("Raptiformica sanguinea", "Formica pratensis", "Formica polyctena", "Formica rufa", "Formica (para)lugubris", "Formica aquilonia", "Formica truncorum", "Coptoformica exsecta", "Coptoformica pressilabris", "Serviformica picea", "Serviformica rufibarbis", "Serviformica cunicularia", "Serviformica lemani", "Serviformica fusca", "Serviformica fuscocinerea", "Serviformica cinerea", "Serviformica selysi")
col_Formica <- c("#BE1E2D", "#D7DF23", "#CBBBA0", "#D45E4D", "#29235C", "#FBB040", "#E6007E", "#006633", "#8DC63F", "#009640", "#00A19A", "#A97C50", "#1B75BC", "#603913", "#27AAE1", "#6D6E71", "#2B3990")

# Compute relative heterozygosity (proportion of heterozygous SNPs)
meta$het <- apply(gt, 2, prop.het)

meta$percent_missing <- apply(gt, 2, count.na) / dim(gt)[1]

#################################################
### Identify rufibarbis individuals with NUMT ###
#################################################

### Function written by Liam Revell (author of phytools):
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

plot(tree.not.pruned)
edgelabels()
# the NUMTs all derive from edge 3548
all.NUMT.nodes <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = tree.not.pruned$edge[3548,2])]
NUMT.tips <- all.NUMT.nodes[which(is.na(all.NUMT.nodes)==F)]

plot(drop.tip(tree.not.pruned, tip = tree.not.pruned$tip.label[which(tree.not.pruned$tip.label %in% NUMT.tips)]))

NUMT.tips <- read.table("rufibarbis_NUMT.txt")$V1
tree_noNUMT <- drop.tip(phy = tree.not.pruned, tip = tree.not.pruned$tip.label[tree.not.pruned$tip.label %in% NUMT.tips], trim.internal = T)

tree <- drop.tip(phy = tree_noNUMT, tip = tree_noNUMT$tip.label[tree_noNUMT$tip.label %in% meta$CATALOGUENUMBER == F], trim.internal = T)


###########
### MDS ###
###########

par(mfrow=c(1,1))
dist_all <- dist(gl)
MDS_all <- isoMDS(d = dist_all, k = 2, maxit = 100)
plot(MDS_all$points[,2] ~ MDS_all$points[,1], col=meta$col.morpho, pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - all markers", cex=10 * sqrt(meta$het))
legend("topright", legend=legend_Formica, pch=15, col=col_Formica, bty="n", ncol = 2)

# With Formica sensu stricto subset
# Subset the vcf
# the first column of the vcf is the format, so we add it and +1 to all column numbers
str.vcf <- vcf[,c(1, which(MDS_all$points[,1] > 10)+1)]
str.gl <- vcfR2genlight(str.vcf)
dist_str <- dist(str.gl)
MDS_str <- isoMDS(d = dist_str, k = 2, maxit = 100)
plot(MDS_str$points[,2] ~ MDS_str$points[,1], col=meta$col.morpho[match(row.names(MDS_str$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Formica s. str.", cex=80 * meta$het[match(row.names(MDS_str$points),  meta$CATALOGUENUMBER)])

sstr <- individuals[which(MDS_all$points[,1] > 10)]
write(sstr, "Formica_s_str.txt")
notsstr <- individuals[which(MDS_all$points[,1] < 10)]
write(notsstr, "Formica_no_s_str.txt")

# With Coptoformica
copt.vcf <- vcf[,c(1, which(MDS_all$points[,1] > -10 & MDS_all$points[,2] < -10)+1)]
copt.gl <- vcfR2genlight(copt.vcf)
dist_copt <- dist(copt.gl)
MDS_copt <- isoMDS(d = dist_copt, k = 2, maxit = 100)
plot(MDS_copt$points[,2] ~ MDS_copt$points[,1], col=meta$col.morpho[match(row.names(MDS_copt$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Coptoformica", cex=80 * meta$het[match(row.names(MDS_copt$points),  meta$CATALOGUENUMBER)])

# With Serviformica cunicularia/rufibarbis
curu.vcf <- vcf[,c(1, which(MDS_all$points[,2] > 10)+1)]
curu.gl <- vcfR2genlight(curu.vcf)
dist_curu <- dist(curu.gl)
MDS_curu <- isoMDS(d = dist_curu, k = 2, maxit = 100)
plot(MDS_curu$points[,2] ~ MDS_curu$points[,1], col=meta$col.morpho[match(row.names(MDS_curu$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Serviformica cunicularia / rufibarbis", cex=80 * meta$het[match(row.names(MDS_curu$points),  meta$CATALOGUENUMBER)])

# With Serviformica lemani/fusca etc.
lem.vcf <- vcf[,c(1, which(MDS_all$points[,1] < -10 & MDS_all$points[,2] < -10)+1)]
lem.gl <- vcfR2genlight(lem.vcf)
dist_lem <- dist(lem.gl)
MDS_lem <- isoMDS(d = dist_lem, k = 2, maxit = 100)
plot(MDS_lem$points[,2] ~ MDS_lem$points[,1], col=meta$col.morpho[match(row.names(MDS_lem$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Serviformica lemani, fusca, etc.", cex=80 * meta$het[match(row.names(MDS_lem$points),  meta$CATALOGUENUMBER)])

# With para/lugu
lug.vcf <- str.vcf[,c(1, which(MDS_str$points[,1] > 0 & MDS_str$points[,2] < 10)+1)]
lug.gl <- vcfR2genlight(lug.vcf)
dist_lug <- dist(lug.gl)
MDS_lug <- isoMDS(d = dist_lug, k = 2, maxit = 100)
plot(MDS_lug$points[,2] ~ MDS_lug$points[,1], col=meta$col.morpho[match(row.names(MDS_lug$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Formica lugubris/paralugubris", cex=80 * meta$het[match(row.names(MDS_lug$points),  meta$CATALOGUENUMBER)])

# With rufa/polyctena
ruf.vcf <- str.vcf[,c(1, which(MDS_str$points[,1] < 0 & MDS_str$points[,2] < 0)+1)]
ruf.gl <- vcfR2genlight(ruf.vcf)
dist_ruf <- dist(ruf.gl)
MDS_ruf <- isoMDS(d = dist_ruf, k = 2, maxit = 100)
plot(MDS_ruf$points[,2] ~ MDS_ruf$points[,1], col=meta$col.morpho[match(row.names(MDS_ruf$points),  meta$CATALOGUENUMBER)], pch=19, las=1, xlab="Dimension 1", ylab="Dimension 2", main="MDS - Formica rufa/polyctena", cex=100 * meta$het[match(row.names(MDS_ruf$points),  meta$CATALOGUENUMBER)])

#################
### admixture ###
#################

### Choose best K
CV <- read.table("admixture/CV_summary_across_rep.txt", h=F)
CV <- CV[order(CV$V1),]
dimnames(CV)[[1]] <- 1:20
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
# we chose K = 13 based on these results
# the replicate with the lowest CV error was replicate 3

k13 <- read.table("admixture/rep2/cleaned.14.Q")
row.names(k13) <- dimnames(gt)[[2]]


par(mfrow=c(10,1), mar=c(2, 2, 0.1, 0.1), oma=c(0.5,0.5,0.1,0.1))
for (i in 1:10) {
  k13 <- read.table(paste("admixture/rep", i, "/cleaned.15.Q", sep=""))
  barplot(t(as.matrix(k13))[,match(tree$tip.label, individuals)], col=c("black", "gold", "tomato", "red", "pink", "green", "darkgreen", "lightblue", "blue", "darkblue", "grey", "brown", "purple", "yellow", "turquoise"), las = 2, cex.names = 0.5, space = 0, border = NA)
}

par(mfrow=c(3,1), mai=c(0.1, 0.6, 0.1, 0.1))
barplot(meta$het[match(tree$tip.label, individuals)], las=1, ylab="Relative heterozygosity", border=NA)
barplot(t(as.matrix(k13))[,match(tree$tip.label, individuals)], las=2, border=NA, space=0, cex.names=0.6, col=c("black", "gold", "tomato", "red", "pink", "green", "darkgreen", "lightblue", "blue", "darkblue", "grey", "brown", "purple", "yellow"))#col=c("#BAB023", "#F9ED32", "#A97C50", "#F15A29", "#EC008C", "#2B3990", "#8DC63F", "#006838", "#1B75BC"))
plot.phylo(tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
legend("bottomleft", legend=legend_Formica, pch=15, col=col_Formica, bty="n")
tiplabels(pch = 15, col=meta$col.morpho[match(tree$tip.label, meta$CATALOGUENUMBER)])

###############################
### admixture sensu stricto ###
###############################

### Choose best K
CV <- read.table("sstr_admixture/CV_summary_across_rep.txt", h=F)
CV <- CV[order(CV$V1),]
dimnames(CV)[[1]] <- 1:10
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
# we chose K = 8 based on these results
# the replicate with the lowest CV error was replicate 3

k <- read.table("sstr_admixture/rep2/cleaned.5.Q")
row.names(k) <- sstr

sstr.tree <- drop.tip(phy = tree, tip = tree$tip.label[tree$tip.label %in% sstr == F], trim.internal = T)

par(mfrow=c(3,1), mai=c(0.1, 0.6, 0.1, 0.1))
barplot(meta$het[match(sstr.tree$tip.label, sstr)], las=1, ylab="Relative heterozygosity", border=NA)
barplot(t(as.matrix(k))[,match(sstr.tree$tip.label, sstr)], las=2, border=NA, space=0, cex.names=0.6, col=c("black", "gold", "tomato", "red", "pink", "green", "darkgreen", "lightblue", "blue", "darkblue", "grey", "brown", "purple", "yellow"))#col=c("#BAB023", "#F9ED32", "#A97C50", "#F15A29", "#EC008C", "#2B3990", "#8DC63F", "#006838", "#1B75BC"))
plot.phylo(sstr.tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
legend("bottomleft", legend=legend_Formica, pch=15, col=col_Formica, bty="n")
tiplabels(pch = 15, col=meta$col.morpho[match(sstr.tree$tip.label, individuals)])


par(mfrow=c(11,1), mar=c(2, 2, 0.1, 0.1), oma=c(0.5,0.5,0.1,0.1))
for (i in 1:10) {
  k <- read.table(paste("sstr_admixture/rep", i, "/cleaned.5.Q", sep=""))
  barplot(t(as.matrix(k))[,match(sstr.tree$tip.label, sstr)], col=c("black", "gold", "tomato", "red", "pink", "green", "darkgreen", "lightblue", "blue", "darkblue", "grey", "brown", "purple", "yellow", "turquoise"), las = 2, cex.names = 0.5, space = 0, border = NA)
}
plot.phylo(sstr.tree, align.tip.label=T, cex=0.4, direction="upwards", root.edge=T)
tiplabels(pch = 15, col=meta$col.morpho[match(sstr.tree$tip.label, individuals)])


order.structure <- function(Q.matrix){
  # number of clusters
  k <- dim(Q.matrix)[2]
  # function to find the maximum value in a vector
  whichismax <- function (x) which(x == max(x, na.rm = T))
  # assign species to each individual by finding which k
  # the individual has most contribution from
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


##################
### genetic ID ###
##################
meta$sp.nucl <- NA
meta$sp.nucl[which(MDS_all$points[,1] > -10 & MDS_all$points[,1] < 0 & MDS_all$points[,2] > -10)] <- "Form_sang"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_str$points[,2] > 10)))] <- "Form_prat"
# F. paralugubris is thought to be a hybrid species between aquilonia and lugubris
# so we consider that the group which clusters with aquilonia is paralugubris.
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_lug$points[,1] > 0.5)))] <- "Form_para"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_lug$points[,1] < 0.5)))] <- "Form_lugu"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_ruf$points[,1] > 1)))] <- "Form_poly"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_ruf$points[,1] < 1)))] <- "Form_rufa"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_lem$points[,2] > -5 & MDS_lem$points[,1] > 0)))] <- "Form_fusc"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_lem$points[,2] > -5 & MDS_lem$points[,1] < 0)))] <- "Form_lema"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_curu$points[,1] < 0)))] <- "Form_cuni"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_curu$points[,1] > 0)))] <- "Form_rufi"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_copt$points[,1] > 0)))] <- "Form_pres"
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% names(which(MDS_copt$points[,1] < 0)))] <- "Form_exse"

# Remove individuals belonging to F. aquilonia and truncorum
meta$sp.nucl[which(meta$CATALOGUENUMBER %in% c("17294", "18702", "Faqui10", "Faqui1", "Faqui5", "Faqui6", "Faqui8"))] <- NA

write.table(meta, "metadata_filled.csv", quote = F, row.names = F, sep=";")

################
### 2nd step ###
################

meta$retained_for_step_2 <- F
meta$retained_for_step_2[which(is.na(meta$sp.nucl) == F)] <- T

Formica.pop <- data.frame(CATALOGUENUMBER=meta$CATALOGUENUMBER[which(meta$retained_for_step_2 == T)], pop = rep("-", sum(meta$retained_for_step_2)))
#Formica.pop <- read.table("2nd_admixture/Formica.pop")
names(Formica.pop) <- c("CATALOGUENUMBER", "pop")
nref = 5
cutoff = 0.9995
pure <- row.names(k13[which(apply(k13, 1, max) > cutoff),])

# for each species, take the nref individuals with lowest missing data among
for (i in (unique(meta$sp.nucl[which(is.na(meta$sp.nucl)==F)]))) {
  Formica.pop$pop[which(Formica.pop$CATALOGUENUMBER %in% meta$CATALOGUENUMBER[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)][order(meta$percent_missing[which(meta$sp.nucl == i & meta$CATALOGUENUMBER %in% pure)])][1:nref])] <- i
}

write.table(Formica.pop, "Formica.pop", row.names = F, col.names = F, quote = F)



##########################
### Introgression rate ###
##########################

### Choose best K again

CV2 <- read.table("2nd_admixture/CV_results.txt")
which(CV2$V4 == min(CV2$V4))
# the replicate with the lowest CV error was replicate 6

k12_2 <- read.table("2nd_admixture/rep11/cleaned.12.Q")
dimnames(k12_2)[[1]] <- Formica.pop$CATALOGUENUMBER


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

k12_renamed <- rename.admixture.columns(k12_2, meta)

barplot(t(as.matrix(order.structure(k12_renamed))), col=c("#00A19A", "#006633", "#A97C50", "#CBBBA0", "#9D1E63", "#D7DF23", "#D45E4D", "#603913", "#29235C", "#1B75BC", "#BE1E2D", "#8DC63F"), las = 2, cex.names = 0.1, space = 0, border = NA)

# select only hybrids
hybrids2nd <- which(apply(k12_renamed, 1, max) < 0.999)

barplot(t(as.matrix(order.structure(k12_renamed[hybrids2nd,]))), col=c("#00A19A", "#006633", "#A97C50", "#CBBBA0", "#9D1E63", "#D7DF23", "#D45E4D", "#603913", "#29235C", "#1B75BC", "#BE1E2D", "#8DC63F"), las = 2, cex.names = 1, space = 0, border = NA)



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

hybridization_rate.matrix <- hybridization.rates(k12_renamed)
write.table(hybridization_rate.matrix, "Formica_hybridization_rates.txt", quote = F)

################
### Figure 1 ###
################

## Step 1: Figure out amount of introgression per species

# create vector for introgression amount
intr_amount <- numeric(dim(k12_renamed)[2])
names(intr_amount) <- dimnames(k12_renamed)[[2]]

# Get maximum contribution of each species
max_admixture <- apply(k12_renamed, 1, max)

# For each species i
for (i in 1:dim(k12_renamed)[2]) {
    # Get which individuals belong to this species
    ismax <- k12_renamed[,i] == max_admixture

    # estimate the amout of introgression in that species as
    # the sum of the values of other columns for individuals of that species
    # divided by the number of individuals of that species
    # i.e. average amount of introgression in the species
    intr_amount[i] <- sum(k12_renamed[ismax,-i]) / (sum(ismax))
}

sort(intr_amount)

## Step 2: plot "own species" in one colour
## and "introgressed" (all species mixed) in another
# 
# # assign a species name to each sample
# whichismax <- function (x) which(x == max(x, na.rm = T))
# species_column <- apply(k12_renamed, 1, whichismax)
# species_ID = colnames(k12_renamed)[species_column]


# create new df with value for "own species" and "introgressed"
for (i in 1:dim(k12_renamed)[[1]]){
  newadmix$introgressed[i] <- sum(k12_renamed[i,]) - max(k12_renamed[i,])
  newadmix$own[i] <- max(k12_renamed[i,])
}


# Compute maximum assignment probability for all individuals
max_Q <- apply(k12_renamed, 1, max)

# Identify the majority cluster name for each individual
major_cluster <- names(intr_amount)[apply(k12_renamed, 1, which.max)]

## Compute sum of probabilities for all non-majority clusters
sum_other <- rowSums(k12_renamed) - max_Q

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
Q.for.figure$pure2[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10,12)])] <- Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10,12)])]
Q.for.figure$pure[which(Q.for.figure$majority_cluster %in% names(sort(intr_amount))[c(2,4,6,8,10,12)])] <- 0


## Step 3: Plot
par(mfrow=c(1,1), mai=c(0.1, 0.5, 0.1, 0.1))
barplot(t(as.matrix(Q.for.figure[,c(1,2,4)])), las = 2, cex.names = 1, space = 0, border = NA, col=c("#A3195B", "grey80", "grey90"))

#########################################
### Genetic distances between species ###
#########################################

# This function returns a symmetrical matrix of N x N dimension, where N is the number of species.
# Each cell contains the average distance (same scale as the tree) between individuals of the two corresponding species.
# tree is an object of class = "phylo".
# meta is a dataframe containing at least two columns: 
#   - CATALOGUENUMBER is the name of each individual, with values matching tree$tip.label
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
write.table(COI_distance.reordered, "Formica_COI_distances.txt", quote = F)

#######################
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

par(mfrow=c(1,1))
plot(hybridization.df$hybridization ~ -hybridization.df$divergence, las=1, xlab="COI divergence", ylab="Admixture proportion inherited from other species", pch=19, col=rgb(0,0,0,0.2))
write.table(hybridization.df, "Formica_hybridization_df.txt", row.names = F, quote = F)



#########################
### Species abundance ###
#########################
table(meta$sp.nucl)
write.table(table(meta$sp.nucl), "spp_abundance.txt", sep = "\t", quote = F)
###


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
write.table(co_occ_BDM, "Formica_co_occurence_BDM.txt", row.names = F, quote = F)


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
write.table(co_occ_plot, "Formica_co_occurence_fine_scale.txt", row.names = F, quote = F)

###################################
### Pi with and without hybrids ###
###################################

# Defining different thresholds for hybrid exclusion.

introgression.thresholds <- c(0.5, 0.9, 0.99, 0.9999)

# Get amount of "pure" genome per individual
maxQ <- apply(k12_2, 1, max, na.rm=T)

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


########################################
### Output data for niche estimation ###
########################################

geo.all <- data.frame(CATALOGUENUMBER=meta$CATALOGUENUMBER, SP=meta$sp.nucl, X=meta$X, Y=meta$Y)
# some samples have "0" as coordinates. Define them as NA
geo.all$X[which(geo.all$X == 0)] <- NA
geo.all$Y[which(geo.all$Y == 0)] <- NA

# retain only pure individuals
# defined as having at least 0.999 assignment to one species in the supervised admixture
pure <- names(maxQ)[which(maxQ > 0.999)]
geo <- subset(geo.all, geo.all$CATALOGUENUMBER %in% pure)

# write to a file
write.table(geo, "Formica_for_niche.csv", sep=";", row.names = F, quote = F)
