setwd("C:/Users/glavanc1/Documents/ants_hybridisation/All_genera/")

library(lme4)
#library(car)
#library(vegan)

#######################################
### Reading and formatting the data ###
#######################################

meta.all <- read.csv("../Databases/Metadata_combined_public_scientific.csv", sep=";", h=T)

Camp.hyb <- read.table("Camponotus_hybridization_df.txt", h=T)
Form.hyb <- read.table("Formica_hybridization_df.txt", h=T)
Lasi.hyb <- read.table("Lasius_hybridization_df.txt", h=T)
Myrm.hyb <- read.table("Myrmica_hybridization_df.txt", h=T)
Tapi.hyb <- read.table("Tapinoma_hybridization_df.txt", h=T)
Temn.hyb <- read.table("Temnothorax_hybridization_df.txt", h=T)
Tetr.hyb <- read.table("Tetramorium_hybridization_df.txt", h=T)

Camp.ab <- read.table("Camponotus_spp_abundance.txt", h=T)
Form.ab <- read.table("Formica_spp_abundance.txt", h=T)
Lasi.ab <- read.table("Lasius_spp_abundance.txt", h=T)
Myrm.ab <- read.table("Myrmica_spp_abundance.txt", h=T)
Tapi.ab <- read.table("Tapinoma_spp_abundance.txt", h=T)
Temn.ab <- read.table("Temnothorax_spp_abundance.txt", h=T)
Tetr.ab <- read.table("Tetramorium_spp_abundance.txt", h=T)

Camp.fineco <- read.table("Camponotus_co_occurence_fine_scale.txt", h=T)
Form.fineco <- read.table("Formica_co_occurence_fine_scale.txt", h=T)
Lasi.fineco <- read.table("Lasius_co_occurence_fine_scale.txt", h=T)
Myrm.fineco <- read.table("Myrmica_co_occurence_fine_scale.txt", h=T)
Tapi.fineco <- read.table("Tapinoma_co_occurence_fine_scale.txt", h=T)
Temn.fineco <- read.table("Temnothorax_co_occurence_fine_scale.txt", h=T)
Tetr.fineco <- read.table("Tetramorium_co_occurence_fine_scale.txt", h=T)

Camp.co <- read.table("Camponotus_co_occurence_BDM.txt", h=T)
Form.co <- read.table("Formica_co_occurence_BDM.txt", h=T)
Lasi.co <- read.table("Lasius_co_occurence_BDM.txt", h=T)
Myrm.co <- read.table("Myrmica_co_occurence_BDM.txt", h=T)
Tapi.co <- read.table("Tapinoma_co_occurence_BDM.txt", h=T)
Temn.co <- read.table("Temnothorax_co_occurence_BDM.txt", h=T)
Tetr.co <- read.table("Tetramorium_co_occurence_BDM.txt", h=T)

Form.date <- read.table("Formica_pubished_dates.txt", h=T)
Lasi.date <- read.table("Lasius_published_dates.txt", h=T)
Myrm.date <- read.table("Myrmica_published_dates.txt", h=T)
Temn.date <- read.table("Temnothorax_published_dates.txt", h=T)

pheno <- read.table("Phenology_days.csv", h=T, sep=";")

pi <- read.table("pi_all_genera.txt", h=T)
dxy <- read.table("dxy_all_genera.txt", h=T)

overlap=read.table("../niche_overlap/overlap.txt", h=T)

add.abundance.to.df <- function(hyb, ab, co, fineco, dates = NULL) {
  hyb2 <- hyb
  hyb2$N_sp1 = NA
  hyb2$abundance_sp1 = NA
  hyb2$N_sp2 = NA
  hyb2$abundance_sp2 = NA
  hyb2$co_occurence = NA
  hyb2$fine_scale_co_occurence = NA
  hyb2$published_date = NA
  for (i in 1:dim(hyb2)[1]) {
    hyb2$N_sp1[i] <- ab$Freq[which(ab$Var1 == hyb2$sp1[i])]
    hyb2$abundance_sp1[i] <- ab$Freq[which(ab$Var1 == hyb2$sp1[i])] / sum(ab$Freq)
    hyb2$N_sp2[i] <- ab$Freq[which(ab$Var1 == hyb2$sp2[i])]
    hyb2$abundance_sp2[i] <- ab$Freq[which(ab$Var1 == hyb2$sp2[i])] / sum(ab$Freq)
    hyb2$co_occurence[i] <- co[which(colnames(co) == hyb2$sp1[i]), which(colnames(co) == hyb2$sp2[i])]
    hyb2$fine_scale_co_occurence[i] <- fineco[which(colnames(fineco) == hyb2$sp1[i]), which(colnames(fineco) == hyb2$sp2[i])]
    if(!missing(dates)) {
      hyb2$published_date[i] <- dates[which(colnames(dates) == hyb2$sp1[i]), which(colnames(dates) == hyb2$sp2[i])]
    }
  }
return(hyb2)
}

Camp.full <- add.abundance.to.df(Camp.hyb, Camp.ab, Camp.co, Camp.fineco)
Form.full <- add.abundance.to.df(Form.hyb, Form.ab, Form.co, Form.fineco, Form.date)
Lasi.full <- add.abundance.to.df(Lasi.hyb, Lasi.ab, Lasi.co, Lasi.fineco, Lasi.date)
Myrm.full <- add.abundance.to.df(Myrm.hyb, Myrm.ab, Myrm.co, Myrm.fineco, Myrm.date)
Tapi.full <- add.abundance.to.df(Tapi.hyb, Tapi.ab, Tapi.co, Tapi.fineco)
Temn.full <- add.abundance.to.df(Temn.hyb, Temn.ab, Temn.co, Temn.fineco, Temn.date)
Tetr.full <- add.abundance.to.df(Tetr.hyb, Tetr.ab, Tetr.co, Tetr.fineco)

# home-made mess to clean Form.full
#Form2 <- Form.full
#Form.full <- data.frame(sp1 = Form2$sp1, sp2 = Form2$sp2, hybridization = Form2$hybridization, divergence = Form2$divergence, N_sp1 = Form2$N_sp1, abundance_sp1 = Form2$abundance_sp1, N_sp2 = Form2$N_sp2, abundance_sp2 = Form2$abundance_sp2, co_occurence = Form2$co_occurence, published_date = Form2$boro_date)

# Combine all datasets
all <- rbind(Camp.full, Form.full, Lasi.full, Myrm.full, Tapi.full, Temn.full, Tetr.full)

all <- subset(all, all$sp1 != all$sp2)
# define a "genus" variable (four first characters of sp1)
all$genus <- substr(all$sp1, 1, 4)
all$genus.factor <- as.factor(all$genus)

# Add dxy, da & pi
all$dxy <- NA
all$pi_sp1 <- NA
all$pi_sp2 <- NA
for (sp1 in levels(as.factor(all$sp1))) {
  for (sp2 in levels(as.factor(all$sp2))) {
    line <- which(all$sp1 == sp1 & all$sp2 == sp2)
    # '%in% c(sp1, sp2)' is because the species pair are only in one order in the dxy file
    dxy_here <- sum(dxy$count_diffs[which(dxy$pop1 %in% c(sp1, sp2) & dxy$pop2 %in% c(sp1, sp2) & dxy$chromosome != "FsiM_PB_v5_scf3")], na.rm = T) / sum(dxy$count_comparisons[which(dxy$pop1 %in% c(sp1, sp2) & dxy$pop2 %in% c(sp1, sp2) & dxy$chromosome != "FsiM_PB_v5_scf3")], na.rm = T)
    pi_1 <- sum(pi$count_diffs[which(pi$pop == sp1 & pi$chromosome != "FsiM_PB_v5_scf3")], na.rm = T) / sum(pi$count_comparisons[which(pi$pop == sp1 & pi$chromosome != "FsiM_PB_v5_scf3")], na.rm = T)
    pi_2 <- sum(pi$count_diffs[which(pi$pop == sp2 & pi$chromosome != "FsiM_PB_v5_scf3")], na.rm = T) / sum(pi$count_comparisons[which(pi$pop == sp2 & pi$chromosome != "FsiM_PB_v5_scf3")], na.rm = T)
    all$dxy[line] <- dxy_here
    all$pi_sp1[line] <- pi_1
    all$pi_sp2[line] <- pi_2
  }
}
all$da <- all$dxy - (all$pi_sp1 + all$pi_sp2) / 2


# Define a colour per genus
all$col <- NA
all$col[which(all$genus == "Camp")] <- "#29235C"
all$col[which(all$genus == "Form")] <- "#A3195B"
all$col[which(all$genus == "Lasi")] <- "#00A19A"
all$col[which(all$genus == "Myrm")] <- "#3BBDDC"
all$col[which(all$genus == "Tapi")] <- "#F39200"
all$col[which(all$genus == "Temn")] <- "#95C11F"
all$col[which(all$genus == "Tetr")] <- "#2D2E83"
colours <- c("#29235C", "#A3195B", "#00A19A", "#3BBDDC", "#F39200", "#95C11F", "#2D2E83")
genera <- c("Camponotus", "Formica", "Lasius", "Myrmica", "Tapinoma", "Temnothorax", "Tetramorium")
# define semitransparent colours
all$col50 <- paste(all$col, "77", sep="")
all$col10 <- paste(all$col, "1A", sep="")


# log-transform hybridization values
all$hybridization.log10 <- log10(all$hybridization)
all$hybridization.log10_nozero <- all$hybridization.log10

# Removing zero values (10^-5, or sometimes very close)
all$hybridization.log10_nozero[which(all$hybridization.log10_nozero < -4.9)] <- NA

### Add phenology overlap as Jaccard similarity index of the two species

### this function computes Jaccard similarity index
jaccard <- function(min1, max1, min2, max2) {
  intersect <- min(max1, max2) - max(min1, min2)
  union <- max(max1, max2) - min(min1, min2)
  jaccard <- intersect / union
  return(jaccard)
}

# Ugly function to add phenology overlap information
# df must have columns "sp1" and "sp2"
# pheno must have a column "sp" coded the same way as sp1 and sp2,
# or with "_1" and "_2" at the end if there are two periods
# and a column "day_min" and a "day_max" with days coded as number (January 1st = 1, December 31st = 365)
add.pheno_overlap <- function(df, pheno) {
  pheno_overlap <- numeric(length=length(df$sp1))
  for (i in 1:length(pheno_overlap)) {
    if (length(which(substr(pheno$sp, 1, 9)== df$sp1[i])) > 1) {
      # the two lines of sp1
      line1.1 <- which(substr(pheno$sp,1,9) == df$sp1[i])[1]
      line1.2 <- which(substr(pheno$sp,1,9) == df$sp1[i])[2]
      line2 <- which(pheno$sp == df$sp2[i])
      # intersect = sum of intersect of the two periods of sp1 with sp2
      intersect <- (min(pheno$day_max[line1.1], pheno$day_max[line2]) - max(pheno$day_min[line1.1], pheno$day_min[line2])) + (min(pheno$day_max[line1.2], pheno$day_max[line2]) - max(pheno$day_min[line1.2], pheno$day_min[line2]))
      union <- max(pheno$day_max[line1.1], pheno$day_max[line1.2], pheno$day_max[line2]) - min(pheno$day_min[line1.1], pheno$day_min[line1.2], pheno$day_min[line2])
      pheno_overlap[i] <- intersect / union
    } else if (length(which(substr(pheno$sp, 1, 9)== df$sp2[i])) > 1) {
      line1 <- which(pheno$sp == df$sp1[i])
      line2.1 <- which(substr(pheno$sp,1,9) == df$sp2[i])[1]
      line2.2 <- which(substr(pheno$sp,1,9) == df$sp2[i])[2]
      # intersect = sum of intersect of sp1 with the two periods of sp2
      intersect <- (min(pheno$day_max[line1], pheno$day_max[line2.1]) - max(pheno$day_min[line1], pheno$day_min[line2.1])) + (min(pheno$day_max[line1], pheno$day_max[line2.2]) - max(pheno$day_min[line1], pheno$day_min[line2.2]))
      union <- max(pheno$day_max[line1], pheno$day_max[line2.1], pheno$day_max[line2.2]) - min(pheno$day_min[line1], pheno$day_min[line2.1], pheno$day_min[line2.2])
      pheno_overlap[i] <- intersect / union
    } else {
      line1 <- which(pheno$sp == df$sp1[i])
      line2 <- which(pheno$sp == df$sp2[i])
      pheno_overlap[i] <- jaccard(min1 = pheno$day_min[line1], max1 = pheno$day_max[line1], min2 = pheno$day_min[line2], max2 = pheno$day_max[line2])
    }
  }
  pheno_overlap[which(pheno_overlap <= 0)] <- 0
  return(pheno_overlap)
}

all$pheno_overlap <- add.pheno_overlap(all, pheno)


## Add niche overlap
#copy lower triangle to upper triangle
overlap[upper.tri(overlap)] <- t(overlap)[upper.tri(overlap)]

all$niche_overlap = NA
for (i in 1:length(all$niche_overlap)) {
  sp1 = all$sp1[i]
  sp2 = all$sp2[i]
  if (sum(c(all$sp1[i] != "Form_poly", all$sp2[i] != "Form_poly")) == 2) {
    all$niche_overlap[i] <- overlap[which(row.names(overlap) == sp1), which(row.names(overlap) == sp2)]
  }
}

###########################################
### Distribution of hybridization rates ###
###########################################

library(ridgeline)
par(mfrow=c(1,1))
with(all[which(is.na(all$hybridization.log10)==F),], ridgeline(x=hybridization.log10, y = genus.factor, bw = 0.1,  palette=colours))
with(all[which(is.na(all$hybridization.log10_nozero)==F),], ridgeline(x=hybridization.log10_nozero, y = genus.factor, bw = 0.15,  palette=colours))

# proportion of species pairs with non-null admixture
sum(all$hybridization.log10 > -4.9, na.rm=T) / sum(is.na(all$hybridization.log10) == F)


### Comparison with grey zone
plot(all$da ~ all$dxy, col=all$col, pch=19, las=1)
abline(a=0, b=1)

all$da.log10 <- log10(all$da)
plot(all$hybridization.log10 ~ all$da.log10, pch=19, xlim=c(-5, -1), col=all$col, las=1, xlab="da (log10)", ylab="Introgression proportion (log10)")
abline(v=log10(0.015), lwd=2, col="grey")


plot(all$hybridization.log10 ~ all$da, pch=19, xlim=c(0.00005, 0.028), col=all$col, las=1, xlab="da (log10)", ylab="Introgression proportion (log10)")
abline(v=0.015, lwd=2, col="grey")
abline(v=0.011, lwd=1, col="lightgrey")
abline(v=0.027, lwd=1, col="lightgrey")

# Proportion of individual species with non-null admixture
species <- as.character(levels(as.factor(all$sp1)))
hyb_per_species <- data.frame(sp = species, hyb = numeric(length(species)))
for (i in species) {
  hyb_per_species$hyb[which(hyb_per_species$sp == i)] <- sum(all$hybridization.log10_nozero[which(all$sp1 == i | all$sp1 == i)], na.rm=T)
}
hyb_per_species$hyb[which(hyb_per_species$hyb < 0 )] <- 1

table(hyb_per_species$hyb)
# Proportion of species that engaged in hybridization among genotyped species
table(hyb_per_species$hyb)[2] / sum(table(hyb_per_species$hyb))

meta.all <- read.csv("/Users/glavanc1/Documents/ants_hybridisation/Databases/Metadata_combined_public_scientific.csv", h=T, sep=";")
# Removing clades not identified to species level,
# categories not matching a species
# Tetramorium alpestre and semilaevae (wrongly identified)
# and Formica aquilonia (added as references but not from the area)
meta.all$SPECIESSUMMARY[which(meta.all$SPECIESSUMMARY %in% c("0", "Camp", "Camp_exotique", "Camp_herc/lign", "Camp_pice ?", "Espèce exotique", "Form", "Form Form", "Form Serv", "Form_fusca/lema", "Form_aqui", "Form_poly/rufa", "Form_lugu/para", "Lasi", "Lasi brun", "Lasi Chtono", "Lasi jaune", "Lasi_alie gr", "Lasi_emar/brun", "Lasi_nige/plat", "Lept", "Lept_gred/musc", "Myrm", "non fourmis", "non identifiable", "Plag", "Pone", "Sole", "Tapi", "Tapi_erra/nige_gr", "Tapi_nige_gr", "Temn", "Temn_nyla gr", "Temn_tube gr", "Temn_unif gr", "Tetr", "Tetr_alpe", "Tetr_caes gr", "Tetr_semi", "tube vide"))] <- NA
# adding Tapinoma magnum and Tapinoma darioi (detected gentically)
all_species <- c(as.character(levels(as.factor(meta.all$SPECIESSUMMARY))), "Tapi_magn", "Tapi_dari")
levels(as.factor(meta.all$SPECIESSUMMARY))

# Proportion of species that engaged in hybridization among all species found in the area
table(hyb_per_species$hyb)[2] / length(all_species)


### Figure 1C:
barplot(matrix(c(length(all_species)-table(hyb_per_species$hyb)[2], table(hyb_per_species$hyb)[2])), beside=F, las=1, col=c("grey90", "grey20"))


### Figure 1D: distribution of introgressed amount per species

# Create a new dataset with one entry per species
# and merge introgression from all "donors" for each "receiving"

per_sp <- data.frame(sp = unique(sort(all$sp1)), introgressed = NA)
for (i in per_sp$sp) {
  line = which(per_sp$sp == i)
  which_are_introgressed <- which(all$sp1 == i & all$hybridization > 1.1e-05) # 1e-5 is "zero" for admixture and a few species have very slightly higher values, which I interpret as artifacts.
  per_sp$introgressed[line] <- sum(all$hybridization[which_are_introgressed])
}
per_sp$genus <- as.factor(substr(per_sp$sp, 1, 4))
per_sp$introgressed.log10 <- log10(per_sp$introgressed)

per_sp_nonzero <- per_sp[which(per_sp$introgressed > 0),]

library(ridgeline)
par(mfrow=c(1,1), mai=c(0.5, 1, 0.1, 0.1))
ridgeline(x=per_sp_nonzero$introgressed.log10, y = per_sp_nonzero$genus, bw = 0.15,  palette=colours)

##############################################

plot(het.ind$het ~ het.ind$introgressed.log10, pch=20, col=het.ind$col.transp, las=1)

#################################################################
### Correlation between COI divergence and species divergence ###
#################################################################

plot(all$divergence ~ all$published_date, col=all$col, pch=20)

par(mfrow=c(3,1), mai=c(0.8, 0.8, 0.1, 0.1))
with(all[which(all$genus == "Form"),], plot(divergence ~ published_date, pch=19, col=col50, las=1, xlab="Species divergence [Mya] (Borowiec et al. 2021)", ylab="COI divergence"))
div_lm_Form <- with(all[which(all$genus == "Form"),], lm(divergence ~ published_date))
summary(div_lm_Form)
abline(div_lm_Form, lty=2)

with(all[which(all$genus == "Myrm"),], plot(divergence ~ published_date, pch=19, col=col50, las=1, xlab="Species divergence [Mya] (Jansen et al. 2010)", ylab="COI divergence"))
div_lm_Myrm <- with(all[which(all$genus == "Myrm"),], lm(divergence ~ published_date))
summary(div_lm_Myrm)
abline(div_lm_Myrm, lty=2)

with(all[which(all$genus == "Temn"),], plot(divergence ~ published_date, pch=19, col=col50, las=1, xlab="Species divergence [Mya] (Prebus 2017)", ylab="COI divergence"))
div_lm_Temn <- with(all[which(all$genus == "Temn"),], lm(divergence ~ published_date))
summary(div_lm_Temn)
abline(div_lm_Temn, lty=2)


##################################
### Introgression ~ divergence ###
##################################

# Trying all possible measures of divergence, all give very concordant results
par(mfrow=c(2,2))
# dxy (what was actually used in Figure 3)
plot(all$hybridization.log10 ~ all$dxy, pch=19, col=all$col, las=1, xlab="dxy", ylab="Introgression proportion (log10)")

# da (estimated mostly to compare with Camille Roux's grey "zone")
plot(all$hybridization.log10 ~ all$da, pch=19, col=all$col, las=1, xlab="da", ylab="Introgression proportion (log10)")

# COI divergence
plot(all$hybridization.log10 ~ all$divergence, pch=19, col=all$col, las=1, xlab="COI divergence", ylab="Introgression proportion (log10)")

# Published species divergence time in My (available only for a subset of species)
plot(all$hybridization.log10 ~ all$published_date, pch=19, col=all$col, las=1, xlab="Published species divergence time [My]", ylab="Introgression proportion (log10)")
legend("topright", legend=genera, col=colours, pch=19, bty="n")



### Testing dxy with randomizations


# number of replicates for each test
nrep = 10000

### Only dxy ###
plot(all$hybridization.log10 ~ all$dxy, pch=19, col=all$col, las=1, xlab="dxy", ylab="Introgression proportion (log10)")
LM.div <- lm(hybridization.log10 ~ dxy, data=all)
summary(LM.div)

newdata.dxy <- seq(from=min(all$dxy), to = max(all$dxy), length=length(all$dxy))
pred.dxy <- predict(LM.div, newdata = data.frame(dxy=newdata.dxy), interval = "confidence", type="response")

polygon(c(rev(newdata.dxy), newdata.dxy), c(rev(pred.dxy[ ,3]), pred.dxy[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.dxy[1], x1 = tail(newdata.dxy, 1), y0 = pred.dxy[1,1], y1 = tail(pred.dxy, 1)[,1])

dxy.perm.t <- numeric(0)

for (i in 1:nrep) {
  new.intro <- all$hybridization.log10
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list
    sp.randomized <- sample(sp.ordered)
    for (j in which(all$genus == genus)) {
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      new.intro[j] <- all$hybridization.log10[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  t <- summary(lm(new.intro ~ all$dxy))$coefficients[,3]
  dxy.perm.t <- c(dxy.perm.t, t[2])
}

hist(dxy.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(dxy.perm.t < summary(LM.div)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.div)$coefficients[,3][2], col="red")
print(paste("p-value = ",  sum(dxy.perm.t < summary(LM.div)$coefficients[,3][2])/ nrep))

#####################################
### Introgression ~ niche overlap ###
#####################################

plot(all$hybridization.log10 ~ all$niche_overlap, pch=19, col=all$col, las=1, xlab="Niche overlap", ylab="Introgression proportion (log10)")
LM.niche <- lm(hybridization.log10 ~ niche_overlap, data=all)
summary(LM.niche)

newdata.niche <- seq(from=min(all$niche_overlap, na.rm = T), to = max(all$niche_overlap, na.rm = T), length=length(all$niche_overlap))
pred.niche <- predict(LM.niche, newdata = data.frame(niche_overlap=newdata.niche), interval = "confidence", type="response")

polygon(c(rev(newdata.niche), newdata.niche), c(rev(pred.niche[ ,3]), pred.niche[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.niche[1], x1 = tail(newdata.niche, 1), y0 = pred.niche[1,1], y1 = tail(pred.niche, 1)[,1])


niche.perm.t <- numeric(0)

for (i in 1:nrep) {
  new.intro <- all$hybridization.log10
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list
    sp.randomized <- sample(sp.ordered)
    for (j in which(all$genus == genus)) {
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      new.intro[j] <- all$hybridization.log10[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  t <- summary(lm(new.intro ~ all$niche_overlap))$coefficients[,3]
  niche.perm.t <- c(niche.perm.t, t[2])
}

hist(niche.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(niche.perm.t > summary(LM.niche)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.niche)$coefficients[,3][2], col="red")
print(paste("p-value = ",  sum(niche.perm.t > summary(LM.niche)$coefficients[,3][2])/ nrep))


#######################################################
### Introgression ~ fine-scale spatial co-occurence ###
#######################################################

par(mfrow=c(1,1))
plot(all$hybridization.log10 ~ all$fine_scale_co_occurence, pch=19, col=all$col, las=1, xlab="Fine-scale spatial overlap", ylab="Introgression proportion (log10)")
LM.fineco <- lm(hybridization.log10 ~ fine_scale_co_occurence, data=all)
summary(LM.fineco)

newdata.fineco <- seq(from=min(all$fine_scale_co_occurence, na.rm = T), to = max(all$fine_scale_co_occurence, na.rm = T), length=length(all$fine_scale_co_occurence))
pred.fineco <- predict(LM.fineco, newdata = data.frame(fine_scale_co_occurence=newdata.fineco), interval = "confidence", type="response")

polygon(c(rev(newdata.fineco), newdata.fineco), c(rev(pred.fineco[ ,3]), pred.fineco[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.fineco[1], x1 = tail(newdata.fineco, 1), y0 = pred.fineco[1,1], y1 = tail(pred.fineco, 1)[,1])


fineco.perm.t <- numeric(0)

for (i in 1:nrep) {
  new.intro <- all$hybridization.log10
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list
    sp.randomized <- sample(sp.ordered)
    for (j in which(all$genus == genus)) {
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      new.intro[j] <- all$hybridization.log10[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  t <- summary(lm(new.intro ~ all$fine_scale_co_occurence))$coefficients[,3]
  fineco.perm.t <- c(fineco.perm.t, t[2])
}

hist(fineco.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(fineco.perm.t > summary(LM.fineco)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.fineco)$coefficients[,3][2], col="red")
print(paste("p-value = ",  sum(fineco.perm.t > summary(LM.fineco)$coefficients[,3][2])/ nrep))



#################################
### Introgression ~ phenology ###
#################################

par(mfrow=c(1,1))
plot(all$hybridization.log10 ~ all$pheno_overlap, pch=19, col=all$col, las=1, xlab="Mating phenology overlap", ylab="Introgression proportion (log10)")
LM.pheno <- lm(hybridization.log10 ~ pheno_overlap, data=all)
summary(LM.pheno)

newdata.pheno <- seq(from=min(all$pheno_overlap), to = max(all$pheno_overlap), length=length(all$pheno_overlap))
pred.pheno <- predict(LM.pheno, newdata = data.frame(pheno_overlap=newdata.pheno), interval = "confidence", type="response")

polygon(c(rev(newdata.pheno), newdata.pheno), c(rev(pred.pheno[ ,3]), pred.pheno[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.pheno[1], x1 = tail(newdata.pheno, 1), y0 = pred.pheno[1,1], y1 = tail(pred.pheno, 1)[,1])


pheno.perm.t <- numeric(0)

for (i in 1:nrep) {
  new.intro <- all$hybridization.log10
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list
    sp.randomized <- sample(sp.ordered)
    for (j in which(all$genus == genus)) {
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      new.intro[j] <- all$hybridization.log10[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  t <- summary(lm(new.intro ~ all$pheno_overlap))$coefficients[,3]
  pheno.perm.t <- c(pheno.perm.t, t[2])
}

hist(pheno.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(pheno.perm.t > summary(LM.pheno)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.pheno)$coefficients[,3][2], col="red")
print(paste("p-value = ",  sum(pheno.perm.t > summary(LM.pheno)$coefficients[,3][2])/ nrep))


##################
#### Figure 3 ####
##################
par(mfrow=c(2,4))

hist(pheno.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(pheno.perm.t > summary(LM.pheno)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.pheno)$coefficients[,3][2], col="red")

hist(niche.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(niche.perm.t > summary(LM.niche)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.niche)$coefficients[,3][2], col="red")

hist(fineco.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(fineco.perm.t > summary(LM.fineco)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.fineco)$coefficients[,3][2], col="red")

hist(dxy.perm.t, xlab="", ylab="", las=1, breaks=20, border=F, xlim=c(-10, 10), main=paste("p-value = ",  sum(dxy.perm.t < summary(LM.div)$coefficients[,3][2])/ 10000))
abline(v = summary(LM.div)$coefficients[,3][2], col="red")

plot(all$hybridization.log10 ~ all$pheno_overlap, pch=19, col=all$col, las=1, xlab="Mating phenology overlap", ylab="Introgression proportion (log10)")
polygon(c(rev(newdata.pheno), newdata.pheno), c(rev(pred.pheno[ ,3]), pred.pheno[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.pheno[1], x1 = tail(newdata.pheno, 1), y0 = pred.pheno[1,1], y1 = tail(pred.pheno, 1)[,1])

plot(all$hybridization.log10 ~ all$niche_overlap, pch=19, col=all$col, las=1, xlab="Ecological niche overlap", ylab="Introgression proportion (log10)")
polygon(c(rev(newdata.niche), newdata.niche), c(rev(pred.niche[ ,3]), pred.niche[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.niche[1], x1 = tail(newdata.niche, 1), y0 = pred.niche[1,1], y1 = tail(pred.niche, 1)[,1])

plot(all$hybridization.log10 ~ all$fine_scale_co_occurence, pch=19, col=all$col, las=1, xlab="Fine-scale spatial overlap", ylab="Introgression proportion (log10)")
polygon(c(rev(newdata.fineco), newdata.fineco), c(rev(pred.fineco[ ,3]), pred.fineco[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.fineco[1], x1 = tail(newdata.fineco, 1), y0 = pred.fineco[1,1], y1 = tail(pred.fineco, 1)[,1])

plot(all$hybridization.log10 ~ all$dxy, pch=19, col=all$col, las=1, xlab="dxy", ylab="Introgression proportion (log10)")
polygon(c(rev(newdata.dxy), newdata.dxy), c(rev(pred.dxy[ ,3]), pred.dxy[ ,2]), col = 'grey80', border = NA)
segments(x0 = newdata.dxy[1], x1 = tail(newdata.dxy, 1), y0 = pred.dxy[1,1], y1 = tail(pred.dxy, 1)[,1])
