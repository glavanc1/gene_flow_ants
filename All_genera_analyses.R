setwd("C:/Users/glavanc1/Documents/ants_hybridisation/All_genera/")

library(lme4)
#library(car)
library(glmmTMB)
library(DHARMa)
library(ridgeline)

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
all$col[which(all$genus == "Camp")] <- "#49196B"
all$col[which(all$genus == "Form")] <- "#A3195B"
all$col[which(all$genus == "Lasi")] <- "#00A19A"
all$col[which(all$genus == "Myrm")] <- "#3BBDDC"
all$col[which(all$genus == "Tapi")] <- "#F39200"
all$col[which(all$genus == "Temn")] <- "#95C11F"
all$col[which(all$genus == "Tetr")] <- "#2D2E83"
colours <- c("#49196B", "#A3195B", "#00A19A", "#3BBDDC", "#F39200", "#95C11F", "#2D2E83")
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
      intersect <- max(0, min(pheno$day_max[line1.1], pheno$day_max[line2]) - max(pheno$day_min[line1.1], pheno$day_min[line2])) + max(0, min(pheno$day_max[line1.2], pheno$day_max[line2]) - max(pheno$day_min[line1.2], pheno$day_min[line2]))
      union <- (pheno$day_max[line1.1] - pheno$day_min[line1.1]) + (pheno$day_max[line1.2] - pheno$day_min[line1.2]) + (pheno$day_max[line2] - pheno$day_min[line2]) - intersect
      pheno_overlap[i] <- intersect / union
    } else if (length(which(substr(pheno$sp, 1, 9)== df$sp2[i])) > 1) {
      line1 <- which(pheno$sp == df$sp1[i])
      line2.1 <- which(substr(pheno$sp,1,9) == df$sp2[i])[1]
      line2.2 <- which(substr(pheno$sp,1,9) == df$sp2[i])[2]
      # intersect = sum of intersect of sp1 with the two periods of sp2
      intersect <- max(0, min(pheno$day_max[line1], pheno$day_max[line2.1]) - max(pheno$day_min[line1], pheno$day_min[line2.1])) + max(0, min(pheno$day_max[line1], pheno$day_max[line2.2]) - max(pheno$day_min[line1], pheno$day_min[line2.2]))
      union <- (pheno$day_max[line1] - pheno$day_min[line1]) + (pheno$day_max[line2.1] - pheno$day_min[line2.1]) + (pheno$day_max[line2.2] - pheno$day_min[line2.2]) - intersect      
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

# Add daily phenology overlap
# Convert times into minutes since midnight
pheno$min.minutes <- 60* as.numeric(substr(pheno$hour_min, 1,2)) + as.numeric(substr(pheno$hour_min, 4,5))
pheno$max.minutes <- 60* as.numeric(substr(pheno$hour_max, 1,2)) + as.numeric(substr(pheno$hour_max, 4,5))

## Add niche overlap
#copy lower triangle to upper triangle
overlap[upper.tri(overlap)] <- t(overlap)[upper.tri(overlap)]

all$overlap.hour <- numeric(length=length(all$sp1))
for (i in 1:length(all$overlap.hour)) {
  line1 <- which(substr(pheno$sp, 1, 9) == all$sp1[i])[1]
  line2 <- which(substr(pheno$sp, 1, 9) == all$sp2[i])[1]
  all$overlap.hour[i] <- jaccard(min1 = pheno$min.minutes[line1], max1 = pheno$max.minutes[line1], min2 = pheno$min.minutes[line2], max2 = pheno$max.minutes[line2])
}

all$overlap.hour[which(all$overlap.hour <= 0)] <- 0
all$phenology_overal.total <- all$pheno_overlap * all$overlap.hour

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

################################
### Mapping rate per species ###
################################

Camp.mapping <- read.table("Campo_mapping_summary.txt", h=T)
Camp.mapping$col = "#49196B"
Form.mapping <- read.table("Form_mapping_summary.txt", h=T)
Form.mapping$col = "#A3195B"
Lasi.mapping <- read.table("Lasi_mapping_summary.txt", h=T)
Lasi.mapping$col = "#00A19A"
Myrm.mapping <- read.table("Myrmica_mapping_summary.txt", h=T)
Myrm.mapping$col = "#3BBDDC"
Tapi.mapping <- read.table("Tapinoma_mapping_summary.txt", h=T)
Tapi.mapping$col = "#F39200"
Temn.mapping <- read.table("Temnothorax_mapping_summary.txt", h=T)
Temn.mapping$col = "#95C11F"
Tetr.mapping <- read.table("Tetramorium_mapping_summary.txt", h=T)
Tetr.mapping$col = "#2D2E83"


map.summary <- rbind(Camp.mapping, Form.mapping, Lasi.mapping, Myrm.mapping, Tapi.mapping, Temn.mapping, Tetr.mapping)

plot(map.summary$mapping.mean ~ map.summary$dxy,
     pch = 20,
     col = map.summary$col,
     las = 1,
     cex = 1.5,
     xlab = "dxy",
     ylab = "Percentage of reads mapping",
     ylim=c(84, 105))

for (i in 1:length(map.summary$sp)) {
  segments(x0 = map.summary$dxy[i],
           x1 = map.summary$dxy[i],
           y0 = map.summary$mapping.mean[i] - map.summary$mapping.SD[i],
           y1 = map.summary$mapping.mean[i] + map.summary$mapping.SD[i],
           col = map.summary$col[i])
}

###########################################
### Distribution of introgression rates ###
###########################################

# proportion of species pairs with non-null admixture
sum(all$hybridization.log10 > -4.9, na.rm=T) / sum(is.na(all$hybridization.log10) == F)

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

# Proportion of species that engaged in hybridization among all species found in the area
meta.all <- read.csv("/Users/glavanc1/Documents/ants_hybridisation/Databases/Metadata_combined_public_scientific.csv", h=T, sep=";")
# Removing clades not identified to species level,
# categories not matching a species
# Tetramorium alpestre and semilaevae (wrongly identified)
# and Formica aquilonia (added as references but not from the area)
meta.all$SPECIESSUMMARY[which(meta.all$SPECIESSUMMARY %in% c("0", "Camp", "Camp_exotique", "Camp_herc/lign", "Camp_pice ?", 
                                                             "Espèce exotique", "Form", "Form Form", "Form Serv", "Form_fusca/lema",
                                                             "Form_aqui", "Form_poly/rufa", "Form_lugu/para", "Lasi", "Lasi brun",
                                                             "Lasi Chtono", "Lasi jaune", "Lasi_alie gr", "Lasi_emar/brun", "Lasi_nige/plat",
                                                             "Lept", "Lept_gred/musc", "Myrm", "non fourmis", "non identifiable", "Plag", "Pone",
                                                             "Sole", "Tapi", "Tapi_erra/nige_gr", "Tapi_nige_gr", "Temn", "Temn_nyla gr",
                                                             "Temn_tube gr", "Temn_unif gr", "Tetr", "Tetr_alpe", "Tetr_caes gr", "Tetr_semi", "tube vide"))] <- NA
# adding Tapinoma magnum and Tapinoma darioi (detected gentically)
all_species <- c(as.character(levels(as.factor(meta.all$SPECIESSUMMARY))), "Tapi_magn", "Tapi_dari")
levels(as.factor(meta.all$SPECIESSUMMARY))
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

par(mfrow=c(1,1), mai=c(0.5, 1, 0.1, 0.1))
ridgeline(x=per_sp_nonzero$introgressed.log10, y = per_sp_nonzero$genus, bw = 0.15,  palette = colours[c(2:4, 6, 7)])

########################################
### Correlates of introgression rate ###
########################################

### transforming the 10^-5 value returned by admixture for absence of introgression into 0
all$hybridization2 <- all$hybridization
all$hybridization2[which(all$hybridization == 1e-05)] <- 0

############
### d_XY ###
############

# dxy needs scaling for the model to converge:
all$dxy_scaled <- as.numeric(scale(all$dxy))

# fitting zero-inflated model with genus as random factor
(zi.dxy <- glmmTMB(hybridization2 ~ dxy_scaled + (1|genus.factor),
                   ziformula = ~ dxy_scaled + (1|genus.factor), 
                   family = beta_family(link = "logit"), 
                   data=all))

summary(zi.dxy)

sim_res <- simulateResiduals(zi.dxy)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)
plotResiduals(sim_res, form = all$dxy)
plotResiduals(sim_res, form = all$genus.factor)

dxy.zi.z <- numeric(0)
dxy.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.dxy <- glmmTMB(perm.hybridization ~ dxy_scaled + (1|genus.factor),
                     ziformula = ~ dxy_scaled + (1|genus.factor),
                     family = beta_family(link = "logit"), 
                     data=all)
  
  # Extract z-value for the Zero-Inflated model
  dxy.zi.z <- c(dxy.zi.z, summary(perm.zi.dxy)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  dxy.cond.z <- c(dxy.cond.z, summary(perm.zi.dxy)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(dxy.zi.z, "dxy.zi.z.txt")
write(dxy.cond.z, "dxy.cond.z.txt")

dxy.zi.z <- read.table("dxy.zi.z.txt")$V1
dxy.cond.z <- read.table("dxy.cond.z.txt")$V1

hist(dxy.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(dxy.zi.z > summary(zi.dxy)$coefficients$zi[2,3], na.rm = T)/ (nrep - sum(is.na(dxy.zi.z)))))
abline(v = summary(zi.dxy)$coefficients$zi[2,3], col="red")

hist(dxy.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(dxy.cond.z < summary(zi.dxy)$coefficients$cond[2,3], na.rm = T)/ (nrep - sum(is.na(dxy.cond.z)))))
abline(v = summary(zi.dxy)$coefficients$cond[2,3], col="red")

########################
### Yearly phenology ###
########################

(zi.pheno.day <- glmmTMB(hybridization2 ~ pheno_overlap + (1|genus.factor),
                         ziformula = ~ pheno_overlap + (1|genus.factor), 
                         family = beta_family(link = "logit"), 
                         data=all))

summary(zi.pheno.day)


sim_res <- simulateResiduals(zi.pheno.day)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)
plotResiduals(sim_res, form = all$dxy)
plotResiduals(sim_res, form = all$genus.factor)

pheno.day.zi.z <- numeric(0)
pheno.day.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.pheno.day <- glmmTMB(perm.hybridization ~ pheno_overlap + (1|genus.factor),
                               ziformula = ~ pheno_overlap + (1|genus.factor), 
                               family = beta_family(link = "logit"), 
                               data=all)
  
  # Extract z-value for the Zero-Inflated model
  pheno.day.zi.z <- c(pheno.day.zi.z, summary(perm.zi.pheno.day)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  pheno.day.cond.z <- c(pheno.day.cond.z, summary(perm.zi.pheno.day)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(pheno.day.zi.z, "pheno.day.zi.z.txt")
write(pheno.day.cond.z, "pheno.day.cond.z.txt")

hist(pheno.day.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.day.zi.z < summary(zi.pheno.day)$coefficients$zi[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.day.zi.z)))))
abline(v = summary(zi.pheno.day)$coefficients$zi[2,3], col="red")

hist(pheno.day.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.day.cond.z > summary(zi.pheno.day)$coefficients$cond[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.day.cond.z)))))
abline(v = summary(zi.pheno.day)$coefficients$cond[2,3], col="red")


##################################
#### Ecological niche overlap ####
##################################

(zi.niche <- glmmTMB(hybridization2 ~ niche_overlap  + (1|genus.factor),
                     ziformula = ~ niche_overlap + (1|genus.factor), 
                     family = beta_family(link = "logit"),
                     data=all))

summary(zi.niche)


sim_res <- simulateResiduals(zi.niche)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)
plotResiduals(sim_res, form = all$zi.niche)
plotResiduals(sim_res, form = all$genus.factor)

niche.zi.z <- numeric(0)
niche.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.niche <- glmmTMB(perm.hybridization ~ niche_overlap + (1|genus.factor),
                           ziformula = ~ niche_overlap + (1|genus.factor), 
                           family = beta_family(link = "logit"), 
                           data=all)
  
  # Extract z-value for the Zero-Inflated model
  niche.zi.z <- c(niche.zi.z, summary(perm.zi.niche)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  niche.cond.z <- c(niche.cond.z, summary(perm.zi.niche)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(niche.zi.z, "niche.zi.z.txt")
write(niche.cond.z, "niche.cond.z.txt")
niche.zi.z <- read.table("niche.zi.z.txt")$V1
niche.cond.z <- read.table("niche.cond.z.txt")$V1

hist(niche.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(niche.zi.z < summary(zi.niche)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(niche.zi.z)))))
abline(v = summary(zi.niche)$coefficients$zi[2,3], col="red")

hist(niche.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(niche.cond.z > summary(zi.niche)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(niche.cond.z)))))
abline(v = summary(zi.niche)$coefficients$cond[2,3], col="red")


####################################
#### Fine-scale spatial overlap ####
####################################

(zi.fine <- glmmTMB(hybridization2 ~ fine_scale_co_occurence  + (1|genus.factor),
                     ziformula = ~ fine_scale_co_occurence + (1|genus.factor), 
                     family = beta_family(link = "logit"),
                     data=all))

summary(zi.fine)


sim_res <- simulateResiduals(zi.fine)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)
plotResiduals(sim_res, form = all$zi.fine)
plotResiduals(sim_res, form = all$genus.factor)

fine.zi.z <- numeric(0)
fine.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.fine <- glmmTMB(perm.hybridization ~ fine_scale_co_occurence + (1|genus.factor),
                           ziformula = ~ fine_scale_co_occurence + (1|genus.factor), 
                           family = beta_family(link = "logit"), 
                           data=all)
  
  # Extract z-value for the Zero-Inflated model
  fine.zi.z <- c(fine.zi.z, summary(perm.zi.fine)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  fine.cond.z <- c(fine.cond.z, summary(perm.zi.fine)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(fine.zi.z, "fine.zi.z.txt")
write(fine.cond.z, "fine.cond.z.txt")
fine.zi.z <- read.table("fine.zi.z.txt")$V1
fine.cond.z <- read.table("fine.cond.z.txt")$V1

hist(fine.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(fine.zi.z < summary(zi.fine)$coefficients$zi[2,3], na.rm = T)/ (nrep - sum(is.na(fine.zi.z)))))
abline(v = summary(zi.fine)$coefficients$zi[2,3], col="red")

hist(fine.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(fine.cond.z > summary(zi.fine)$coefficients$cond[2,3], na.rm = T)/ (nrep - sum(is.na(fine.cond.z)))))
abline(v = summary(zi.fine)$coefficients$cond[2,3], col="red")

################
### Figure 3 ###
################
# as a function (written by claude)
all$hyb01 <- as.numeric(all$hybridization2 > 0)
nz <- all$hybridization2 > 0

band <- function(x, lo, hi)
  polygon(c(x, rev(x)), c(lo, rev(hi)), col = "grey95", border = NA)

plot_pair <- function(mod, var, xlab, back = identity, xlim = NULL) {
  # prediction grid on the model's scale
  nd <- data.frame(seq(min(all[[var]], na.rm = TRUE), max(all[[var]], na.rm = TRUE),
                       length.out = 200))
  names(nd) <- var
  nd$genus.factor <- NA
  
  pz <- predict(mod, newdata = nd, type = "zlink", re.form = NA, se.fit = TRUE)
  pc <- predict(mod, newdata = nd, type = "link",  re.form = NA, se.fit = TRUE)
  
  xg <- back(nd[[var]])    # grid on plotting scale
  xo <- back(all[[var]])   # data on plotting scale
  if (is.null(xlim)) xlim <- range(xo, na.rm = TRUE)
  
  # Panel A: probability of hybridization
  plot(xo, all$hyb01, type = "n", xlim = xlim, ylim = c(0, 1),
       xlab = xlab, ylab = "P(hybridization)")
  band(xg, 1 - plogis(pz$fit + 1.96 * pz$se.fit), 1 - plogis(pz$fit - 1.96 * pz$se.fit))
  points(jitter(xo), all$hyb01, pch = "|", col = all$col)
  lines(xg, 1 - plogis(pz$fit), lwd = 2, col="grey60")
  
  # Panel B: non-zero values (log y)
  lo <- plogis(pc$fit - 1.96 * pc$se.fit)
  hi <- plogis(pc$fit + 1.96 * pc$se.fit)
  plot(xo[nz], all$hybridization2[nz], type = "n", log = "y", xlim = xlim,
       ylim = range(c(all$hybridization2[nz], lo, hi), na.rm = TRUE),
       xlab = xlab, ylab = "Hybridization (non-zero)")
  band(xg, lo, hi)
  points(xo[nz], all$hybridization2[nz], pch = 16, col = all$col[nz], cex=1.5)
  lines(xg, plogis(pc$fit), lwd = 2, col="grey60")
}

mods <- list(
  list(zi.pheno.day, "pheno_overlap",           "Phenology overlap",        xlim = c(0, 1)),
  list(zi.niche,     "niche_overlap",           "Niche overlap",            xlim = c(0, 1)),
  list(zi.fine,      "fine_scale_co_occurence", "Fine-scale co-occurrence", xlim = c(0, 1)),
  list(zi.dxy,       "dxy_scaled",              "dxy",
       function(z) z * sd(all$dxy) + mean(all$dxy)), xlim = c(0, 0.25)
)

par(mfcol = c(2, 4), mar = c(4, 4, 1, 1))
for (m in mods) do.call(plot_pair, m)
par(mfcol = c(1, 1))

############################
### Estimate effect size ###
############################

# Effect sizes for a zero-inflated beta model
#   mod   : fitted glmmTMB model
#   var   : name of the predictor exactly as used in the model formula
#   lo,hi : two predictor values to compare (on the scale the model uses)
#   label : text shown in the output table
# ------------------------------------------------------------------
eff_size <- function(mod, var, lo, hi, label) {
  
  # Small data frame with the two predictor values to compare.
  # The column must have the same name as the predictor in the model.
  nd <- data.frame(c(lo, hi)); names(nd) <- var
  # genus.factor = NA + re.form = NA (below) -> predictions for an
  # "average genus", i.e. ignoring the genus random effect
  nd$genus.factor <- NA
  
  # Model predictions on the internal (logit) scale, with standard errors:
  #   type = "zlink" -> zero-inflation part (log-odds of a ZERO)
  #   type = "link"  -> conditional part (logit of the mean non-zero value)
  pz <- predict(mod, newdata = nd, type = "zlink", re.form = NA, se.fit = TRUE)
  pc <- predict(mod, newdata = nd, type = "link",  re.form = NA, se.fit = TRUE)
  
  # Size of the step between the two predictor values
  d  <- hi - lo
  
  # Slope of the predictor (b) and its standard error (s) in each part:
  #   fixef()  gives the estimates
  #   vcov()   gives their variances; sqrt of the diagonal = standard errors
  bz <- fixef(mod)$zi[var];   sz <- sqrt(diag(vcov(mod)$zi))[var]
  bc <- fixef(mod)$cond[var]; sc <- sqrt(diag(vcov(mod)$cond))[var]
  
  data.frame(
    model        = label,
    
    # --- Zero-inflation part: does hybridization happen at all? ---
    # plogis() converts log-odds to a probability (here: probability of zero);
    # 1 - that = probability of hybridizing
    P_hyb_lo     = 1 - plogis(pz$fit[1]),
    P_hyb_hi     = 1 - plogis(pz$fit[2]),
    
    # --- Conditional part: how much hybridization, when it happens? ---
    # Predicted mean of the non-zero values at lo and hi
    mean_nz_lo   = plogis(pc$fit[1]),
    mean_nz_hi   = plogis(pc$fit[2]),

    row.names = NULL
  )
}


# dxy: compare minimum vs maximum observed values.
# Scale the minimum and maximum observed values since dxy is scaled in the model
dxy_lo <- (min(all$dxy, na.rm = TRUE) - mean(all$dxy, na.rm = TRUE)) / sd(all$dxy, na.rm = TRUE)
dxy_hi <- (max(all$dxy, na.rm = TRUE) - mean(all$dxy, na.rm = TRUE)) / sd(all$dxy, na.rm = TRUE)

# Run for the four models and stack results into one table.
# Overlap predictors: compare no overlap (0) vs full overlap (1).
res <- rbind(
  eff_size(zi.dxy,       "dxy_scaled",              dxy_lo, dxy_hi, "dxy"),
  eff_size(zi.pheno.day, "pheno_overlap",           min(all$pheno_overlap, na.rm = T), max(all$pheno_overlap, na.rm = T),           "Phenology overlap"),
  eff_size(zi.niche,     "niche_overlap",           min(all$niche_overlap, na.rm = T), max(all$niche_overlap, na.rm = T),           "Niche overlap"),
  eff_size(zi.fine,      "fine_scale_co_occurence", min(all$fine_scale_co_occurence, na.rm = T), max(all$fine_scale_co_occurence, na.rm = T),           "Fine-scale co-occurrence")
)

# Table 1
print(res, digits = 3)

####################################################################
### For supplementary figure: with hourly and combined phenology ###
####################################################################

#### hourly phenology ####


(zi.pheno.hour <- glmmTMB(hybridization2 ~ overlap.hour + (1|genus.factor),
                         ziformula = ~ overlap.hour + (1|genus.factor), 
                         family = beta_family(link = "logit"), 
                         data=all))

summary(zi.pheno.hour)


sim_res <- simulateResiduals(zi.pheno.hour)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)

pheno.hour.zi.z <- numeric(0)
pheno.hour.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.pheno.hour <- glmmTMB(perm.hybridization ~ overlap.hour + (1|genus.factor),
                               ziformula = ~ overlap.hour + (1|genus.factor), 
                               family = beta_family(link = "logit"), 
                               data=all)
  
  # Extract z-value for the Zero-Inflated model
  pheno.hour.zi.z <- c(pheno.hour.zi.z, summary(perm.zi.pheno.hour)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  pheno.hour.cond.z <- c(pheno.hour.cond.z, summary(perm.zi.pheno.hour)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(pheno.hour.zi.z, "pheno.hour.zi.z.txt")
write(pheno.hour.cond.z, "pheno.hour.cond.z.txt")

hist(pheno.hour.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.hour.zi.z < summary(zi.pheno.hour)$coefficients$zi[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.hour.zi.z)))))
abline(v = summary(zi.pheno.hour)$coefficients$zi[2,3], col="red")

hist(pheno.hour.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.hour.cond.z > summary(zi.pheno.hour)$coefficients$cond[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.hour.cond.z)))))
abline(v = summary(zi.pheno.hour)$coefficients$cond[2,3], col="red")


#### Combined phenology ####
# Defined as the product of hourly * daily

(zi.pheno.tot <- glmmTMB(hybridization2 ~ phenology_overal.total + (1|genus.factor),
                          ziformula = ~ phenology_overal.total + (1|genus.factor), 
                          family = beta_family(link = "logit"), 
                          data=all))

summary(zi.pheno.tot)


sim_res <- simulateResiduals(zi.pheno.tot)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)

pheno.tot.zi.z <- numeric(0)
pheno.tot.cond.z <- numeric(0)

nrep = 10000

start.time <- Sys.time()
for (i in 1:nrep) {
  # Make a new column in the data frame for the permutated introgression values
  all$perm.hybridization <- all$hybridization2
  # Loop across genera
  for (genus in levels(all$genus.factor)) {
    # get the list of species from this genus
    sp.ordered <- (levels(unique(factor(all$sp1[which(all$genus == genus)]))))
    # randomize the species list of this genus
    sp.randomized <- sample(sp.ordered)
    # Within species of this genus
    for (j in which(all$genus == genus)) {
      # permute sp 1
      new.sp1 <- sp.randomized[which(sp.ordered==all$sp1[j])]
      # permute sp 2 the same way
      new.sp2 <- sp.randomized[which(sp.ordered==all$sp2[j])]
      # permute the introgression values accordingly
      all$perm.hybridization[j] <- all$hybridization2[which(all$sp1 == new.sp1 & all$sp2 == new.sp2)]
    }
  }
  # Run the model with the new permuted introgression values
  perm.zi.pheno.tot <- glmmTMB(perm.hybridization ~ phenology_overal.total + (1|genus.factor),
                                ziformula = ~ phenology_overal.total + (1|genus.factor), 
                                family = beta_family(link = "logit"), 
                                data=all)
  
  # Extract z-value for the Zero-Inflated model
  pheno.tot.zi.z <- c(pheno.tot.zi.z, summary(perm.zi.pheno.tot)$coefficients$zi[2,3])
  # Extract z-value for the Conditional model
  pheno.tot.cond.z <- c(pheno.tot.cond.z, summary(perm.zi.pheno.tot)$coefficients$cond[2,3])
}

end.time <- Sys.time()
end.time - start.time

write(pheno.tot.zi.z, "pheno.tot.zi.z.txt")
write(pheno.tot.cond.z, "pheno.tot.cond.z.txt")

hist(pheno.tot.zi.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.tot.zi.z < summary(zi.pheno.tot)$coefficients$zi[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.tot.zi.z)))))
abline(v = summary(zi.pheno.tot)$coefficients$zi[2,3], col="red")

hist(pheno.tot.cond.z, xlab="", ylab="", las=1, border=F, main=paste("p-value = ",  sum(pheno.tot.cond.z > summary(zi.pheno.tot)$coefficients$cond[2,3], na.rm = T)/ (nrep - sum(is.na(pheno.tot.cond.z)))))
abline(v = summary(zi.pheno.tot)$coefficients$cond[2,3], col="red")



plot_pair <- function(mod, var, xlab, back = identity, xlim = NULL) {
  # prediction grid on the model's scale
  nd <- data.frame(seq(min(all[[var]], na.rm = TRUE), max(all[[var]], na.rm = TRUE),
                       length.out = 200))
  names(nd) <- var
  nd$genus.factor <- NA
  
  pz <- predict(mod, newdata = nd, type = "zlink", re.form = NA, se.fit = TRUE)
  pc <- predict(mod, newdata = nd, type = "link",  re.form = NA, se.fit = TRUE)
  
  xg <- back(nd[[var]])    # grid on plotting scale
  xo <- back(all[[var]])   # data on plotting scale
  if (is.null(xlim)) xlim <- range(xo, na.rm = TRUE)
  
  # Panel A: probability of hybridization
  plot(xo, all$hyb01, type = "n", xlim = xlim, ylim = c(0, 1),
       xlab = xlab, ylab = "P(hybridization)")
  band(xg, 1 - plogis(pz$fit + 1.96 * pz$se.fit), 1 - plogis(pz$fit - 1.96 * pz$se.fit))
  points(jitter(xo, amount = 0.1), all$hyb01, pch = "|", col = all$col)
  lines(xg, 1 - plogis(pz$fit), lwd = 2, col="grey60")
  
  # Panel B: non-zero values (log y)
  lo <- plogis(pc$fit - 1.96 * pc$se.fit)
  hi <- plogis(pc$fit + 1.96 * pc$se.fit)
  plot(xo[nz], all$hybridization2[nz], type = "n", log = "y", xlim = xlim,
       ylim = range(c(all$hybridization2[nz], lo, hi), na.rm = TRUE),
       xlab = xlab, ylab = "Hybridization (non-zero)")
  band(xg, lo, hi)
  points(xo[nz], all$hybridization2[nz], pch = 16, col = all$col[nz], cex=1.5)
  lines(xg, plogis(pc$fit), lwd = 2, col="grey60")
}

mods <- list(
  list(zi.pheno.hour, "overlap.hour", "Phenology overlap (daily)", xlim = c(0, 1)),
  list(zi.pheno.tot, "phenology_overal.total", "Phenology overlap (combined)", xlim = c(0, 1))
)
par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
for (m in mods) do.call(plot_pair, m)
par(mfcol = c(1, 1))


#############################################
### d_a: comparisons with the "grey zone" ###
#############################################

# dxy needs scaling for the model to converge:
all$da_scaled <- as.numeric(scale(all$da))

# fitting zero-inflated model with genus as random factor
(zi.da <- glmmTMB(hybridization2 ~ da_scaled + (1|genus.factor),
                   ziformula = ~ da_scaled + (1|genus.factor), 
                   family = beta_family(link = "logit"), 
                   data=all))

summary(zi.da)

sim_res <- simulateResiduals(zi.da)

plot(sim_res)                      # QQ-plot + residual-vs-predicted
testZeroInflation(sim_res)         # is there STILL excess zero problem?
testDispersion(sim_res)            # over/under-dispersion
testOutliers(sim_res)
plotResiduals(sim_res, form = all$da)
plotResiduals(sim_res, form = all$genus.factor)

mods <- list(
  list(zi.da,       "da_scaled",              "da",
       function(z) z * sd(all$da) + mean(all$da), xlim = c(0, 0.027))
)


### Figure S25
par(mfcol = c(2, 1), mar = c(4, 4, 1, 1))
for (m in mods) do.call(plot_pair, m)
par(mfcol = c(1, 1))

abline(v=0.015, lwd=2, col="grey")
abline(v=0.011, lwd=1, col="lightgrey")
abline(v=0.027, lwd=1, col="lightgrey")


### Figure S27
### Distribution of permutation z-values

par(mfcol = c(4, 3), mar = c(2, 0.5, 3, 0.5))

hist(pheno.day.zi.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("Mating phenology\noverlap (seasonal)\np-value = ",  round(sum(pheno.day.zi.z < summary(zi.pheno.day)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(pheno.day.zi.z))), 4)))
abline(v = summary(zi.pheno.day)$coefficients$zi[2,3], col="red")

hist(pheno.day.cond.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("p-value = ",  round(sum(pheno.day.cond.z > summary(zi.pheno.day)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(pheno.day.cond.z))), 4)))
abline(v = summary(zi.pheno.day)$coefficients$cond[2,3], col="red")

hist(niche.zi.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("Ecological niche overlap\np-value = ",  round(sum(niche.zi.z < summary(zi.niche)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(niche.zi.z))), 4)))
abline(v = summary(zi.niche)$coefficients$zi[2,3], col="red")

hist(niche.cond.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("p-value = ",  round(sum(niche.cond.z > summary(zi.niche)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(niche.cond.z))), 4)))
abline(v = summary(zi.niche)$coefficients$cond[2,3], col="red")


hist(fine.zi.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("Fine-scale spatial overlap\np-value = ",  round(sum(fine.zi.z < summary(zi.fine)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(fine.zi.z))), 4)))
abline(v = summary(zi.fine)$coefficients$zi[2,3], col="red")

hist(fine.cond.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("p-value = ",  round(sum(fine.cond.z > summary(zi.fine)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(fine.cond.z))), 4)))
abline(v = summary(zi.fine)$coefficients$cond[2,3], col="red")


hist(dxy.zi.z, xlab="", ylab="", xlim=c(-2.5, 6), las=1, border=F, yaxt="n", main=paste("Genetic divergence\np-value = ",  round(sum(dxy.zi.z > summary(zi.dxy)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(dxy.zi.z))), 4)))
abline(v = summary(zi.dxy)$coefficients$zi[2,3], col="red")

hist(dxy.cond.z, xlab="", ylab="", las=1, xlim=c(-3, 0.5), border=F, yaxt="n", main=paste("p-value = ",  round(sum(dxy.cond.z < summary(zi.dxy)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(dxy.cond.z))), 4)))
abline(v = summary(zi.dxy)$coefficients$cond[2,3], col="red")


hist(pheno.hour.zi.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("Mating phenology\noverlap (daily)\np-value = ",  round(sum(pheno.hour.zi.z < summary(zi.pheno.hour)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(pheno.hour.zi.z))), 4)))
abline(v = summary(zi.pheno.hour)$coefficients$zi[2,3], col="red")

hist(pheno.hour.cond.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("p-value = ",  round(sum(pheno.hour.cond.z > summary(zi.pheno.hour)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(pheno.hour.cond.z))), 4)))
abline(v = summary(zi.pheno.hour)$coefficients$cond[2,3], col="red")


hist(pheno.tot.zi.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("Mating phenology\noverlap (combined)\np-value = ",  round(sum(pheno.tot.zi.z < summary(zi.pheno.tot)$coefficients$zi[2,3], na.rm = T)/ (length(!is.na(pheno.tot.zi.z))), 4)))
abline(v = summary(zi.pheno.tot)$coefficients$zi[2,3], col="red")

hist(pheno.tot.cond.z, xlab="", ylab="", las=1, border=F, yaxt="n", main=paste("p-value = ",  round(sum(pheno.tot.cond.z > summary(zi.pheno.tot)$coefficients$cond[2,3], na.rm = T)/ (length(!is.na(pheno.tot.cond.z))), 4)))
abline(v = summary(zi.pheno.tot)$coefficients$cond[2,3], col="red")

