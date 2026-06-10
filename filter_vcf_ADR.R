# this script takes two arguments: 
# $1 is a vcf file 
#Ã$2 is a threshold value, which defines the percentile of the expected ADR distribution
# that is used as cutoff for above which genotypes are corrected.
# it takes values between 0 and 1 (default = 0.01, meaning the top 1%)

### Part 1: load and format data

# get arguments from the script
# "trailingonly = T" makes sure only user-specified arguments are used
args <- commandArgs(trailingOnly = T)

print("argument 1:")
print(args[1])
print("argument 2:")
print(args[2])

# Install or load required packages
packages <- c("vcfR")
install.packages(setdiff(packages, rownames(installed.packages())), repos = "http://cran.us.r-project.org")
library(vcfR)

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

# get the name of the vcf from the first given argument
vcf_name <- args[1]

# get the threshold above which to filter putative contaminants
# defined as 1 - second given argument because we want the 1 - xth percentile
# 0.99 by default
threshold = 0.99
threshold = 1 - as.numeric(args[2])

# remove ".vcf" or ".recode.vcf" at the end of the vcf file name
# if there are other "." in the name, only what's in front will be retained
name = unlist(strsplit(x = vcf_name, split = "[.]"))[1]

# load the vcf file in R
vcf <- read.vcfR(vcf_name)

# extract genotypes from the vcf
gt <- extract.gt(vcf, "GT")

# extract depth values from the vcf
dp <- extract.gt(vcf, "DP", as.numeric=T)

print(paste("dimensions of gt:", dim(gt)))

### Part two: Visualize the distribution of Allelic Depth Ratio (ADR) across samples

visualize_ADR_distribution <- function(vcf = vcf) {
  #This function produces a plot of the distribution of Alleic Depth Ratio (ADR)
  # across samples. It takes one argument:
  # "vcf" is a vcf file, read in R with vcfR::read.vcfR()
  
  #get allelic depth for the reference allele at heterozygous sites
  ad.ref <<- extract.gt(vcf, "AD", as.numeric = T)
  ad.ref[which(dp < 11)] <- NA
  ad.ref.het <<- ad.ref
  ad.ref.het[which(gt %in% c("0/0", "1/1"))] <- NA
  
  # and the alternative allele
  ad.alt <<- dp - ad.ref
  ad.alt[which(dp < 11)] <- NA
  ad.alt.het <<- ad.alt
  ad.alt.het[which(gt %in% c("0/0", "1/1"))] <- NA
  
  # computes allelic depth ratio (ADR), defined as
  # depth at the allele supported by most reads divided by total depth at the locus.
  # takes two matrices with allelic depths (integer) as input
  get.ad.ratio <<- function(ad.ref, ad.alt) return(pmax(ad.ref, ad.alt) / mapply(psum, ad.ref, ad.alt, na.rm=T))
  
  # use the function on the whole dataset
  # <<- makes it a global variable (accessible outside the function)
  AD.ratio <<- get.ad.ratio(as.matrix(ad.ref.het), as.matrix(ad.alt.het))
  
  # visualize the distribution of ADR at heterozygous loci per individual
#  png(paste("ADR_distribution_per_sample", name, ".png", sep=""), width = 600, height = 600)
#  par(mfrow=c(1,1), mai=c(1, 0.2, 0.2, 0.2))
#  plot(NULL, xlim=c(0.45, 1), ylim=c(0,12), xlab="Allelic depth ratio (heterozygous alleles only)", ylab="", yaxt="n", las=1)
#  for (i in 1:dim(AD.ratio)[[2]]) {
#    points(density(AD.ratio[,i][which(AD.ratio[,i] != 1)], na.rm=T) , type="l", col="#00000030")
#  }
#  dev.off()
}

print("visualizing AD ratio distribution...")
visualize_ADR_distribution(vcf)

print("step one complete. Dimensions of AD.ratio:")
print(dim(AD.ratio))  
### Part 3: estimate an ADR cutoff above which to correct genotypes

simulate_ADR <- function(gt = gt, dp = dp, AD.ratio = AD.ratio) {
  # This function simulates expected ADR
  # given the observed distribution of depths at heterozygous genotypes
  # assuming no contamination and no biases
  
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
  
  # visualize overall observed distribution of ADR at heterozygous loci
  # and expected distribution based on depth at these loci
  png(paste("ADR_observed_expected", name, ".png", sep=""), width = 600, height = 600)
  
  # Plot a histogram of observed ADR at heterozygous genotypes
  hist(AD.ratio, freq=F, las=1, main="", xlab="")
  
  # probably not the most efficient way but it works.
  # applies sim.ad.ref on observed depth values at heterozygous loci
  sim.ad.ref <- unlist(lapply(X = dp.het, FUN = sim.ad))
  
  # create a dataframe with simulated AD at both alleles
  ad.df <- data.frame(A=sim.ad.ref, B=dp.het - sim.ad.ref)
  
  # compute ADR for the simulated AD values
  sim.ad.ratio <- get.ad.ratio(apply(ad.df, 1, max), apply(ad.df, 1, min))
  
  # plot it as a density
  points(density(sim.ad.ratio, bw = 0.01) , type="l", col="#1D71B8")
  
  # define the cutoff for further filtering
  # as the "thresholdth" percentile of the simulated distribution
  cutoff <<- sort(sim.ad.ratio)[threshold*length(sim.ad.ratio)]
  abline(v=cutoff, col="red", lw=2)
  dev.off()
}

simulate_ADR(gt = gt, dp = dp, AD.ratio = AD.ratio)

### Part 4: filtering the vcf

vcf.filtering <- function(vcf=vcf, cutoff=cutoff, gt=gt, dp=dp, ad.ref, ad.alt) {
  # copying gt
  gt.filtered <- gt
  # copying dp
  dp.filtered <- dp
  dp.filtered[which(is.na(gt.filtered))] <- NA
  # copying ad.ref
  ad.ref.filtered <- ad.ref
  ad.ref.filtered[which(is.na(gt.filtered))] <- NA
  # creating ad.alt to be filtered
  ad.alt.filtered <- dp - ad.ref
  ad.alt.filtered[which(is.na(gt.filtered))] <- NA
  
  # identifying cells of the genotype matrix with ADR > cutoff
  # and assigning them to either "conta.alt" or "conta.ref"
  # in "conta.alt", the alternative allele has the lower coverage and is assumed to be the result of a contamination
  # in "conta.ref", the reference allele has the lower coverage and is assumed to be the result of a contamination
  
  conta.alt <- which(AD.ratio > cutoff & ad.ref.filtered > ad.alt.filtered)
  conta.ref <- which(AD.ratio > cutoff & ad.alt.filtered > ad.ref.filtered)
  
  # correcting genotype values
  gt.filtered[conta.alt] <- "0/0"
  gt.filtered[conta.ref] <- "1/1"
  
  # correcting depth values
  dp.filtered[conta.alt] <- ad.ref.filtered[conta.alt]
  dp.filtered[conta.ref] <- ad.alt.filtered[conta.ref]
  
  # correcting allelic depth values
  ad.alt.filtered[conta.alt] <- 0
  ad.ref.filtered[conta.ref] <- 0
  
  # copying the original vcf
  vcf.filtered <- vcf
  
  # modifying the "FORMAT" field
  vcf.filtered@gt[,1] <- "GT:DP:AD"
  
  # creating the other columns of the vcf@gt by pasting GT, DP and AD
  # and replacing into vcf.filtered to keep all metadata
  vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]] <- paste(gt.filtered, dp.filtered, paste(ad.ref.filtered, ad.alt.filtered, sep=","), sep=":")
  vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]][which(vcf.filtered@gt[,2:dim(vcf.filtered@gt)[2]] == "NA:NA:NA,NA")] <- "./.:.:.,."

  return (vcf.filtered@gt)
}

vcf.filtered <- vcf.filtering(vcf=vcf, cutoff=cutoff, gt=gt, dp=dp, ad.ref, ad.alt)

write.vcf(vcf.filtered, file = paste(name, "_ADR_filtered.vcf", sep=""))
