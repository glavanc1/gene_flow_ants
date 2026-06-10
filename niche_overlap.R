library(terra)
library(dplyr)
library(ecospat)

#load occ data
occ<-read.csv2("Z:/FBM/DEE/aguisan/opfourmis/D2c/50_introgression/all_sp_occurences.csv",header=FALSE)
names(occ)<-c("site","species","X","Y","source")
table(occ$species, occ$source)
sp.list<-unique(occ$species)

######################## bias correction

#clara
cl8<-rast("Z:/FBM/DEE/aguisan/opfourmis/D2c/40_MSc/results/bias_corrected_predictions/cluster_selection_CLUSTERFILES/CLARA_PCA_clusters_8.tif")
plot(cl8)

occ_cl8<-cbind(occ,extract(cl8,occ[,3:4])[,2])
names(occ_cl8)[6]<-"cl8"
occ_cl8<-na.exclude(occ_cl8) # remove 9 samples falling outside of cl8


# loop through species

occ_all<-data.frame(matrix(nrow=0,ncol=6)) #table where corrected occurrence are stored
names(occ_all)<-c("site","species","X","Y","source","cl8")
pval<-c()
nocc_pub<-c()
nocc_sci<-c()
nocc_all<-c()

for(i in 1:length(sp.list)){
  sp<-sp.list[i]
  occi<-filter(occ_cl8,species==sp)
  tbl<-table(occi$source,occi$cl8)
  
  # Chi2 test
  if(nrow(tbl)==1) pval[i]<-1
  if (nrow(tbl)==2){
    xsq<-chisq.test(tbl)
    pval[i]<-round(xsq$p.value,3)
  }
  
  #store number of occurrence
  nocc_pub[i]<-sum(tbl[which(row.names(tbl)=="public"),])
  nocc_sci[i]<-sum(tbl[which(row.names(tbl)=="scientific"),])
  
  # biais correction
  if(pval[i]<0.05){

    if(sum(tbl["scientific",])<5) {
      nocc_all[i]<-0
      next()
      }
    tbl<-tbl[,which(tbl["scientific",]>0)] #remove cl categories where no scientific occ
    fsci<-tbl[2,]/sum(tbl[2,])
    fpub<-round(sum(tbl[1,])*fsci)
    for (cl in names(fpub)){
      fpubi<-fpub[which(names(fpub)==cl)]

      occipub<-filter(occi,source=="public") %>% filter(cl8==cl)
      min<-min(nrow(occipub),fpubi)
      occipubcor<-occipub[sample(1:nrow(occipub),min),]
      occiscicor<-filter(occ_cl8,species==sp) %>% filter(cl8==cl) %>% filter(source=="scientific")
      #add data to db
      occ_all<-rbind(occ_all,occiscicor,occipubcor)
    }
  }
  if(pval[i]>=0.05){
    occiall<-filter(occ_cl8,species==sp)
    occ_all<-rbind(occ_all,occiall)
  }
  #store number of occurrence
  nocc_all[i]<-nrow(filter(occ_all,species==sp))
}

sp.list_pval<-data.frame(cbind(sp.list,nocc_pub,nocc_sci,nocc_all,pval))
sum(sp.list_pval$pval<0.05)/length(sp.list)


######################## overlap calculation

#load env data
pca_rast<-rast("pca_rast.tif")

scores.clim<-spatSample(pca_rast,10000,na.rm=TRUE)

cb<-combn(sp.list_pval$sp.list[-10],2)
cb<-rbind(cb,rep(NA,ncol(cb)))
row.names(cb)<-c("sp1","sp2","D")

# species loop
for (i in 1:ncol(cb)){
  scores.sp1<-extract(pca_rast,filter(occ_all,species==cb[1,i])[,3:4],ID=FALSE)
  z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1, R = 100)
  scores.sp2<-extract(pca_rast,filter(occ_all,species==cb[2,i])[,3:4],ID=FALSE)
  z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2, R = 100)
  cb[3,i]<-ecospat.niche.overlap(z1,z2,cor=FALSE)$D
  cat(".")
}

m<-matrix(nrow=39,ncol=39)
colnames(m)<-sp.list_pval$sp.list[-10]
rownames(m)<-sp.list_pval$sp.list[-10]

for (i in 1:ncol(cb)){
  m[which(rownames(m)==cb[2,i]),which(colnames(m)==cb[1,i])]<-as.numeric(cb[3,i])
}

corrplot::corrplot(m,is.cor=FALSE)

write.table(m,"overlap.txt",sep="\t")
