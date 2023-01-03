#'---
#'title: Validation in GSE141256
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))
suppressMessages(library(binom))
suppressMessages(library(limma))
suppressMessages(library(sva))
suppressMessages(library(pamr))
suppressMessages(library(GEOquery))
suppressMessages(library(GEOmetadb))

suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(rafalib))


options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))
suppressMessages(source(here("scripts","00_heat_scree_plot_generic.R")))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))


gse <- getGEO("GSE141256", GSEMatrix = TRUE)

GSE141256_meta_450K<-pData(gse[[2]]) # 1 is expression 2 is 450K
GSE141256_meta_EPIC<-pData(gse[[3]]) # 3 is EPIC

GSE141256_meta_450K<-GSE141256_meta_450K[,c(1,2,8,10:13,21,24,33)]
GSE141256_meta_EPIC<-GSE141256_meta_EPIC[,c(1,2,8,10:13,21,24,33)]
identical(colnames(GSE141256_meta_EPIC), colnames(GSE141256_meta_450K))


#cd data/published_organoids/GSE141256
#wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141256/suppl/GSE141256_RAW.tar'
#tar -xvf GSE141256_RAW.tar

path<-"data/published_organoids/GSE141256"

GSE141256_meta_450K$Assay.Name<-paste(GSE141256_meta_450K$geo_accession,"_",GSE141256_meta_450K$title, sep="")
GSE141256_meta_EPIC$Assay.Name<-paste(GSE141256_meta_EPIC$geo_accession,"_",GSE141256_meta_EPIC$title, sep="")

GSE141256_meta_450K$array.id.path <- file.path(here(path), GSE141256_meta_450K$Assay.Name)
GSE141256_meta_EPIC$array.id.path <- file.path(here(path), GSE141256_meta_EPIC$Assay.Name)

GSE141256_meta_450K$sentrix_ID<-sapply(1:nrow(GSE141256_meta_450K), function(x){
  strsplit(GSE141256_meta_450K$title[x], "_")[[1]][1]
})

GSE141256_meta_EPIC$sentrix_ID<-sapply(1:nrow(GSE141256_meta_EPIC), function(x){
  strsplit(GSE141256_meta_EPIC$title[x], "_")[[1]][1]
})


#' ### Normalize 450K Arrays
rgset_450k <- read.metharray(GSE141256_meta_450K$array.id.path, verbose = FALSE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina405K <- preprocessFunnorm(rgset_450k, sex=GSE141256_meta_450K$characteristics_ch1.2)
GSE141256_beta_450K<-getBeta(MSet.illumina405K)

#Detection pvalue analysis
avg_detPval <- colMeans(detectionP(rgset_450k))
GSE141256_meta_450K$det_pval<-avg_detPval


#' ### Normalize EPIC Arrays
rgset_EPIC <- read.metharray(GSE141256_meta_EPIC$array.id.path, verbose = FALSE)
class(rgset_EPIC)
reset = updateObject(rgset_EPIC)
class(rgset_EPIC)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illuminaEPIC <- preprocessFunnorm(rgset_EPIC, sex=GSE141256_meta_EPIC$characteristics_ch1.2)
GSE141256_beta_EPIC<-getBeta(MSet.illuminaEPIC)

#Detection pvalue analysis
avg_detPval <- colMeans(detectionP(rgset_EPIC))
GSE141256_meta_EPIC$det_pval<-avg_detPval


#' ### 450K QC and Probe Filtering 

#' ### Probe Filtering 
#' #### Sex Chromosomes 
anno_450k<-anno_450k[match(rownames(GSE141256_beta_450K),anno_450k$IlmnID),]

          # GSE141256_beta_450K<-GSE141256_beta_450K[which(!(anno_450k$CHR%in%c('X','Y'))),] #485512
          # filt_sex<-nrow(GSE141256_beta_450K)
          # print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))

#' #### Cross-hybridizing probes and polymorphic probes. 
#' Some probes have been found to cross-hybridize with other chromosomes (Price et al. 2013 *Epigenetics*).
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
price<-read.table(here("data","GPL16304-47833.txt"), sep='\t', header=T, skip=22)
price<-price[match(rownames(GSE141256_beta_450K),price$ID),]

cross_hyb<-price[which(price$XY_Hits=="XY_YES" | price$Autosomal_Hits=="A_YES"),]
GSE141256_beta_450K<-GSE141256_beta_450K[which(!(rownames(GSE141256_beta_450K)%in%cross_hyb$ID)),]
filt_cross<-nrow(GSE141256_beta_450K)
print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))

#' Polymorphic probes
SnpatCpG<-price[which(price$Target.CpG.SNP!=""),] # 20696
GSE141256_beta_450K<-GSE141256_beta_450K[which(!(rownames(GSE141256_beta_450K)%in%SnpatCpG$ID)),]
filt_poly<-nrow(GSE141256_beta_450K)
print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))


#' #### Probe filtering based on detection pvalue and detection over background (NA)
#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(GSE141256_beta_450K), function(y) length(which(is.na(GSE141256_beta_450K[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(GSE141256_beta_450K)*0.05))
GSE141256_beta_450K<-GSE141256_beta_450K[na_count_probe_good,]
filt_bead<-nrow(GSE141256_beta_450K)
print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))


#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_450k)
failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

GSE141256_beta_450K<-GSE141256_beta_450K[which(!(rownames(GSE141256_beta_450K)%in%bad_det_p)),]
GSE141256_beta_450K<-GSE141256_beta_450K[,which(!(colnames(GSE141256_beta_450K)%in%bad_det_psamp))]

filt_detp<-nrow(GSE141256_beta_450K)
print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))






#' ### EPIC QC and Probe Filtering 
#' ### Probe Filtering 
# SNP probes should already be removed
GSE141256_beta_EPIC <- GSE141256_beta_EPIC[!grepl("rs",rownames(GSE141256_beta_EPIC)), ]
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))

        #' #' #### Sex Chromosomes
        #' anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(GSE141256_beta_EPIC),]
        #' identical(rownames(GSE141256_beta_EPIC),anno_EPIC$IlmnID)
        #' GSE141256_beta_EPIC <- GSE141256_beta_EPIC[!anno_EPIC$CHR%in%c("X", "Y"), ]
        #' filt_sex<-nrow(GSE141256_beta_EPIC)
        #' print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes. 
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#' "43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv
cross_reactive<-read.csv(here("data", "13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[which(!(rownames(GSE141256_beta_EPIC)%in%cross_reactive$PROBE)),]
filt_cross<-nrow(GSE141256_beta_EPIC)
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))


#'For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites." and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes" but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("data", "13059_2016_1066_MOESM4_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(polymorphic$PROBE))," polymorphic probes (genetic variants at targeted CpG sites).", sep=""))

baseext<-read.csv(here("data", "13059_2016_1066_MOESM5_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(baseext$PROBE))," polymorphic probes (single base extension sites for Infinium Type I probes).", sep=""))

GSE141256_beta_EPIC<-GSE141256_beta_EPIC[which(!(rownames(GSE141256_beta_EPIC)%in%c(polymorphic$PROBE, baseext$PROBE))),]
filt_poly<-nrow(GSE141256_beta_EPIC)
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))


#' #### Probe filtering based on detection pvalue and detection over background (NA)
#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(GSE141256_beta_EPIC), function(y) length(which(is.na(GSE141256_beta_EPIC[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(GSE141256_beta_EPIC)*0.05))
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[na_count_probe_good,]
filt_bead<-nrow(GSE141256_beta_EPIC)
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))


#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_EPIC)
detP<-detP[,which(colnames(detP)%in%GSE141256_meta_EPIC$Assay.Name)]
identical(colnames(detP),GSE141256_meta_EPIC$Assay.Name)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

GSE141256_beta_EPIC<-GSE141256_beta_EPIC[which(!(rownames(GSE141256_beta_EPIC)%in%bad_det_p)),]
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[,which(!(colnames(GSE141256_beta_EPIC)%in%bad_det_psamp))]
filt_detp<-nrow(GSE141256_beta_EPIC)

print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))




#' ### Combine the EPIC and 450K data
GSE141256_beta_450K<-GSE141256_beta_450K[which(rownames(GSE141256_beta_450K)%in%rownames(GSE141256_beta_EPIC)),]
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[which(rownames(GSE141256_beta_EPIC)%in%rownames(GSE141256_beta_450K)),]
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[match(rownames(GSE141256_beta_450K),rownames(GSE141256_beta_EPIC)),]
identical(rownames(GSE141256_beta_EPIC),rownames(GSE141256_beta_450K))
GSE141256_beta_combo<-cbind(GSE141256_beta_450K,GSE141256_beta_EPIC)

identical(colnames(GSE141256_meta_EPIC),colnames(GSE141256_meta_450K))
GSE141256_meta_combo<-rbind(GSE141256_meta_450K, GSE141256_meta_EPIC)

identical(colnames(GSE141256_beta_combo), GSE141256_meta_combo$Assay.Name)

print(paste("With combining the EPIC and 450K there are ",ncol(GSE141256_beta_combo)," samples and ",nrow(GSE141256_beta_combo)," CpGs",sep=""))


#' Restructure meta
GSE141256_meta_combo$cell_type<-sapply(1:nrow(GSE141256_meta_combo), function(x){strsplit(GSE141256_meta_combo$characteristics_ch1[x],": ")[[1]][[2]]})
GSE141256_meta_combo$age<-sapply(1:nrow(GSE141256_meta_combo), function(x){strsplit(GSE141256_meta_combo$characteristics_ch1.1[x],": ")[[1]][[2]]})
GSE141256_meta_combo$sex<-sapply(1:nrow(GSE141256_meta_combo), function(x){strsplit(GSE141256_meta_combo$characteristics_ch1.2[x],": ")[[1]][[2]]})
GSE141256_meta_combo$batch<-sapply(1:nrow(GSE141256_meta_combo), function(x){strsplit(GSE141256_meta_combo$characteristics_ch1.3[x],": ")[[1]][[2]]})
GSE141256_meta_combo<-GSE141256_meta_combo[,c(1:3,8,9,11,13:18)]





#'## Combat array EPIC and 450K

#' impute 0 and 1
GSE141256_beta_combo[GSE141256_beta_combo==0]<-0.01
GSE141256_beta_combo[GSE141256_beta_combo==1]<-0.99

#' impute NA
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
GSE141256_beta_combo<-t(imputeMedianv3(GSE141256_beta_combo))

Mval<-function(beta) log2(beta/(1-beta))
edata = apply(GSE141256_beta_combo, 1, Mval) # need mvalues for combat
edata = as.data.frame(edata)
edata = t(edata)

batch = GSE141256_meta_combo$platform_id
combat_GSE141256_mval = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)

#' Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_GSE141256_Beta = apply(combat_GSE141256_mval, 1, betas) # need mvalues for combat
combat_GSE141256_Beta = as.data.frame(combat_GSE141256_Beta)
combat_GSE141256_Beta = t(combat_GSE141256_Beta)
combat_GSE141256_Beta<-as.data.frame(combat_GSE141256_Beta)

combat_GSE141256_Beta<-t(imputeMedianv3(combat_GSE141256_Beta))

save(combat_GSE141256_Beta, GSE141256_meta_combo, file=paste(here("data"),"/GSE141256_beta_organoids_withXY.RData",sep=""))

#load(here("data","GSE141256_beta_organoids_withXY.RData"))

