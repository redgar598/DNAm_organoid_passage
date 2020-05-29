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


grid.arrange(ggplot(GSE141256_meta_450K)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
               geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
                          position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
               xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE),
             
             ggplot(GSE141256_meta_EPIC)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
               geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
                          position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
               xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE))


ggsave(here("figs","GSE141256_detection_pvalue_organoids.pdf"), 
	grid.arrange(ggplot(GSE141256_meta_450K)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
	  geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
	             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
	  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE),

	 ggplot(GSE141256_meta_EPIC)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
	  geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
	             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
	  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)),
	 width=6, height=5)


ggsave(here("figs/jpeg","GSE141256_detection_pvalue_organoids.jpeg"),
	grid.arrange(ggplot(GSE141256_meta_450K)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
	  geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
	             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
	  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE),

	 ggplot(GSE141256_meta_EPIC)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
	  geom_point(aes(as.factor(sentrix_ID), det_pval, group=geo_accession, fill=as.factor(sentrix_ID)), shape=21, color="black",
	             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
	  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)),
	   width=6, height=5)





#' Beta distribution before and after normalization
# 450K
GSE141256_beta_450K_raw<-getBeta(rgset_450k)

Beta_melted<- melt(GSE141256_beta_450K)
Beta_melted_raw<- melt(GSE141256_beta_450K_raw)

Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE141256_meta_450K, by.x="ID", by.y="Assay.Name")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,GSE141256_meta_450K, by.x="ID", by.y="Assay.Name")
beta_dis_450k<-ggplot(Beta_Plot, aes(Beta, group=as.character(geo_accession), color=as.character(description)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")
beta_dis_450k_raw<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(geo_accession), color=as.character(description)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")



# EPIC
GSE141256_beta_EPIC_raw<-getBeta(rgset_EPIC)

Beta_melted<- melt(GSE141256_beta_EPIC)
Beta_melted_raw<- melt(GSE141256_beta_EPIC_raw)

Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE141256_meta_EPIC, by.x="ID", by.y="Assay.Name")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,GSE141256_meta_EPIC, by.x="ID", by.y="Assay.Name")
beta_dis_EPIC<-ggplot(Beta_Plot, aes(Beta, group=as.character(geo_accession), color=as.character(description)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")
beta_dis_EPIC_raw<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(geo_accession), color=as.character(description)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")

grid.arrange(beta_dis_450k_raw,beta_dis_450k,beta_dis_EPIC_raw,beta_dis_EPIC)

ggsave(here("figs","GSE141256_beta_distribution.pdf"),grid.arrange(beta_dis_450k_raw,beta_dis_450k,beta_dis_EPIC_raw,beta_dis_EPIC),  w=10, h=5)
ggsave(here("figs/jpeg","GSE141256_beta_distribution.jpeg"),grid.arrange(beta_dis_450k_raw,beta_dis_450k,beta_dis_EPIC_raw,beta_dis_EPIC), w=10, h=5)



#' ### 450K QC and Probe Filtering 

#' #### Confirm ID with SNPs and cluster by DNAm
# Clustering By sampling site
# remove rows with NAs
GSE141256_450K_betas_cluster<-GSE141256_beta_450K[complete.cases(GSE141256_beta_450K),]

d <- dist(t(GSE141256_450K_betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(GSE141256_meta_450K$description), cex=1.5)

pdf(here('figs','GSE141256_cluster_whole450K_organoid.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(GSE141256_meta_450K$description), cex=1.5)
dev.off()

#' Genotyping Probes
SNPs <- getSnpBeta(rgset_450k)
SNPs<-SNPs[complete.cases(SNPs),]

SNPs<-SNPs[,which(colnames(SNPs)%in%GSE141256_meta_450K$Assay.Name)]
identical(colnames(SNPs),GSE141256_meta_450K$Assay.Name)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(as.character(GSE141256_meta_450K$description)), cex=1.5)

pdf(here('figs','GSE141256_cluster_snps_450K.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(as.character(GSE141256_meta_450K$description)), cex=1.5)
dev.off()



#' #### by sex too
#' #### 450K annotation from illumina
# https://emea.support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
anno_450k<-read.csv(here("data","HumanMethylation450_15017482_v1-2.csv"), skip=7)
anno_450k<-anno_450k[match(rownames(GSE141256_beta_450K),anno_450k$IlmnID),]

identical(rownames(GSE141256_beta_450K),anno_450k$IlmnID)

GSE141256_450K_sex<-GSE141256_beta_450K[which(anno_450k$CHR%in%c('X','Y')),]

d <- dist(t(GSE141256_450K_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(GSE141256_meta_450K$characteristics_ch1.2), cex=1.5)

pdf(here('figs','GSE141256_cluster_sex_450K.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_450K$geo_accession, lab.col=as.fumeric(GSE141256_meta_450K$characteristics_ch1.2), cex=1.5)
dev.off()




#' ### Probe Filtering 
#' #### Sex Chromosomes 
anno_450k<-anno_450k[match(rownames(GSE141256_beta_450K),anno_450k$IlmnID),]

GSE141256_beta_450K<-GSE141256_beta_450K[which(!(anno_450k$CHR%in%c('X','Y'))),] #485512
filt_sex<-nrow(GSE141256_beta_450K)
print(paste("Samples available: ",ncol(GSE141256_beta_450K),"Probes available: ",nrow(GSE141256_beta_450K),sep=""))


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

#' #### Confirm ID with SNPs and cluster by DNAm
# Clustering By sample site
# remove rows with NAs
GSE141256_EPIC_betas_cluster<-GSE141256_beta_EPIC[complete.cases(GSE141256_beta_EPIC),]

d <- dist(t(GSE141256_EPIC_betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(GSE141256_meta_EPIC$description), cex=1.5)


pdf(here('figs','GSE141256_cluster_wholeEPIC_organoid.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(GSE141256_meta_EPIC$description), cex=1.5)
dev.off()

#Genotyping Probes
SNPs <- getSnpBeta(rgset_EPIC)
SNPs<-SNPs[complete.cases(SNPs),]

SNPs<-SNPs[,which(colnames(SNPs)%in%GSE141256_meta_EPIC$Assay.Name)]
identical(colnames(SNPs),GSE141256_meta_EPIC$Assay.Name)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(as.character(GSE141256_meta_EPIC$description)), cex=1.5)

pdf(here('figs','GSE141256_cluster_snps_EPIC.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(as.character(GSE141256_meta_EPIC$description)), cex=1.5)
dev.off()



#' #### by sex too
#' Using the cg ID to chromosome annotation from illumina 
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)
GSE141256_beta_EPIC<-GSE141256_beta_EPIC[which(rownames(GSE141256_beta_EPIC)%in%anno_EPIC$IlmnID),]

anno_EPIC<-anno_EPIC[match(rownames(GSE141256_beta_EPIC),anno_EPIC$IlmnID),]
identical(rownames(GSE141256_beta_EPIC),anno_EPIC$IlmnID)

GSE141256_EPIC_sex<-GSE141256_beta_EPIC[which(anno_EPIC$CHR%in%c('X','Y')),]

d <- dist(t(GSE141256_EPIC_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(GSE141256_meta_EPIC$characteristics_ch1.2), cex=1.5)

pdf(here('figs','GSE141256_cluster_sex_EPIC.pdf'), width=30)
myplclust(hc, labels=GSE141256_meta_EPIC$geo_accession, lab.col=as.fumeric(GSE141256_meta_EPIC$characteristics_ch1.2), cex=1.5)
dev.off()


#' ### Probe Filtering 
# SNP probes should already be removed
GSE141256_beta_EPIC <- GSE141256_beta_EPIC[!grepl("rs",rownames(GSE141256_beta_EPIC)), ]
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))

#' #### Sex Chromosomes
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(GSE141256_beta_EPIC),]
identical(rownames(GSE141256_beta_EPIC),anno_EPIC$IlmnID)
GSE141256_beta_EPIC <- GSE141256_beta_EPIC[!anno_EPIC$CHR%in%c("X", "Y"), ]
filt_sex<-nrow(GSE141256_beta_EPIC)
print(paste("Samples available: ",ncol(GSE141256_beta_EPIC),"\nProbes available: ",nrow(GSE141256_beta_EPIC),sep=""))


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



#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(GSE141256_beta_combo))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

GSE141256_meta_combo$age<-as.numeric(GSE141256_meta_combo$age)

meta_categorical <- GSE141256_meta_combo[, c(3,5,7,9,11,12)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(GSE141256_meta_combo[, c(8, 10)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Sample Site", "Array Type","Sentrix ID", "Cell Type","Sex","Batch")
colnames(meta_continuous) <- c("det_pval", "age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

ggsave(here("figs/GSE141256_heat_scree_before_combat.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)
ggsave(here("figs/jpeg","GSE141256_heat_scree_before_combat.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)




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

save(combat_GSE141256_Beta, GSE141256_meta_combo, file=paste(here("data"),"/GSE141256_beta_organoids.RData",sep=""))

#load(here("data","GSE141256_beta_organoids.RData"))






# PCA on batch corrected data
pca_res <- prcomp(t(combat_GSE141256_Beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

meta_categorical <- GSE141256_meta_combo[, c(3,5,7,9,11,12)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(GSE141256_meta_combo[, c(8, 10)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Sample Site", "Array Type","Sentrix ID", "Cell Type","Sex","Batch")
colnames(meta_continuous) <- c("det_pval", "age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

ggsave(here("figs/GSE141256_heat_scree_after_combat.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)
ggsave(here("figs/jpeg","GSE141256_heat_scree_after_combat.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)






## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, GSE141256_meta_combo, by="Assay.Name")

pc2<-ggplot(Loadings_meta, aes(PC1, PC2, fill=source_name_ch1))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (27%)")+ylab("PC2 (26%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#1a9850", "cornflowerblue","#e1e1ff","#1347a4","#000192"))

pc1<-ggplot(Loadings_meta, aes(PC1, PC2, fill=cell_type))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (27%)")+ylab("PC2 (26%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#35978f","#8c510a","#c7eae5","#dfc27d"))

grid.arrange(pc1, pc2)




#' ## Overall Variance Across most Variable CpGs with Passage
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

## only in spheroids
GSE141256_meta_combo_spheroids<-GSE141256_meta_combo[grep("spheroid", GSE141256_meta_combo$cell_type),]
combat_GSE141256_Beta_spheroids<-combat_GSE141256_Beta[,which(colnames(combat_GSE141256_Beta)%in%GSE141256_meta_combo_spheroids$Assay.Name)]

GSE141256_mval_combo = apply(combat_GSE141256_Beta_spheroids, 1, Mval) # need mvalues for combat
GSE141256_mval_combo = as.data.frame(GSE141256_mval_combo)
GSE141256_mval_combo = t(GSE141256_mval_combo)

ref_range_dnam<-sapply(1:nrow(GSE141256_mval_combo), function(x) Variation(GSE141256_mval_combo[x,]))
GSE141256_beta_VeryVariable<-combat_GSE141256_Beta_spheroids[which(ref_range_dnam>=2.75),]

print(paste("There are ",nrow(GSE141256_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))


#' ## Additonal Sample Information 
#' Shared by Leanne Jones
DNAm_GEOcodes<-read.csv(here("data/published_organoids/GSE141256","Lewis.etal.DNAmSamples_Organoids_GEOcodes-1.csv"))

DNAm_GEOcodes$GSM<-gsub(" ","",as.character(DNAm_GEOcodes$GEO.ID))
DNAm_GEOcodes$GSM<-sapply(1:nrow(DNAm_GEOcodes), function(x){paste(strsplit(DNAm_GEOcodes$GSM[x],"")[[1]][1:10], collapse="")})
intersect(DNAm_GEOcodes$GSM, GSE141256_meta_combo_spheroids$geo_accession)

GSE141256_meta_combo_spheroids<-merge(GSE141256_meta_combo_spheroids, DNAm_GEOcodes[,c(1,7,10)], by.x="geo_accession", by.y="GSM")
GSE141256_meta_combo_spheroids$passage_factor<-as.factor(as.character(GSE141256_meta_combo_spheroids$Passage.Number))
levels(GSE141256_meta_combo_spheroids$passage_factor)<-c(11,2,3,3,4,5,7)
GSE141256_meta_combo_spheroids$passage_numeric<-as.numeric(as.character(GSE141256_meta_combo_spheroids$passage_factor))
GSE141256_meta_combo_spheroids$passage.or.numeric.factor <- factor(GSE141256_meta_combo_spheroids$passage_factor, levels = c(11,7,5,4,3,2))


#' ### PCA on just organoids
pca_res <- prcomp(t(combat_GSE141256_Beta_spheroids))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, GSE141256_meta_combo_spheroids, by="Assay.Name")

pc1<-ggplot(Loadings_meta, aes(PC1, PC2, fill=source_name_ch1))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab(paste("PC1 (",round(Importance[1]*100,0),"%)", sep=""))+ylab(paste("PC2 (",round(Importance[2]*100,0),"%)", sep=""))+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#1a9850", "cornflowerblue","#e1e1ff","#1347a4","#000192"))

pc2<-ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.numeric.factor)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab(paste("PC2 (",round(Importance[2]*100,0),"%)", sep=""))+ylab(paste("PC3 (",round(Importance[3]*100,0),"%)", sep=""))+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))+
  scale_fill_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")


grid.arrange(pc1, pc2)

ggsave(here("figs/GSE141256_organoid_PCA.pdf"),grid.arrange(pc1, pc2),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","GSE141256_organoid_PCA.jpeg"), grid.arrange(pc1, pc2),width = 7.5, height = 5)


# beta plot variable CpGs
Beta_melted<- melt(GSE141256_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE141256_meta_combo_spheroids, by.x="ID", by.y="Assay.Name")
Beta_Plot$passage.or.numeric.factor <- factor(Beta_Plot$passage_factor, levels = c(11,7,5,4,3,2))

ggplot(Beta_Plot, aes(Beta,  color=passage.or.numeric.factor))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=rev(c("#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1","#5E4FA2")), name="Passage\nNumber")

ggsave(here("figs/GSE141256_variable_CpGs.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","GSE141256_variable_CpGs.jpeg"), width = 7.5, height = 5)


#' ### Beta distributions paired by individual
GSE141256_meta_combo_spheroids_paired<-GSE141256_meta_combo_spheroids[grep("F|H",GSE141256_meta_combo_spheroids$Sample.ID),]
GSE141256_meta_combo_spheroids_paired$individual<-as.factor(as.character(GSE141256_meta_combo_spheroids_paired$Sample.ID))
levels(GSE141256_meta_combo_spheroids_paired$individual)<-c("F","F","H","H","J","K")

GSE141256.organoid_paired<-do.call(rbind,lapply(1:length(unique(GSE141256_meta_combo_spheroids_paired$individual)), function(x){
  sample<-as.character(unique(GSE141256_meta_combo_spheroids_paired$individual)[x])
  samp<-GSE141256_meta_combo_spheroids_paired[GSE141256_meta_combo_spheroids_paired$individual==sample,]
  samp<-samp[order(samp$passage_numeric),]
  samp$hilo<-as.factor(as.character(samp$passage_numeric))
  
  if(length(levels(samp$hilo))==2){levels(samp$hilo)<-c("lower","higher")}else{
    if(length(levels(samp$hilo))==3){levels(samp$hilo)<-c("lower","higher","highest")}else{
      if(length(levels(samp$hilo))==4){levels(samp$hilo)<-c("lowest","lower","higher","highest")}else{samp$hilo<-NA}
    }
  }
  samp
}))

GSE141256_beta_VeryVariable_paird<-GSE141256_beta_VeryVariable[,which(colnames(GSE141256_beta_VeryVariable)%in%GSE141256.organoid_paired$Assay.Name)]

Beta_melted<- melt(GSE141256_beta_VeryVariable_paird)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,GSE141256.organoid_paired, by.x="ID", by.y="Assay.Name")
Beta_Plot$individual<-as.character(Beta_Plot$individual)
GSE141256.organoid_paired$individual<-as.character(GSE141256.organoid_paired$individual)

labels<-as.data.frame(tapply(GSE141256.organoid_paired$passage_numeric, GSE141256.organoid_paired$individual, function(x) paste(x, collapse=", ")))
colnames(labels)<-"passge"
labels$individual<-rownames(labels)

ggplot()+
  geom_density(aes(Beta,color=hilo, group=ID),Beta_Plot, size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values = c ("#9ecae1", "#225ea8", "#081d58"), name="Relative\nPassage\nLevel within\nPatient")+facet_wrap(~individual, nrow=1)+
  geom_text(aes(0.75, 2.75, label=passge), data=labels, color="grey20")+th+theme(strip.text = element_text(size = 10),
                                                                                 axis.text=element_text(size=4),
                                                                                 panel.spacing = unit(0.7, "lines"))+th+
  scale_x_continuous(breaks = c(0,0.5,1))

ggsave(here("figs","GSE141256_paired_beta.pdf"),width = 5, height = 2.2)
ggsave(here("figs/jpeg","GSE141256_paired_beta.jpeg"),width = 5, height = 2.2)





#' # Thresholding Trimodality
GSE141256_meta_combo_spheroids$thresholded_prior_ratio<-sapply(1:nrow(GSE141256_meta_combo_spheroids), function(x){
  print(x)
  converted<-as.numeric(round(GSE141256_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  passage_threshold_params(converted, counts, res)
})

GSE141256_meta_combo_spheroids$thresholded_ratio_max<-F
GSE141256_meta_combo_spheroids$thresholded_ratio_max[which(GSE141256_meta_combo_spheroids$thresholded_prior_ratio>1)]<-T

percent_passing<-round((tapply(GSE141256_meta_combo_spheroids$thresholded_ratio_max, GSE141256_meta_combo_spheroids$passage_numeric, sum)/tapply(GSE141256_meta_combo_spheroids$geo_accession, GSE141256_meta_combo_spheroids$passage_numeric, length))*100,2)
passed_num<-tapply(GSE141256_meta_combo_spheroids$thresholded_ratio_max, GSE141256_meta_combo_spheroids$passage_numeric, sum)
org_numer<-tapply(GSE141256_meta_combo_spheroids$geo_accession, GSE141256_meta_combo_spheroids$passage_numeric, length)

df<-data.frame(passage=names(percent_passing), passing=percent_passing, pro_passing=percent_passing/100, count=org_numer, passed_num=passed_num)
df<-cbind(df,(binom.confint(df$passed_num, df$count, method="exact", conf.level=0.95)))
df$upper<-df$upper*100
df$lower<-df$lower*100

ggplot(df, aes(as.numeric(as.character(passage)), passing))+
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey70", width=.25)+
  geom_line(color="grey20")+geom_point(size=1.25,shape=21,color="black",aes(fill=passage))+xlab("Passage")+
  ylab("Samples with Trimodal\nDistribution (%)")+theme_bw()+theme(axis.title = element_text(size=10))+
  scale_fill_manual(values=c("#5E4FA2", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4","#4CA5B1"),name="Passage\nNumber", guide=F)+
  scale_x_continuous(breaks=c(2,3,4,5,7,11))

ggsave(here("figs","GSE141256_Mixture_model_ratio_threshold_maximize.pdf"), width=3, height=2)



#' ### Plot all samples
plts_paired<-lapply(1:nrow(GSE141256_meta_combo_spheroids), function(x){#1:nrow(epic.organoid)
  print(x)
  passage<-paste("passage: ",GSE141256_meta_combo_spheroids$passage_numeric[x],"\nIndividual: ", GSE141256_meta_combo_spheroids$Sample.ID[x],"\nRatio I/H: ",round(GSE141256_meta_combo_spheroids$thresholded_prior_ratio[x],2), sep="")
  converted<-as.numeric(round(GSE141256_beta_VeryVariable[,x]*1000,0))
  counts<-rep(1000, length(converted))
  res = em(converted, counts, .41, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
  draw_fit_params_gg(converted, counts, res,passage)
})

plts_paired_order<-plts_paired[order(GSE141256_meta_combo_spheroids$passage_numeric)]

pdf(here("figs","GSE141256_organoids_thresholding_all_samples.pdf"))
plts_paired_order
dev.off()





#' ## Differential methylation with passage 
mod<-model.matrix(~ 0 + passage_numeric, data=GSE141256_meta_combo_spheroids)
fit <- lmFit(combat_GSE141256_Beta_spheroids, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
beta<-combat_GSE141256_Beta_spheroids

passage_db<-sapply(1:nrow(beta), function(x){
  sampleinfo_cpg<-GSE141256_meta_combo_spheroids
  sampleinfo_cpg$beta<-as.numeric(beta[x,])
  
  fit<-lm(beta ~ passage_numeric, data=sampleinfo_cpg)
  pval<-summary(fit)$coef["passage_numeric","Pr(>|t|)"]
  slope<-fit$coefficients[2]
  
  (min(GSE141256_meta_combo_spheroids$passage_numeric)*slope) - (max(GSE141256_meta_combo_spheroids$passage_numeric)*slope)})

passage_GSE141256<-data.frame(p.value=ebfit$p.value[,"passage_numeric"], CpG=rownames(beta), db=passage_db)

# Adjust P values
passage_GSE141256$p_adjusted<-p.adjust(passage_GSE141256$p.value, method="BH")

diff_CpG_dbGSE141256<-passage_GSE141256[which(passage_GSE141256$p_adjusted<0.05 & abs(passage_GSE141256$db)>0.15),] #25086
diff_CpG_db_hypoGSE141256<-diff_CpG_dbGSE141256$CpG[which((diff_CpG_dbGSE141256$db)>=0.15)] #  17419
diff_CpG_db_hyperGSE141256<-diff_CpG_dbGSE141256$CpG[which((diff_CpG_dbGSE141256$db)<=(-0.15))] #  7667


#' ### load original organoid passage CpGs
load(here("data","beta_organoids.RData"))
pvals_long<-read.csv(here("data","Heteroskedactiy_pvalues_FDR_1000iter.csv"), header=T)
pvals_long[,1]<-NULL
colnames(pvals_long)<-c("BP_count","diff_count","mean_db","fdr_BP","fdr_diff")
pvals_long$CpG<-rownames(organoid_beta)
pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001

hetero_CpG<-rownames(organoid_beta)[which(pvals_long$BP_fdr<0.05)]
print(paste("CpGs with significant (adjusted p<0.05) heteroskedactiy: ", length(hetero_CpG), sep=""))
diff_CpG<-rownames(organoid_beta)[which(pvals_long$diff_fdr<0.05)]
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),]
print(paste("CpGs with significant (adjusted p<0.05; delta beta >0.05) differential methylation: ", nrow(diff_CpG_db), sep=""))

diff_CpG_db_hypo<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)>=0.15)] #  11772
diff_CpG_db_hyper<-diff_CpG_db$CpG[which((diff_CpG_db$mean_db)<=(-0.15))] #  5214


#' How many differential CpGs could overlap? 450K vs EPIC
print(paste("Of the ", nrow(diff_CpG_db), " CpGs differential with passage in the original organoids, ",
            length(diff_CpG_db$CpG[which(diff_CpG_db$CpG%in%rownames(combat_GSE141256_Beta_spheroids))]), " are in the GSE141256 data", sep=""))

diff_CpG_db_hypo_overlap<-diff_CpG_db_hypo[which(diff_CpG_db_hypo%in%rownames(combat_GSE141256_Beta_spheroids))]
diff_CpG_db_hyper_overlap<-diff_CpG_db_hyper[which(diff_CpG_db_hyper%in%rownames(combat_GSE141256_Beta_spheroids))]

diff_CpG_db_hypoGSE141256_overlap<-diff_CpG_db_hypoGSE141256[which(diff_CpG_db_hypoGSE141256%in%pvals_long$CpG)]
diff_CpG_db_hyperGSE141256_overlap<-diff_CpG_db_hyperGSE141256[which(diff_CpG_db_hyperGSE141256%in%pvals_long$CpG)]

print(paste("Of the ",length(diff_CpG_db_hypo_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hypoGSE141256_overlap, diff_CpG_db_hypo_overlap))," are also hypo in the GSE141256 cohort (",
            round((length(intersect(diff_CpG_db_hypoGSE141256_overlap, diff_CpG_db_hypo_overlap))/length(diff_CpG_db_hypo_overlap))*100,2),"%)",sep=""))

print(paste("Of the ",length(diff_CpG_db_hyper_overlap)," hypo CpGs also on the 450K ",
            length(intersect(diff_CpG_db_hyperGSE141256_overlap, diff_CpG_db_hyper_overlap))," are also hypo in the GSE141256 cohort (",
            round((length(intersect(diff_CpG_db_hyperGSE141256_overlap, diff_CpG_db_hyper_overlap))/length(diff_CpG_db_hyper_overlap))*100,2),"%)",sep=""))


#' ### delta beta directionality plot
plt_db_direction<-merge(pvals_long[,c(3,6)], passage_GSE141256, by="CpG")# 380776

plt_db_direction$sig<-"Not Significant"
plt_db_direction$sig[which(plt_db_direction$CpG%in%c(intersect(diff_CpG_db_hypoGSE141256_overlap, diff_CpG_db_hypo_overlap),intersect(diff_CpG_db_hyperGSE141256_overlap, diff_CpG_db_hyper_overlap)))]<-"Significant\nIn Both\nCohorts"

ggplot(plt_db_direction, aes(mean_db, db))+geom_point(aes(color=sig, alpha=sig),shape=19)+th+theme_bw()+
  scale_color_manual(values=c("lightgrey", "cornflowerblue"), name="Significant\nWith Passage")+
  scale_alpha_manual(values=c(0.25,1), guide=F)+
  geom_hline(yintercept=c(-0.15,0.15), color="grey60")+geom_vline(xintercept=c(-0.15,0.15), color="grey60")+
  ylim(-0.8,0.8)+xlim(-0.8,0.8)+xlab("Original Organoid\nPassage Delta Beta")+ylab("GSE141256 Organoid\nPassage Delta Beta")+
  stat_smooth(method="lm", se=F, color="black")


#ggsave(here("figs","GSE141256_db_directionality.pdf"), width=5, height=3.75)
ggsave(here("figs/jpeg","GSE141256_db_directionality.jpeg"), width=5, height=3.75)

print(paste("Correlation of delta betas between cohorts: ", round(cor(plt_db_direction$db, plt_db_direction$mean_db),2), sep=""))


#' ### representative CpGs
epic.organoid_minimal<-epic.organoid[,c(2, 14, 17)]
colnames(epic.organoid_minimal)[1]<-"Assay.Name"
epic.organoid_minimal$cohort<-"Original Organoids"

GSE141256_meta_combo_spheroids_minimal<-GSE141256_meta_combo_spheroids[,c(6,13,16)]
GSE141256_meta_combo_spheroids_minimal$cohort<-"MTAB-4957 Organoids"
colnames(GSE141256_meta_combo_spheroids_minimal)[2:3]<-c("sample_ID","passage.or.rescope.no_numeric")

sample_info_both<-rbind(GSE141256_meta_combo_spheroids_minimal,epic.organoid_minimal)


plt_hetero_GSE141256<-function(CpGs, legend, axislab, title){
  betas<-melt(cbind(combat_GSE141256_Beta_spheroids[CpGs,],organoid_beta[CpGs,]))
  organoid_plt<-merge(sample_info_both, betas, by.x="Assay.Name",by.y="Var2")
  
  p<-ggplot(organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_manual(values=pass_col,name="Passage\nNumber", drop=T)+
    facet_grid(cohort~Var1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))+
    xlim(1,16)
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}}

plt_hetero_GSE141256(c("cg25402228","cg22009751"))

ggsave(here("figs","Passage_differential_CpGs_GSE141256.pdf"),width = 4.75, height = 4)
ggsave(here("figs/jpeg","Passage_differential_CpGs_GSE141256.jpeg"), width = 4.75, height = 4)



#'## R Session Info
sessionInfo()
