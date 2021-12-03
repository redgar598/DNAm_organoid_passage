#'---
#'title: Normalization and quality control of DNAm arrays
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(rafalib)
library(scales)
library(grid)
library(here)





options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("scripts","00_pretty_plots.R"))
suppressMessages(source(here("scripts","00_heat_scree_plot_generic.R")))


#' ### Load Data
# Load the cohort meta data Felicity Payne shared
path<-"data/original_organoids"
epic.organoid<-read.csv(here(path, "epic_organoid_MethylationArraySamples_31Jul19.csv"))



#' ### Normalize DNAm Arrays
here(path)
epic.organoid$array.id.path <- file.path(here(path), epic.organoid$array.id)
# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$sex)
organoid_beta<-getBeta(MSet.illumina)

print(paste("Samples available: ",ncol(organoid_beta),"; Probes available: ",nrow(organoid_beta),sep=""))

#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
epic.organoid$det_pval<-avg_detPval

ggplot(epic.organoid)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix_ID), det_pval, group=sample_ID, fill=as.factor(sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("figs","detection_pvalue_organoids.pdf"), width=6, height=5)
ggsave(here("figs/jpeg","detection_pvalue_organoids.jpeg"), width=6, height=5)


#' Generally any sample with a average detection p value > 0.05 are removed, but none were in this cohort.
print(paste("Number of samples: ", nrow(epic.organoid), sep=""))
epic.organoid<-epic.organoid[which(epic.organoid$det_pval<0.005),]

#'Number of samples after removal of high detection p value: `r nrow(epic.organoid)`
#'
#' Normalize raw again but this time without bad samples
# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$sex)
organoid_beta<-getBeta(MSet.illumina)

print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' Beta distribution before and after normalization

# extract raw beta values for plotting
beta_raw<-getBeta(rgset_organoid)
identical(colnames(beta_raw),epic.organoid$array.id)

Beta_melted<- melt(organoid_beta)
Beta_melted_raw<- melt(beta_raw)

# Remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

# Add meta data
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,epic.organoid, by.x="ID", by.y="array.id")

ggplot(Beta_Plot_raw, aes(Beta, group=as.character(ID), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("figs","beta_distribution_raw_EPIC_organoid.pdf"), w=10, h=5)
ggsave(here("figs/jpeg","beta_distribution_raw_EPIC_organoid.jpeg"), w=10, h=5)


ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("figs","beta_distribution_EPIC_organoid.pdf"), w=10, h=5)
ggsave(here("figs/jpeg","beta_distribution_EPIC_organoid.jpeg"), w=10, h=5)





#'#### Confirm individuals ID with SNPs probes and clustering by DNAm

# remove rows with NAs
Betas_cluster<-organoid_beta[complete.cases(organoid_beta),]

d <- dist(t(Betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sample.site), cex=1.5)

pdf(here("figs","cluster_wholeEPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sample.site), cex=1.5)
dev.off()


#' #### Genotyping Probes
SNPs <- getSnpBeta(rgset_organoid)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA

SNPs<-SNPs[,which(colnames(SNPs)%in%epic.organoid$array.id)]
identical(colnames(SNPs),epic.organoid$array.id)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(as.character(epic.organoid$case.no)), cex=1.5)

pdf(here("figs","cluster_snps_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(as.character(epic.organoid$case.no)), cex=1.5)
dev.off()



#' #### Sex clustering
#' Using the cg ID to chromosome annotation from illumina 
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%anno_EPIC$IlmnID),]
anno_EPIC<-anno_EPIC[match(rownames(organoid_beta),anno_EPIC$IlmnID),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)

organoid_beta_sex<-organoid_beta[which(anno_EPIC$CHR%in%c('X','Y')),]

d <- dist(t(organoid_beta_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sex), cex=1.5)

pdf(here("figs","cluster_sex_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sex), cex=1.5)
dev.off()



#' #### Remove samples which do not cluster correctly
#' Samples 287_SC and T036_SC do not cluster as expected based on meta data.  
#' 287_SC is likely a tissue mislabel the wrong tissue as it is labeled SC and clusters with TI. 
#' T036_SC (passage 10 sample) clusters with the wrong sex as labelled and clustering on the SNPs probes it does not cluster with T036_SC

epic.organoid<-epic.organoid[which(!(epic.organoid$array.id%in%c("203548970031_R03C01","203548970036_R03C01"))),]
organoid_beta<-organoid_beta[,which(!(colnames(organoid_beta)%in%c("203548970031_R03C01","203548970036_R03C01")))]
identical(colnames(organoid_beta), as.character(epic.organoid$array.id))

print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' ### Probe Filtering 
# SNP probes should already be removed
organoid_beta <- organoid_beta[!grepl("rs",rownames(organoid_beta)), ]
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))

#' #### Sex Chromosomes 
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(organoid_beta),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)
organoid_beta <- organoid_beta[!anno_EPIC$CHR%in%c("X", "Y"), ]

filt_sex<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes. 
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#' "43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv
cross_reactive<-read.csv(here("data", "13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%cross_reactive$PROBE)),]

filt_cross<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#'For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites." and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes" but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("data", "13059_2016_1066_MOESM4_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(polymorphic$PROBE))," polymorphic probes (genetic variants at targeted CpG sites).", sep=""))

baseext<-read.csv(here("data", "13059_2016_1066_MOESM5_ESM.csv"), stringsAsFactors = F)
print(paste("Filtering ",length(unique(baseext$PROBE))," polymorphic probes (single base extension sites for Infinium Type I probes).", sep=""))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%c(polymorphic$PROBE, baseext$PROBE))),]

filt_poly<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe filtering based on detection pvalue and detection over background (NA)

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(organoid_beta), function(y) length(which(is.na(organoid_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(organoid_beta)*0.05))
organoid_beta<-organoid_beta[na_count_probe_good,]

filt_bead<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_organoid)
detP<-detP[,which(colnames(detP)%in%epic.organoid$array.id)]
identical(colnames(detP),epic.organoid$array.id)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%bad_det_p)),]
organoid_beta<-organoid_beta[,which(!(colnames(organoid_beta)%in%bad_det_psamp))]
identical(colnames(organoid_beta), as.character(epic.organoid$array.id))

filt_detp<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe attrition plot
df<-data.frame(sample_num_remaining=c(866238,865918,865859,filt_sex,filt_cross,filt_poly,filt_bead,filt_detp),
               filter=c("EPIC Probe Number","Missing Annotation Data","Removal of SNP Probes",
                        "Removal of X and Y chromosome probes","Removal of Cross Reactive Probes",
                        "Removal of Polymorphic Probes", "Removal of Probes with Beadcount <3\nin 5 % of Samples",
                        "Removal of Probes with 1 % of samples\nwith a detection p-value greater than 0.05"))
df$sample_num_lost<-c(0,sapply(2:nrow(df), function(x) df$sample_num_remaining[x-1]-df$sample_num_remaining[x]))

df$filter<-factor(df$filter, rev(df$filter))

ggplot(df)+
  geom_bar(aes(filter,-sample_num_remaining), stat="identity", fill="grey70", color="black")+
  geom_bar(aes(filter,sample_num_lost), stat="identity",fill="darkred", color="black")+
  geom_text(aes(x=filter, y=-min(sample_num_remaining)/2,  label=comma(sample_num_remaining)))+
  geom_text(aes(x=filter, y=max(sample_num_lost)/1.5,  label=comma(sample_num_lost)))+
  geom_hline(yintercept=0)+
  coord_flip()+theme_bw()+ylab("")+xlab("")+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "grey20", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(position = "top")

ggsave(here("figs","probe_attrition_EPIC_organoid.pdf"), width = 8, height = 3)
ggsave(here("figs/jpeg","probe_attrition_EPIC_organoid.jpeg"), width = 8, height = 3)



#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

# Restructure meta for PCA
epic.organoid$sentrix_ID<-as.factor(epic.organoid$sentrix_ID)
epic.organoid$case.no<-as.factor(epic.organoid$case.no)
epic.organoid$sex<-as.factor(epic.organoid$sex)
epic.organoid$sample.site<-as.factor(epic.organoid$sample.site)

epic.organoid$passage.or.rescope.no_numeric<-as.factor(as.character(epic.organoid$passage.or.rescope.no))
levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,11,14,16,2,3,4,6,7,8,2)
epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))

#' #### Cohort distribution
meta_categorical <- epic.organoid[, c(1,5,9,11)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid[, c(10, 17)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.", "Sample Site","Sex","Sentrix ID")
colnames(meta_continuous) <- c("Age", "Passage")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how many PCs to display in plot?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.4, 2.75))

ggsave(here("figs","heat_scree_EPIC_organoid.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.4, 2.75)),width = 9, height = 6)
ggsave(here("figs/jpeg","heat_scree_EPIC_organoid.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.4, 2.75)),width = 9, height = 6)



#' #### heat scree simplified
meta_categorical <- epic.organoid[, c(5,9)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid[, c(10, 17)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Sample Site","Gender")
colnames(meta_continuous) <- c("Age", "Passage")

ord<-c(1,4,3,2)
# how many PCs to display in plot?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.4, 2.75))

ggsave(here("figs","heat_scree_EPIC_organoid_simplified.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 2.8)),width = 9, height = 5)
ggsave(here("figs/jpeg","heat_scree_EPIC_organoid_simplified.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 2.8)),width = 9, height = 5)



#' ### PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, epic.organoid, by="array.id")

#' #### Sample Site
ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+fillscale_sampsite+xlab(paste("PC1 (",round(Importance[1]*100,0),"%)", sep=""))+ylab(paste("PC2 (",round(Importance[2]*100,0),"%)", sep=""))+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))

ggsave(here("figs","PC1_PC2_organoid.pdf"), width = 6.5, height = 5)
ggsave(here("figs/jpeg","PC1_PC2_organoid.jpeg"), width = 6.5, height = 5)


#' #### Passage
pc_plt<-ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_line(aes(PC2,PC3, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(shape=21,size=3, color="black")+#
  theme_bw()+xlab(paste("PC2 (",round(Importance[2]*100,0),"%)", sep=""))+ylab(paste("PC3 (",round(Importance[3]*100,0),"%)", sep=""))+th+theme(axis.text = element_text(size=12),axis.title = element_text(size=14))+
  scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")#+ #scale_color_manual(values=c("white", "black"))

legend<-ggplot(epic.organoid, aes(as.factor(-passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+
  coord_flip()+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+ylab("Sample\nNumber")+th

r <- ggplot() + theme_void()

grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(7,1))

ggsave(here("figs","PC2_PC3_organoid.pdf"), grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.45,1.25,0.19)), ncol=2, widths=c(6,1)),width = 7, height = 5)
ggsave(here("figs/jpeg","PC2_PC3_organoid.jpeg"), grid.arrange(pc_plt,arrangeGrob(r,legend,r, heights=c(0.6,1.25,0.4)), ncol=2, widths=c(8,1)),width = 8, height = 6)





#' ### Correlation of PC2 and passage
cor_coef<-suppressWarnings(cor.test(Loadings_meta$PC2, Loadings_meta$passage.or.rescope.no_numeric, method="spearman"))

ggplot(Loadings_meta, aes(passage.or.rescope.no_numeric, PC2))+geom_smooth(method="lm", se=F, color="black")+
  geom_line(aes(passage.or.rescope.no_numeric, PC2, group=sample_ID), color="lightgrey")+#, color=sampling.time.point
  geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=3)+th+theme_bw()+
  scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+annotate("text", x=12, y=25, label = paste("r = ",round(cor_coef$estimate, 3), sep=""), size=5)

ggsave(here("figs","PC2_passage_organoid.pdf"), width = 7.5, height = 6)
ggsave(here("figs/jpeg","PC2_passage_organoid.jpeg"), width = 7.5, height = 6)



#' ### Save intermediate object for further analysis
save(organoid_beta, epic.organoid, file=paste(here("data"),"/beta_organoids.RData",sep=""))
write.csv(organoid_beta, here("data", "beta_organoids.csv"))
write.csv(epic.organoid, here("data", "meta_organoids.csv"))



#'## R Session Info
sessionInfo()
