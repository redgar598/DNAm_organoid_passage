#'---
#'title: Adjust from cell composistion
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#load sleuth library
suppressMessages({
  library(sleuth)
  library(testit)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(here)
})
options(stringsAsFactors = FALSE)
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))



###########################
### Cell type deconvolution
###########################
load(here("rna_seq/data","ibd_RNAseq_sleuth.RData"))
sampleinfo_rnaseq<-sampleinfo_noOutlier
sampleinfo_rnaseq$sample_ID<-paste(sampleinfo_rnaseq$case.no, sampleinfo_rnaseq$sample.site)

load(file=here("DNAm/data","DNAm_organoid_counts.RData"))
sampleinfo_DNAm<-sampleinfo_DNAm_countsConstrained

sampleinfo_DNAm<-sampleinfo_DNAm[which(is.na(sampleinfo_DNAm$passage.or.rescope.no)),]
sampleinfo_DNAm$sample_ID<-paste(sampleinfo_DNAm$case.no, sampleinfo_DNAm$sample.site)


sampleinfo_rnaseq<-merge(sampleinfo_rnaseq, sampleinfo_DNAm[,c("sample_ID","organoid", "tcell")], by="sample_ID")

sampleinfo_rnaseq<-sampleinfo_rnaseq[match(colnames(mat_noOutlier), sampleinfo_rnaseq$sample),]
identical(colnames(mat_noOutlier), sampleinfo_rnaseq$sample)

## color PCA by count
pca_res_noOutlier <- prcomp(t(mat_noOutlier))
vars <- pca_res_noOutlier$sdev^2
Importance<-vars/sum(vars)

Loadings_meta_noOutlier<-as.data.frame(pca_res_noOutlier$x)
Loadings_meta_noOutlier$sample<-rownames(Loadings_meta_noOutlier)
Loadings_meta_noOutlier<-merge(Loadings_meta_noOutlier, sampleinfo_rnaseq, by="sample")
ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=organoid))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+
  scale_fill_gradientn(colours = c("#0b4899","cornflowerblue","white"), values=c(0,0.5,1),
                       name="Organoid Estimate")+xlab("PC1 (38%)")+ylab("PC2 (12%)")

ggsave(here("rna_seq/figs","PC1_PC2_RNAseq_organoid.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq_organoid.jpeg"), width = 7.5, height = 6)



#########################
### Cell composition adjustment
#########################
# Impute NAs
identical(colnames(mat_noOutlier), sampleinfo_rnaseq$sample)

#imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
#Blood_beta_imputed<-t(imputeMedianv3(Blood_beta))
avebeta.lm<-apply(mat_noOutlier, 1, function(x){
  lm(x~sampleinfo_rnaseq$organoid)
})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(mat_noOutlier)
ibd_rnaseq_adjusted<-residuals+matrix(apply(mat_noOutlier, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

save(ibd_rnaseq_adjusted,sampleinfo_rnaseq, file=here("rna_seq/data","ibd_adjusted_rnaseq.RData"))


## PCA after adjustement
pca_res_noOutlier <- prcomp(t(ibd_rnaseq_adjusted))
vars <- pca_res_noOutlier$sdev^2
Importance<-vars/sum(vars)

Loadings_meta_noOutlier<-as.data.frame(pca_res_noOutlier$x)
Loadings_meta_noOutlier$sample<-rownames(Loadings_meta_noOutlier)
Loadings_meta_noOutlier<-merge(Loadings_meta_noOutlier, sampleinfo_rnaseq, by="sample")
ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=organoid))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+
  scale_fill_gradientn(colours = c("#0b4899","cornflowerblue","white"), values=c(0,0.5,1),
                       name="Organoid Estimate")+xlab("PC1 (75%)")+ylab("PC2 (9%)")
ggsave(here("rna_seq/figs","PC1_PC2_RNAseq_organoid_adjusted.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq_organoid_adjusted.jpeg"), width = 7.5, height = 6)

ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+fillscale_sampsite+th+xlab("PC1 (75%)")+ylab("PC2 (9%)")
ggsave(here("rna_seq/figs","PC1_PC2_RNAseq_adjusted.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq_adjusted.jpeg"), width = 7.5, height = 6)

ggplot(Loadings_meta_noOutlier, aes(PC1, PC2, fill=diagnosis, color=sample.site))+geom_point(shape=21,size=3)+theme_bw()+
  fillscale_diagnosis+th+scale_color_manual(values=c("grey40","black","white"), name="Sample Site")+xlab("PC1 (75%)")+ylab("PC2 (9%)")
ggsave(here("rna_seq/figs","PC1_PC2_RNAseq_diagnosis_adjusted.pdf"), width = 7.5, height = 6)
ggsave(here("rna_seq/figs/jpeg","PC1_PC2_RNAseq_diagnosis_adjusted.jpeg"), width = 7.5, height = 6)

# ## moar PCs
# round(Importance, digits=2)[1:10]
# PC23<-ggplot(Loadings_meta_noOutlier, aes(PC2, PC3, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
#   fillscale_sampsite+th+xlab("PC2 (12%)")+ylab("PC3 (9%)")
# PC34<-ggplot(Loadings_meta_noOutlier, aes(PC3, PC4, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
#   fillscale_sampsite+th+xlab("PC3 (9%)")+ylab("PC4 (7%)")
# 
# ggsave("../../figs/PC234_RNAseq_adjusted.pdf", grid.arrange(PC23, PC34), width = 7.5, height = 12)
# ggsave("../../figs/jpeg/PC234_RNAseq_adjusted.jpeg", grid.arrange(PC23, PC34), width = 7.5, height = 12)
# 



# heat scree
sampleinfo_rnaseq$case.no<-as.factor(sampleinfo_rnaseq$case.no)
Loadings_noOutlier<-as.data.frame(pca_res_noOutlier$x)

meta_categorical <- sampleinfo_rnaseq[, c(2,6,7,12,13)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo_rnaseq[, c(14)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.", "Diagnosis", "Sample Site","Inflammation","Sex")
colnames(meta_continuous) <- c("Age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings_noOutlier, Importance, 3.3, 1.5))

ggsave(here("rna_seq/figs","heat_scree_RNAseq_noOutlier_adjusted.pdf"), suppressWarnings(heat_scree_plot(Loadings_noOutlier, Importance, 3.3, 1.5)),width = 9, height = 6)
ggsave(here("rna_seq/figs/jpeg","heat_scree_RNAseq_noOutlier_adjusted.jpeg"), suppressWarnings(heat_scree_plot(Loadings_noOutlier, Importance, 3.3, 1.5)),width = 9, height = 6)



#'## R Session Info
sessionInfo()

